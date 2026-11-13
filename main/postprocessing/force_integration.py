"""Direct TRI3 cylinder forces and velocity-gradient calculations.

The fluid force on the cylinder, per unit span, is

    F' = -∮_Γ sigma · n_f ds

where n_f is the outward normal of the fluid domain. On the cylinder hole
that normal points from the fluid into the cylinder. With

    sigma = -p I + mu (grad u + (grad u)^T)

the pressure part is +∮ p n_f ds and the viscous part is -∮ tau · n_f ds.

Linear TRI3 velocity has a constant gradient inside each triangle, so the
viscous stress on a boundary edge is taken from that one parent triangle.
Nodal linear pressure is integrated exactly by the edge average. Both use the saved linear finite-element fields directly.
"""

from __future__ import annotations

import numpy as np


def triangle_geometry(points, triangles):
    """Signed areas and constant TRI3 shape-function gradients.

    points: (N, 2), triangles: (M, 3) integer connectivity.
    grad_n[e, a, j] is dN_a / dx_j on element e. The formulas use the signed
    area, so a clockwise triangle still produces the correct derivatives.
    """
    triangles = np.asarray(triangles, dtype=np.int64)
    coords = np.asarray(points, dtype=float)[triangles]
    edge_1 = coords[:, 1] - coords[:, 0]
    edge_2 = coords[:, 2] - coords[:, 0]
    twice_area = edge_1[:, 0] * edge_2[:, 1] - edge_1[:, 1] * edge_2[:, 0]
    x = coords[:, :, 0]
    y = coords[:, :, 1]
    dndx = np.column_stack(
        (y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1])
    ) / twice_area[:, None]
    dndy = np.column_stack(
        (x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0])
    ) / twice_area[:, None]
    gradients = np.stack((dndx, dndy), axis=-1)
    return 0.5 * twice_area, gradients


def velocity_gradients(velocity, triangles, gradients):
    """Constant velocity gradient on every triangle.

    gradient[e, i, j] = du_i / dx_j.
    """
    nodal = np.asarray(velocity, dtype=float)[np.asarray(triangles)]
    return np.einsum("eai,eaj->eij", nodal, gradients)


def parent_elements(triangles, edges):
    """Map each unique boundary edge to its single adjacent triangle."""
    triangles = np.asarray(triangles, dtype=np.int64)
    edges = np.asarray(edges, dtype=np.int64)
    owners = {}
    for element, nodes in enumerate(triangles):
        pairs = ((nodes[0], nodes[1]), (nodes[1], nodes[2]), (nodes[2], nodes[0]))
        for left, right in pairs:
            key = (int(left), int(right)) if left < right else (int(right), int(left))
            if key in owners:
                owners[key] = -1
            else:
                owners[key] = element
    parents = np.empty(len(edges), dtype=np.int64)
    missing = []
    for index, (left, right) in enumerate(edges):
        key = (int(left), int(right)) if left < right else (int(right), int(left))
        parent = owners.get(key, None)
        if parent is None or parent < 0:
            missing.append(key)
            parents[index] = -1
        else:
            parents[index] = parent
    if missing:
        raise ValueError(
            f"{len(missing)} boundary edges do not have a unique fluid parent."
        )
    return parents


def outward_unit_normals(points, triangles, edges, parents):
    """Fluid-domain outward unit normals, independent of stored edge order.

    The normal is rotated from the edge tangent and then flipped so it points
    from the parent-triangle centroid toward the edge. On an exterior wall
    that direction leaves the fluid. On the cylinder hole it enters the hole.
    """
    points = np.asarray(points, dtype=float)
    edges = np.asarray(edges, dtype=np.int64)
    parents = np.asarray(parents, dtype=np.int64)
    centroids = points[np.asarray(triangles)[parents]].mean(axis=1)
    ends = points[edges]
    midpoints = 0.5 * (ends[:, 0] + ends[:, 1])
    tangents = ends[:, 1] - ends[:, 0]
    normals = np.column_stack((tangents[:, 1], -tangents[:, 0]))
    inward = np.einsum("ij,ij->i", midpoints - centroids, normals) < 0.0
    normals[inward] *= -1.0
    lengths = np.linalg.norm(tangents, axis=1)
    if np.any(lengths <= 0.0):
        raise ValueError("A boundary edge has zero length.")
    return normals / lengths[:, None], lengths


def boundary_geometry(points, triangles, edges, gradients):
    parents = parent_elements(triangles, edges)
    normals, lengths = outward_unit_normals(points, triangles, edges, parents)
    return {
        "edges": np.asarray(edges, dtype=np.int64),
        "parents": parents,
        "normals": normals,
        "lengths": lengths,
        "gradients": gradients[parents],
    }


def _edge_velocity_gradient(velocity, geometry):
    nodal = np.asarray(velocity, dtype=float)[geometry["triangle_nodes"]]
    return np.einsum("eai,eaj->eij", nodal, geometry["gradients"])


def pressure_force(pressure, geometry):
    """Exact integral of a nodal-linear pressure against a straight edge."""
    edges = geometry["edges"]
    average = 0.5 * (
        np.asarray(pressure, dtype=float)[edges[:, 0]]
        + np.asarray(pressure, dtype=float)[edges[:, 1]]
    )
    return np.sum(
        average[:, None] * geometry["normals"] * geometry["lengths"][:, None],
        axis=0,
    )


def viscous_stress(gradients, mu):
    return float(mu) * (gradients + np.swapaxes(gradients, -1, -2))


def viscous_force_from_gradients(gradients, geometry, mu):
    traction = np.einsum("eij,ej->ei", viscous_stress(gradients, mu), geometry["normals"])
    return -np.sum(traction * geometry["lengths"][:, None], axis=0)


def direct_viscous_force(velocity, geometry, mu):
    """Viscous force using the constant gradient of each parent triangle."""
    return viscous_force_from_gradients(
        _edge_velocity_gradient(velocity, geometry), geometry, mu
    )


def element_vorticity(element_gradients):
    """Out-of-plane vorticity duy/dx - dux/dy."""
    return element_gradients[:, 1, 0] - element_gradients[:, 0, 1]


def coefficient_scale(density, speed, diameter):
    scale = 0.5 * float(density) * float(speed) ** 2 * float(diameter)
    if scale == 0.0:
        raise ValueError("Force-coefficient scale is zero.")
    return scale


def closed_normal_sum(geometry):
    return np.sum(geometry["normals"] * geometry["lengths"][:, None], axis=0)


def uniform_pressure_force(geometry, pressure):
    return pressure_force(np.full(int(geometry["edges"].max()) + 1, pressure), geometry)


def prepare_geometry(points, triangles, boundary_edges):
    """Build volume and boundary operators used by every time step."""
    signed_area, gradients = triangle_geometry(points, triangles)
    if np.any(signed_area <= 0.0):
        raise ValueError("The volume mesh contains a clockwise or degenerate triangle.")
    boundaries = {}
    for name, edges in boundary_edges.items():
        geometry = boundary_geometry(points, triangles, edges, gradients)
        geometry["triangle_nodes"] = np.asarray(triangles, dtype=np.int64)[
            geometry["parents"]
        ]
        boundaries[name] = geometry
    return {
        "points": np.asarray(points, dtype=float),
        "triangles": np.asarray(triangles, dtype=np.int64),
        "signed_area": signed_area,
        "gradients": gradients,
        "boundaries": boundaries,
    }


def validate_force_operator(prepared, cylinder_name="cylinder", radius=0.5):
    """Analytic checks on the cylinder operator before any CFD value is trusted."""
    geometry = prepared["boundaries"][cylinder_name]
    length = float(np.sum(geometry["lengths"]))
    analytic = 2.0 * np.pi * float(radius)
    normal_sum = closed_normal_sum(geometry)
    gauge = uniform_pressure_force(geometry, 1.0)
    zero_velocity = np.zeros((len(prepared["points"]), 2))
    viscous = direct_viscous_force(zero_velocity, geometry, 0.01)
    order = np.arange(len(geometry["edges"]))
    rng = np.random.default_rng(7)
    rng.shuffle(order)
    shuffled = {
        key: (value[order] if isinstance(value, np.ndarray) and len(value) == len(order) else value)
        for key, value in geometry.items()
    }
    sample_pressure = prepared["points"][:, 0]
    original = pressure_force(sample_pressure, geometry)
    reordered = pressure_force(sample_pressure, shuffled)
    reversed_edges = geometry["edges"][:, ::-1]
    reversed_geometry = boundary_geometry(
        prepared["points"],
        prepared["triangles"],
        reversed_edges,
        prepared["gradients"],
    )
    reversed_geometry["triangle_nodes"] = prepared["triangles"][
        reversed_geometry["parents"]
    ]
    reversed_force = pressure_force(sample_pressure, reversed_geometry)
    checks = {
        "boundary_length": length,
        "analytic_circumference": analytic,
        "relative_length_error": abs(length - analytic) / analytic,
        "normal_sum": normal_sum.tolist(),
        "normal_sum_norm": float(np.linalg.norm(normal_sum)),
        "unit_pressure_force": gauge.tolist(),
        "unit_pressure_force_norm": float(np.linalg.norm(gauge)),
        "zero_gradient_viscous_force": viscous.tolist(),
        "reorder_force_difference": float(np.linalg.norm(original - reordered)),
        "reversed_edge_force_difference": float(np.linalg.norm(original - reversed_force)),
        "edge_count": int(len(geometry["edges"])),
        "unique_parents": int(len(np.unique(geometry["parents"]))),
    }
    tolerances = {
        "relative_length_error": 1.0e-3,
        "normal_sum_norm": 1.0e-8,
        "unit_pressure_force_norm": 1.0e-8,
        "zero_gradient_viscous_norm": 1.0e-12,
        "reorder_force_difference": 1.0e-10,
        "reversed_edge_force_difference": 1.0e-10,
    }
    failures = []
    if checks["relative_length_error"] > tolerances["relative_length_error"]:
        failures.append("boundary length")
    if checks["normal_sum_norm"] > tolerances["normal_sum_norm"]:
        failures.append("closed normal sum")
    if checks["unit_pressure_force_norm"] > tolerances["unit_pressure_force_norm"]:
        failures.append("uniform pressure cancellation")
    if float(np.linalg.norm(viscous)) > tolerances["zero_gradient_viscous_norm"]:
        failures.append("zero viscous force")
    if checks["reorder_force_difference"] > tolerances["reorder_force_difference"]:
        failures.append("edge reordering")
    if checks["reversed_edge_force_difference"] > tolerances["reversed_edge_force_difference"]:
        failures.append("reversed edge connectivity")
    if checks["edge_count"] != checks["unique_parents"]:
        failures.append("parent coverage")
    checks["tolerances"] = tolerances
    checks["failures"] = failures
    checks["passed"] = not failures
    return checks
