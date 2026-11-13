#!/usr/bin/env python3
"""2D flow past cylinder benchmark mesh generation.

Mesh is generated using a rectangle minus a unit cylinder with Gmsh Frontal-Delaunay triangles.

Combine near-wall distance refinement with a graded wake size field, then
export the volume mesh and boundary edges as solver-ready vtu/vtp files.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict
import json
from pathlib import Path
import time

import gmsh
import numpy as np
import vtk
from vtk.util.numpy_support import (
    numpy_to_vtk,
    numpy_to_vtkIdTypeArray,
    vtk_to_numpy,
)

from .mesh_parameters import (
    BOUNDARY_NAMES, CENTER, MESH_ALGORITHM,
    MESH_ALGORITHM_NAME, RADIUS, Domain, MeshParameters,
)


# Target sizes are interpolated along x/D and blended into the far-field size.
SIZE_SCHEDULE = (
    ("far_upstream_size", -15.0),
    ("upstream_size", -5.0),
    ("approach_size", -2.0),
    ("initial_wake_size", 0.5),
    ("initial_wake_size", 5.0),
    ("middle_wake_size", 10.0),
    ("extended_wake_size", 15.0),
    ("transition_size", 20.0),
    ("outlet_size", 25.0),
)

GEOMETRY_TOLERANCE = 1e-5


def centerline_knots(parameters: MeshParameters) -> list[tuple[float, float]]:
    return [(x_value, float(getattr(parameters, name))) for name, x_value in SIZE_SCHEDULE]


def wake_size_expression(parameters: MeshParameters) -> str:
    """Interpolate streamwise targets and smoothly blend into the transverse far field."""
    knots = centerline_knots(parameters)
    along_x = f"{knots[0][1]:.8g}"
    for (x0, h0), (x1, h1) in zip(knots, knots[1:]):
        fraction = f"max(0,min(1,((x)-({x0:.8g}))/({(x1 - x0):.8g})))"
        along_x = f"({along_x}+({(h1 - h0):.8g})*({fraction}))"
    band = (
        f"max(0,min(1,(abs(y)-({parameters.wake_half_width:.8g}))"
        f"/({parameters.wake_blend_width:.8g})))"
    )
    smooth = f"(({band})*({band})*(3-2*({band})))"
    return f"(({along_x})*(1-({smooth}))+({parameters.far_size:.8g})*({smooth}))"


def _int32_array(values: np.ndarray, name: str):
    result = numpy_to_vtk(np.ascontiguousarray(values, dtype=np.int32), deep=True)
    result.SetName(name)
    return result


def _cell_array(connectivity: np.ndarray) -> vtk.vtkCellArray:
    connectivity = np.ascontiguousarray(connectivity, dtype=np.int64)
    count, width = connectivity.shape
    offsets = np.arange(0, (count + 1) * width, width, dtype=np.int64)
    cells = vtk.vtkCellArray()
    cells.SetData(
        numpy_to_vtkIdTypeArray(offsets, deep=True),
        numpy_to_vtkIdTypeArray(connectivity.ravel(), deep=True),
    )
    return cells


def _write_xml(writer, path: Path, data) -> None:
    writer.SetFileName(str(path))
    writer.SetInputData(data)
    writer.SetDataModeToAppended()
    writer.EncodeAppendedDataOff()
    writer.SetCompressorTypeToNone()
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write {path}")


def _classify_curve(tag: int, domain: Domain) -> str:
    """Identify a boundary from coordinates sampled along its geometric curve."""
    tmin, tmax = gmsh.model.getParametrizationBounds(1, tag)
    t0 = float(np.asarray(tmin, dtype=float).reshape(-1)[0])
    t1 = float(np.asarray(tmax, dtype=float).reshape(-1)[0])
    parameters = np.linspace(t0, t1, 24, endpoint=False)
    raw = gmsh.model.getValue(1, tag, np.ascontiguousarray(parameters, dtype=np.float64))
    samples = np.asarray(raw, dtype=np.float64).reshape(-1, 3)
    xy = samples[:, :2]
    radius = np.linalg.norm(xy, axis=1)
    if np.max(np.abs(radius - RADIUS)) < GEOMETRY_TOLERANCE:
        return "cylinder"
    tests = (
        ("inlet", np.abs(xy[:, 0] - domain.x_min)),
        ("outlet", np.abs(xy[:, 0] - domain.x_max)),
        ("top", np.abs(xy[:, 1] - domain.y_max)),
        ("bottom", np.abs(xy[:, 1] - domain.y_min)),
    )
    for name, deviation in tests:
        if np.max(deviation) < GEOMETRY_TOLERANCE:
            return name
    bbox = gmsh.model.getBoundingBox(1, tag)
    raise RuntimeError(f"Cannot classify OCC curve {tag}: bbox={bbox}")


def build_geometry(domain: Domain) -> tuple[int, dict[str, list[int]]]:
    """Subtract the cylinder and identify inlet, outlet, walls, and cylinder curves."""
    gmsh.model.add("cylinder_2d")
    rectangle = gmsh.model.occ.addRectangle(
        domain.x_min, domain.y_min, 0.0,
        domain.x_max - domain.x_min, domain.y_max - domain.y_min,
    )
    disk = gmsh.model.occ.addDisk(*CENTER, 0.0, RADIUS, RADIUS)
    fluid, _ = gmsh.model.occ.cut([(2, rectangle)], [(2, disk)])
    gmsh.model.occ.synchronize()
    if len(fluid) != 1 or fluid[0][0] != 2:
        raise RuntimeError(f"Expected one fluid surface, found {fluid}.")
    surface = fluid[0][1]
    groups = {name: [] for name in BOUNDARY_NAMES}
    for _, tag in gmsh.model.getBoundary([(2, surface)], oriented=False):
        groups[_classify_curve(tag, domain)].append(tag)
    if any(not tags for tags in groups.values()):
        raise RuntimeError(f"Incomplete boundary classification: {groups}")
    return surface, groups


def apply_size_field(parameters: MeshParameters, cylinder_curves: list[int]) -> dict:
    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "CurvesList", cylinder_curves)
    gmsh.model.mesh.field.setNumber(distance, "Sampling", parameters.distance_sampling)

    threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMin", parameters.wall_size)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMax", parameters.near_wall_outer_size)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", parameters.near_wall_distance)
    gmsh.model.mesh.field.setNumber(threshold, "StopAtDistMax", 1)

    expression = wake_size_expression(parameters)
    background = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(background, "F", expression)

    minimum = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold, background])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
    return {
        "wake_expression": expression,
        "centerline_knots": centerline_knots(parameters),
    }


def configure_meshing(parameters: MeshParameters) -> None:
    """Use deterministic, first-order triangular meshing with explicit size limits."""
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.option.setNumber("Mesh.MaxNumThreads1D", 1)
    gmsh.option.setNumber("Mesh.MaxNumThreads2D", 1)
    gmsh.option.setNumber("Mesh.Algorithm", MESH_ALGORITHM)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.RecombineAll", 0)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.8 * parameters.wall_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", parameters.far_size)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromParametricPoints", 0)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
    gmsh.option.setNumber("Mesh.RandomFactor", 1e-9)


def _edge_keys(edges: np.ndarray) -> np.ndarray:
    left = np.minimum(edges[:, 0], edges[:, 1]).astype(np.uint64)
    right = np.maximum(edges[:, 0], edges[:, 1]).astype(np.uint64)
    return (left << np.uint64(32)) | right


def _signed_areas(coords: np.ndarray, tris: np.ndarray) -> np.ndarray:
    xy = coords[:, :2]
    p0 = xy[tris[:, 0]]
    p1 = xy[tris[:, 1]]
    p2 = xy[tris[:, 2]]
    return 0.5 * ((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))


def _node_indices(node_tags: np.ndarray, element_nodes: np.ndarray) -> np.ndarray:
    indices = np.searchsorted(node_tags, element_nodes)
    if np.any(indices >= len(node_tags)) or not np.array_equal(node_tags[indices], element_nodes):
        raise RuntimeError("An element references a node outside the fluid mesh.")
    return indices


def extract_mesh(surface: int, groups: dict[str, list[int]], domain: Domain) -> dict:
    node_tags, raw_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.asarray(node_tags, dtype=np.int64)
    coords = np.asarray(raw_coords, dtype=np.float64).reshape(-1, 3)
    order = np.argsort(node_tags)
    node_tags = node_tags[order]
    coords = coords[order]

    element_types, _, element_nodes = gmsh.model.mesh.getElements(2, surface)
    triangles = []
    for element_type, nodelist in zip(element_types, element_nodes):
        if element_type != 2:
            raise RuntimeError(f"Volume mesh contains non-TRI3 element type {element_type}")
        triangles.append(np.asarray(nodelist, dtype=np.int64).reshape(-1, 3))
    if not triangles:
        raise RuntimeError("Gmsh did not produce any triangles")
    tri_tags = np.vstack(triangles)
    tris = _node_indices(node_tags, tri_tags).astype(np.int64)

    used = np.zeros(len(node_tags), dtype=bool)
    used[tris.ravel()] = True
    if not np.all(used):
        compact = np.cumsum(used) - 1
        tris = compact[tris]
        coords = coords[used]
        node_tags = node_tags[used]
    if np.max(np.abs(coords[:, 2])) > 1e-8:
        raise RuntimeError(f"Mesh nodes leave the z=0 plane: max |z|={np.max(np.abs(coords[:, 2]))}")
    coords[:, 2] = 0.0
    bounds = [coords[:, 0].min(), coords[:, 0].max(), coords[:, 1].min(), coords[:, 1].max()]
    expected = [domain.x_min, domain.x_max, domain.y_min, domain.y_max]
    if not np.allclose(bounds, expected, atol=1e-8, rtol=0):
        raise RuntimeError(f"Mesh bounds {bounds} do not match {expected}.")

    areas = _signed_areas(coords, tris)
    if np.any(areas == 0.0) or np.any(~np.isfinite(areas)):
        raise RuntimeError("Mesh contains a degenerate triangle")
    clockwise = areas < 0.0
    if np.any(clockwise):
        tris = tris.copy()
        tris[clockwise] = tris[clockwise][:, [0, 2, 1]]
        areas = np.abs(areas)

    boundary = {}
    for name, curves in groups.items():
        pieces = []
        for curve in curves:
            curve_types, _, curve_nodes = gmsh.model.mesh.getElements(1, curve)
            for element_type, nodelist in zip(curve_types, curve_nodes):
                if element_type != 1:
                    raise RuntimeError(f"{name} contains non-LINE2 element type {element_type}")
                pieces.append(np.asarray(nodelist, dtype=np.int64).reshape(-1, 2))
        if not pieces:
            raise RuntimeError(f"Boundary {name} has no line elements")
        raw = np.vstack(pieces)
        boundary[name] = _node_indices(node_tags, raw).astype(np.int64)

    topology = _topology(tris)
    exterior = topology["exterior"]
    exterior_keys = _edge_keys(exterior)
    key_order = np.argsort(exterior_keys)
    sorted_exterior_keys = exterior_keys[key_order]
    grouped = []
    parents = {}
    for name in BOUNDARY_NAMES:
        keys = _edge_keys(boundary[name])
        if len(np.unique(keys)) != len(keys):
            raise RuntimeError(f"Boundary {name} contains a duplicate edge")
        positions = np.searchsorted(sorted_exterior_keys, keys)
        if np.any(positions >= len(sorted_exterior_keys)) or not np.array_equal(
            sorted_exterior_keys[positions], keys
        ):
            raise RuntimeError(f"Boundary {name} contains an edge that is not an exterior triangle edge")
        matched_parents = topology["exterior_parents"][key_order[positions]]
        boundary[name] = _orient_outward(coords, boundary[name], matched_parents, tris)
        parents[name] = matched_parents
        grouped.append(keys)
    all_keys = np.concatenate(grouped)
    if len(all_keys) != len(exterior_keys) or len(np.unique(all_keys)) != len(all_keys):
        raise RuntimeError("Boundary groups do not partition the exterior edges")

    area = float(areas.sum())
    relative_area_error = (area - domain.fluid_area) / domain.fluid_area
    if abs(relative_area_error) > 1e-4:
        raise RuntimeError(f"Triangulated area error {relative_area_error:.3e} exceeds 1e-4")
    if _component_count(len(tris), topology["neighbors"]) != 1:
        raise RuntimeError("Fluid mesh has more than one connected component")
    _require_single_cylinder_loop(coords, boundary["cylinder"])

    return {
        "coords": coords,
        "tris": tris.astype(np.int64),
        "boundary": boundary,
        "parents": parents,
        "area": area,
        "relative_area_error": float(relative_area_error),
        "n_exterior": int(len(exterior)),
    }


def _topology(tris: np.ndarray) -> dict:
    edges = np.vstack((tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]])).astype(np.int64)
    owners = np.tile(np.arange(len(tris), dtype=np.int64), 3)
    keys = _edge_keys(edges)
    order = np.argsort(keys, kind="mergesort")
    keys = keys[order]
    edges = edges[order]
    owners = owners[order]
    change = np.flatnonzero(np.diff(keys)) + 1
    starts = np.concatenate(([0], change))
    counts = np.diff(np.concatenate((starts, [len(keys)])))
    if np.any((counts != 1) & (counts != 2)):
        raise RuntimeError("Mesh edge incidence is not 1 or 2")
    interior = counts == 2
    exterior_mask = counts == 1
    neighbors = np.column_stack((owners[starts[interior]], owners[starts[interior] + 1]))
    exterior_edges = np.sort(edges[starts[exterior_mask]], axis=1)
    exterior_parents = owners[starts[exterior_mask]]
    return {
        "exterior": exterior_edges,
        "exterior_parents": exterior_parents,
        "neighbors": neighbors,
    }


def _orient_outward(coords, lines, parents, tris) -> np.ndarray:
    """Orient boundary-edge normals away from their adjacent fluid triangles."""
    oriented = lines.copy()
    p0 = coords[oriented[:, 0], :2]
    p1 = coords[oriented[:, 1], :2]
    tangent = p1 - p0
    normal = np.column_stack((tangent[:, 1], -tangent[:, 0]))
    midpoint = 0.5 * (p0 + p1)
    centroid = coords[tris[parents], :2].mean(axis=1)
    dots = np.einsum("ij,ij->i", normal, midpoint - centroid)
    flip = dots < 0.0
    oriented[flip] = oriented[flip, ::-1]
    return oriented


def _component_count(n_elem: int, neighbors: np.ndarray) -> int:
    parent = np.arange(n_elem)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for left, right in neighbors:
        root_left = find(int(left))
        root_right = find(int(right))
        if root_left != root_right:
            parent[root_right] = root_left
    return len({find(index) for index in range(n_elem)})


def _require_single_cylinder_loop(coords: np.ndarray, lines: np.ndarray) -> None:
    """A closed, connected cylinder boundary is required for the force integral."""
    adjacency = defaultdict(set)
    for left, right in lines:
        adjacency[int(left)].add(int(right))
        adjacency[int(right)].add(int(left))
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise RuntimeError("Cylinder boundary is not a closed loop.")
    visited = set()
    pending = [next(iter(adjacency))]
    while pending:
        node = pending.pop()
        if node not in visited:
            visited.add(node)
            pending.extend(adjacency[node] - visited)
    if len(visited) != len(adjacency):
        raise RuntimeError("Cylinder boundary has disconnected loops.")
    radius_error = np.abs(np.linalg.norm(coords[list(visited), :2], axis=1) - RADIUS)
    if radius_error.max() > 1e-8:
        raise RuntimeError(f"Cylinder node leaves r={RADIUS} by {radius_error.max()}.")


def write_outputs(output_dir: Path, mesh: dict) -> None:
    mesh_dir = output_dir / "mesh"
    surface_dir = mesh_dir / "mesh-surfaces"
    surface_dir.mkdir(parents=True, exist_ok=True)
    coords = mesh["coords"]
    tris = mesh["tris"]
    keys = [f"{x:.16e},{y:.16e}," for x, y in coords[:, :2]]
    if len(set(keys)) != len(keys):
        raise RuntimeError("Distinct nodes collide in the svMultiPhysics x,y coordinate hash")

    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np.ascontiguousarray(coords, dtype=np.float64), deep=True))
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.SetCells(vtk.VTK_TRIANGLE, _cell_array(tris))
    grid.GetPointData().AddArray(_int32_array(np.arange(1, len(coords) + 1, dtype=np.int32), "GlobalNodeID"))
    grid.GetCellData().AddArray(_int32_array(np.arange(1, len(tris) + 1, dtype=np.int32), "GlobalElementID"))
    _write_xml(vtk.vtkXMLUnstructuredGridWriter(), mesh_dir / "mesh-complete.mesh.vtu", grid)

    for name in BOUNDARY_NAMES:
        lines = mesh["boundary"][name]
        parents = mesh["parents"][name]
        used, inverse = np.unique(lines, return_inverse=True)
        local = inverse.reshape(-1, 2)
        poly = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        points.SetData(numpy_to_vtk(np.ascontiguousarray(coords[used], dtype=np.float64), deep=True))
        poly.SetPoints(points)
        poly.SetLines(_cell_array(local))
        poly.GetPointData().AddArray(_int32_array(used + 1, "GlobalNodeID"))
        poly.GetCellData().AddArray(_int32_array(parents + 1, "GlobalElementID"))
        _write_xml(vtk.vtkXMLPolyDataWriter(), surface_dir / f"{name}.vtp", poly)


def verify_written_files(output_dir: Path, mesh: dict) -> None:
    """Reload the VTU and VTP files and compare them with the in-memory mesh."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(output_dir / "mesh" / "mesh-complete.mesh.vtu"))
    reader.Update()
    grid = reader.GetOutput()
    types = vtk_to_numpy(grid.GetCellTypes())
    if grid.GetNumberOfPoints() != len(mesh["coords"]) or grid.GetNumberOfCells() != len(mesh["tris"]):
        raise RuntimeError("Reloaded volume mesh has different dimensions")
    if not np.all(types == vtk.VTK_TRIANGLE):
        raise RuntimeError("Reloaded volume mesh is not exclusively VTK_TRIANGLE")
    coords = vtk_to_numpy(grid.GetPoints().GetData())
    connectivity = vtk_to_numpy(grid.GetCells().GetConnectivityArray()).reshape(-1, 3)
    node_ids = vtk_to_numpy(grid.GetPointData().GetArray("GlobalNodeID"))
    element_ids = vtk_to_numpy(grid.GetCellData().GetArray("GlobalElementID"))
    if node_ids.dtype != np.int32 or not np.array_equal(node_ids, np.arange(1, len(coords) + 1, dtype=np.int32)):
        raise RuntimeError("Volume GlobalNodeID is not contiguous Int32 1..N")
    if element_ids.dtype != np.int32 or not np.array_equal(
        element_ids, np.arange(1, len(mesh["tris"]) + 1, dtype=np.int32)
    ):
        raise RuntimeError("Volume GlobalElementID is not contiguous Int32 1..Ne")
    if not np.array_equal(coords, mesh["coords"]) or not np.array_equal(connectivity, mesh["tris"]):
        raise RuntimeError("Reloaded volume coordinates or connectivity changed")

    for name in BOUNDARY_NAMES:
        surface_reader = vtk.vtkXMLPolyDataReader()
        surface_reader.SetFileName(str(output_dir / "mesh" / "mesh-surfaces" / f"{name}.vtp"))
        surface_reader.Update()
        poly = surface_reader.GetOutput()
        if poly.GetNumberOfPolys() or poly.GetNumberOfStrips() or poly.GetNumberOfVerts():
            raise RuntimeError(f"{name}.vtp contains cells other than lines")
        lines = poly.GetLines()
        offsets = vtk_to_numpy(lines.GetOffsetsArray())
        if len(offsets) < 2 or not np.all(np.diff(offsets) == 2):
            raise RuntimeError(f"{name}.vtp is not exclusively LINE2")
        global_nodes = vtk_to_numpy(poly.GetPointData().GetArray("GlobalNodeID"))
        global_elements = vtk_to_numpy(poly.GetCellData().GetArray("GlobalElementID"))
        if global_nodes.dtype != np.int32 or global_elements.dtype != np.int32:
            raise RuntimeError(f"{name}.vtp ID arrays are not Int32")
        local = vtk_to_numpy(lines.GetConnectivityArray()).reshape(-1, 2)
        volume_lines = global_nodes[local].astype(np.int64) - 1
        parents = global_elements.astype(np.int64) - 1
        if not np.array_equal(volume_lines, mesh["boundary"][name]) or not np.array_equal(
            parents, mesh["parents"][name]
        ):
            raise RuntimeError(f"Reloaded {name}.vtp does not match the generated boundary")
        surface_xyz = vtk_to_numpy(poly.GetPoints().GetData())
        if not np.array_equal(surface_xyz, mesh["coords"][global_nodes.astype(np.int64) - 1]):
            raise RuntimeError(f"{name}.vtp coordinates do not match the volume nodes")


def generate(parameters: MeshParameters, output_dir: Path, domain: Domain = Domain()) -> dict:
    """Generate and verify mesh files; keep generation records alongside them."""
    output_dir = Path(output_dir).expanduser().resolve()
    mesh_dir = output_dir / "mesh"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    from .plot_mesh import plot_mesh_figures

    started = time.perf_counter()
    with (mesh_dir / "mesh_generation.log").open("w") as log:
        try:
            gmsh.initialize()
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.logger.start()
            try:
                surface, groups = build_geometry(domain)
                fields = apply_size_field(parameters, groups["cylinder"])
                configure_meshing(parameters)
                gmsh.model.mesh.generate(2)
                mesh = extract_mesh(surface, groups, domain)
            finally:
                messages = [message.replace(str(output_dir) + "/", "") for message in gmsh.logger.get()]
                log.write("\n".join(messages) + "\n")
                gmsh.logger.stop()
                gmsh.finalize()
            write_outputs(output_dir, mesh)
            verify_written_files(output_dir, mesh)
            figures = plot_mesh_figures(mesh_dir, mesh)
            summary = {
                "output_dir": ".",
                "domain": asdict(domain),
                "gmsh_version": gmsh.__version__,
                "figures": figures,
                "algorithm": {"id": MESH_ALGORITHM, "name": MESH_ALGORITHM_NAME},
                "parameters": asdict(parameters),
                **fields,
                "nodes": len(mesh["coords"]),
                "triangles": len(mesh["tris"]),
                "exterior_edges": mesh["n_exterior"],
                "boundary_edges": {name: len(edges) for name, edges in mesh["boundary"].items()},
                "area": mesh["area"],
                "analytical_area": domain.fluid_area,
                "relative_area_error": mesh["relative_area_error"],
                "gmsh_warnings_or_errors": [line for line in messages if "Warning" in line or "Error" in line],
                "elapsed_seconds": time.perf_counter() - started,
            }
            (mesh_dir / "mesh_generation_report.json").write_text(json.dumps(summary, indent=2) + "\n")
            log.write(f"Verified {summary['nodes']} nodes and {summary['triangles']} triangles.\n")
            return summary
        except Exception as error:
            reason = str(error).replace(str(output_dir), ".")
            log.write(f"Generation failed: {reason}\n")
            raise


if __name__ == "__main__":
    raise SystemExit("Use python main/run/run_case.py --prepare-only (with optional case/domain/size arguments).")
