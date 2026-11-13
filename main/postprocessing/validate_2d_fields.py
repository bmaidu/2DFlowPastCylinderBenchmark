"""Read solver inputs, discover snapshots, and validate their mesh and fields."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from .force_integration import prepare_geometry


def _text(parent, tag):
    element = parent.find(tag)
    if element is None or element.text is None:
        raise ValueError(f"Solver XML is missing <{tag}>.")
    return " ".join(element.text.split())


def read_solver_config(path):
    root = ET.parse(path).getroot()
    general = root.find("GeneralSimulationParameters")
    equation = root.find("Add_equation")
    mesh = root.find("Add_mesh")
    viscosity = equation.find("Viscosity")
    faces = {
        face.get("name"): _text(face, "Face_file_path")
        for face in mesh.findall("Add_face")
    }
    boundaries = {}
    for bc in equation.findall("Add_BC"):
        boundaries[bc.get("name")] = {
            "type": _text(bc, "Type"),
            "value": float(_text(bc, "Value")) if bc.find("Value") is not None else None,
            "profile": _text(bc, "Profile") if bc.find("Profile") is not None else None,
            "effective_direction": (
                _text(bc, "Effective_direction")
                if bc.find("Effective_direction") is not None
                else None
            ),
        }
    return {
        "xml": Path(path),
        "dimensions": int(_text(general, "Number_of_spatial_dimensions")),
        "n_steps": int(_text(general, "Number_of_time_steps")),
        "dt": float(_text(general, "Time_step_size")),
        "vtk_increment": int(_text(general, "Increment_in_saving_VTK_files")),
        "save_start": int(_text(general, "Start_saving_after_time_step")),
        "results_dir": _text(general, "Save_results_in_folder"),
        "vtk_prefix": _text(general, "Name_prefix_of_saved_VTK_files"),
        "density": float(_text(equation, "Density")),
        "viscosity": float(_text(viscosity, "Value")),
        "mesh_path": _text(mesh, "Mesh_file_path"),
        "faces": faces,
        "boundaries": boundaries,
    }


def read_unstructured(path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    if reader.GetErrorCode():
        raise RuntimeError(f"VTK could not read {path}.")
    grid = reader.GetOutput()
    if grid is None or grid.GetNumberOfPoints() == 0:
        raise RuntimeError(f"{path} contains no points.")
    return grid


def read_polydata(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    if reader.GetErrorCode():
        raise RuntimeError(f"VTK could not read {path}.")
    return reader.GetOutput()


def point_coordinates(grid):
    return np.asarray(vtk_to_numpy(grid.GetPoints().GetData()), dtype=float)


def triangle_connectivity(grid):
    cells = grid.GetCells()
    offsets = vtk_to_numpy(cells.GetOffsetsArray())
    if len(offsets) < 2 or np.any(np.diff(offsets) != 3):
        raise ValueError("A result cell is not a 3-node triangle.")
    connectivity = vtk_to_numpy(cells.GetConnectivityArray()).reshape(-1, 3).astype(np.int64)
    types = vtk_to_numpy(grid.GetCellTypes())
    if not np.all(types == vtk.VTK_TRIANGLE):
        raise ValueError("A result cell is not VTK_TRIANGLE.")
    return connectivity


def array_as_matrix(array):
    values = np.asarray(vtk_to_numpy(array), dtype=float)
    if values.ndim == 1:
        return values.reshape(-1, 1)
    return values


def stored_time(grid):
    array = grid.GetFieldData().GetArray("TimeValue")
    if array is None:
        return None
    values = np.asarray(vtk_to_numpy(array), dtype=float).ravel()
    return None if len(values) == 0 else float(values[0])


def require_point_field(grid, name, width):
    array = grid.GetPointData().GetArray(name)
    if array is None:
        raise ValueError(f"Result is missing point array {name}.")
    values = array_as_matrix(array)
    if values.shape[1] < width:
        raise ValueError(f"{name} has {values.shape[1]} components, expected {width}.")
    if not np.all(np.isfinite(values[:, :width])):
        raise ValueError(f"{name} contains a nonfinite value.")
    return values[:, :width] if width > 1 else values[:, 0]


def read_snapshot(path, reference_points, reference_triangles):
    grid = read_unstructured(path)
    if grid.GetNumberOfPoints() != len(reference_points) or grid.GetNumberOfCells() != len(reference_triangles):
        raise ValueError(f"{path.name} does not match the volume mesh size.")
    points = point_coordinates(grid)
    if points.shape != reference_points.shape:
        raise ValueError(f"{path.name} point array does not match the volume mesh.")
    coordinate_error = float(np.max(np.abs(points - reference_points)))
    if coordinate_error > 1.0e-3:
        raise ValueError(f"{path.name} coordinates differ from the mesh by {coordinate_error:.3e}.")
    triangles = triangle_connectivity(grid)
    if not np.array_equal(triangles, reference_triangles):
        raise ValueError(f"{path.name} connectivity does not match the volume mesh.")
    velocity = require_point_field(grid, "Velocity", 2)
    pressure = require_point_field(grid, "Pressure", 1)
    if len(velocity) != len(reference_points) or len(pressure) != len(reference_points):
        raise ValueError(f"{path.name} field length does not match the mesh.")
    return velocity, pressure, stored_time(grid), coordinate_error


def discover_results(directory, prefix, dt, stable_seconds, max_step):
    directory = Path(directory)
    records = []
    rejected = []
    for path in directory.glob(f"{prefix}_*.vtu"):
        suffix = path.stem[len(prefix) + 1 :]
        if not suffix.isdigit():
            rejected.append({"file": str(path), "reason": "step suffix is not an integer"})
            continue
        step = int(suffix)
        if max_step is not None and step > int(max_step):
            continue
        age = path.stat().st_mtime
        records.append({"step": step, "path": path, "mtime": age, "bytes": path.stat().st_size})
    if not records:
        return [], rejected
    import time

    now = time.time()
    stable = []
    for item in records:
        if stable_seconds > 0 and now - item["mtime"] < stable_seconds:
            rejected.append(
                {
                    "file": str(item["path"]),
                    "step": item["step"],
                    "reason": f"newer than the {stable_seconds:g} s stability window",
                }
            )
            continue
        stable.append(item)
    stable.sort(key=lambda item: item["step"])
    steps = [item["step"] for item in stable]
    duplicates = sorted({step for step in steps if steps.count(step) > 1})
    if duplicates:
        raise ValueError(f"Duplicate result steps: {duplicates}")
    for item in stable:
        item["filename_time"] = item["step"] * float(dt)
    return stable, rejected


def expected_steps(config):
    start = config["vtk_increment"]
    if config["save_start"] >= start:
        start = config["vtk_increment"] * int(
            np.ceil(config["save_start"] / config["vtk_increment"])
        )
    return list(range(start, config["n_steps"] + 1, config["vtk_increment"]))


def line_edges(surface, volume_points):
    array = surface.GetPointData().GetArray("GlobalNodeID")
    if array is None:
        raise ValueError("Boundary VTP is missing GlobalNodeID.")
    global_ids = np.asarray(vtk_to_numpy(array), dtype=np.int64)
    if np.any(global_ids < 1) or np.any(global_ids > len(volume_points)):
        raise ValueError("Boundary GlobalNodeID is outside the volume node range.")
    local_to_volume = global_ids - 1
    surface_points = point_coordinates(surface)[:, :2]
    mismatch = np.max(np.abs(surface_points - volume_points[local_to_volume, :2]))
    edges = []
    for index in range(surface.GetNumberOfCells()):
        cell = surface.GetCell(index)
        if cell.GetCellType() != vtk.VTK_LINE or cell.GetNumberOfPoints() != 2:
            raise ValueError(f"Boundary cell {index} is not LINE2.")
        edges.append(
            [
                int(local_to_volume[cell.GetPointId(0)]),
                int(local_to_volume[cell.GetPointId(1)]),
            ]
        )
    return np.asarray(edges, dtype=np.int64), float(mismatch)


def load_case_geometry(project, config):
    mesh_path = project / config["mesh_path"]
    grid = read_unstructured(mesh_path)
    points = point_coordinates(grid)
    triangles = triangle_connectivity(grid)
    boundary_edges = {}
    mismatches = {}
    for name, relative in config["faces"].items():
        surface = read_polydata(project / relative)
        edges, mismatch = line_edges(surface, points)
        boundary_edges[name] = edges
        mismatches[name] = mismatch
    prepared = prepare_geometry(points[:, :2], triangles, boundary_edges)
    prepared["coordinate_z"] = points[:, 2]
    prepared["boundary_coordinate_mismatch"] = mismatches
    prepared["n_points"] = int(len(points))
    prepared["n_triangles"] = int(len(triangles))
    return prepared
