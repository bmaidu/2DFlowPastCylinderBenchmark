"""Plot the full domain, cylinder close-up, and near/middle wake mesh."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

FIGURE_NAMES = (
    "full_domain.png", "near_cylinder.png", "near_middle_wake.png",
)

# Close-up windows use coordinates normalized by the cylinder diameter.
VIEWS = (
    ("near_cylinder.png", "Near cylinder", (-2., 2.), (-2., 2.), .22),
    ("near_middle_wake.png", "Near and middle wake", (-2., 12.), (-4., 4.), .12),
)


def load_mesh(mesh_dir):
    path = Path(mesh_dir) / "mesh-complete.mesh.vtu"
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    if reader.GetErrorCode() or grid.GetNumberOfPoints() == 0:
        raise ValueError(f"Cannot read mesh: {path}")
    return {"coords": vtk_to_numpy(grid.GetPoints().GetData()),
            "tris": vtk_to_numpy(grid.GetCells().GetConnectivityArray()).reshape(-1, 3)}


def plot_mesh_figures(mesh_dir, mesh=None):
    """Render one case; full-domain limits follow its actual mesh bounds."""
    mesh_dir = Path(mesh_dir)
    mesh = load_mesh(mesh_dir) if mesh is None else mesh
    points = np.asarray(mesh["coords"])[:, :2]
    triangles = np.asarray(mesh["tris"], dtype=np.int64)
    output = mesh_dir / "figures"
    output.mkdir(parents=True, exist_ok=True)
    edges = np.vstack((triangles[:, [0,1]], triangles[:, [1,2]], triangles[:, [2,0]]))
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    segments = points[edges]
    lower, upper = segments.min(axis=1), segments.max(axis=1)
    x_bounds = (float(points[:, 0].min()), float(points[:, 0].max()))
    y_bounds = (float(points[:, 1].min()), float(points[:, 1].max()))
    views = (("full_domain.png", "Full domain", x_bounds, y_bounds, .06), *VIEWS)
    for filename, title, x_span, y_span, linewidth in views:
        figure, axis = plt.subplots(figsize=(8.2, 5.8), constrained_layout=True)
        x_view = (max(x_span[0], x_bounds[0]), min(x_span[1], x_bounds[1]))
        y_view = (max(y_span[0], y_bounds[0]), min(y_span[1], y_bounds[1]))
        visible = (upper[:, 0] >= x_span[0]) & (lower[:, 0] <= x_span[1])
        visible &= (upper[:, 1] >= y_span[0]) & (lower[:, 1] <= y_span[1])
        axis.add_collection(LineCollection(segments[visible], colors="0.15", linewidths=linewidth))
        axis.set_xlim(*x_view)
        axis.set_ylim(*y_view)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlabel("x/D")
        axis.set_ylabel("y/D")
        axis.set_title(f"{title}\n{len(points):,} nodes; {len(triangles):,} triangles")
        figure.savefig(output / filename, dpi=180)
        plt.close(figure)

    return [str(Path("figures") / name) for name in FIGURE_NAMES]
