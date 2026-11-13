"""Geometry and sizing shared by case configuration and mesh generation."""
from dataclasses import dataclass, fields
import math

RADIUS = 0.5
CENTER = (0.0, 0.0)
MESH_ALGORITHM = 6
MESH_ALGORITHM_NAME = "Frontal-Delaunay"

BOUNDARY_NAMES = ("inlet", "outlet", "top", "bottom", "cylinder")


@dataclass(frozen=True)
class Domain:
    x_min: float = -15.0
    x_max: float = 25.0
    y_min: float = -12.0
    y_max: float = 12.0

    def __post_init__(self):
        # Dataclasses call this after construction to validate preset and custom bounds.
        values = (self.x_min, self.x_max, self.y_min, self.y_max)
        if not all(isinstance(value, (int, float)) and math.isfinite(value) for value in values):
            raise ValueError("Domain bounds must be finite numbers.")
        if not (self.x_min < -RADIUS and self.x_max > RADIUS
                and self.y_min < -RADIUS and self.y_max > RADIUS):
            raise ValueError("Domain must strictly contain the unit-diameter cylinder.")

    @property
    def fluid_area(self):
        return (self.x_max - self.x_min) * (self.y_max - self.y_min) - math.pi * RADIUS**2


@dataclass(frozen=True)
class MeshParameters:
    """Regional target sizes and refinement distances, in cylinder diameters."""
    wall_size: float = 0.0055
    near_wall_distance: float = 1.0
    near_wall_outer_size: float = 0.080
    far_upstream_size: float = 0.68
    upstream_size: float = 0.36
    approach_size: float = 0.080
    initial_wake_size: float = 0.028
    middle_wake_size: float = 0.038
    extended_wake_size: float = 0.055
    transition_size: float = 0.11
    outlet_size: float = 0.22
    wake_half_width: float = 3.0
    wake_blend_width: float = 1.0
    far_size: float = 0.70
    algorithm: int = MESH_ALGORITHM
    distance_sampling: int = 960

    def __post_init__(self):
        for name in SIZE_FIELDS:
            value = getattr(self, name)
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.wall_size > self.near_wall_outer_size:
            raise ValueError("wall_size must not exceed near_wall_outer_size.")
        if any(getattr(self, name) > self.far_size for name in SIZE_FIELDS if name.endswith("_size")):
            raise ValueError("Regional size targets must not exceed far_size (the global maximum).")
        if self.algorithm != MESH_ALGORITHM:
            raise ValueError("This benchmark uses Frontal-Delaunay (algorithm 6).")
        if not isinstance(self.distance_sampling, int) or self.distance_sampling < 2:
            raise ValueError("distance_sampling must be an integer >= 2.")


SIZE_FIELDS = tuple(field.name for field in fields(MeshParameters)
                    if field.name not in ("algorithm", "distance_sampling"))
