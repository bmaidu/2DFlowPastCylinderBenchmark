"""Configuration and field checks; no CFD or mesh generation is required."""
from dataclasses import asdict
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "main"))

from mesh.generate_mesh import centerline_knots, refinement_bands, wake_size_expression
from mesh.mesh_parameters import Domain, MESH_OVERRIDE_FIELDS, MeshParameters, SIZE_FIELDS, VERTICAL_FIELDS
from run.case_config import CONFIG_PATH, resolve, validate
from run.run_case import SIZE_FIELDS as CLI_SIZE_FIELDS, parser

CROSS_OVERRIDES = {
    "domain_x": [-25.0, 45.0], "domain_y": [-25.0, 25.0],
    "wall_size": 0.0055, "near_wall_distance": 1.0, "near_wall_outer_size": 0.08,
    "far_upstream_size": 0.028, "upstream_size": 0.028, "approach_size": 0.028,
    "initial_wake_size": 0.028, "middle_wake_size": 0.028, "extended_wake_size": 0.028,
    "transition_size": 0.028, "outlet_size": 0.028,
    "wake_half_width": 3.0, "wake_blend_width": 1.0, "far_size": 0.175,
    "vertical_size": 0.055, "vertical_half_width": 3.0, "vertical_blend_width": 1.0,
    "distance_sampling": 960,
}


def field_value(parameters, x, y):
    return eval(wake_size_expression(parameters),
                {"__builtins__": {}, "min": min, "max": max, "abs": abs, "x": x, "y": y})


class ConfigurationTests(unittest.TestCase):
    def test_legacy_presets_load_with_vertical_band_disabled(self):
        base = json.loads(CONFIG_PATH.read_text())
        for preset, mesh in base["presets"].items():
            with self.subTest(preset=preset):
                config = resolve(preset)
                self.assertIsNone(config["mesh"]["vertical_size"])
                for key, value in mesh.items():
                    self.assertEqual(config["mesh"][key], value)
                validate(config)

    def test_original_size_field_list_is_preserved(self):
        mesh = json.loads(CONFIG_PATH.read_text())["presets"]["coarse"]
        self.assertEqual(set(SIZE_FIELDS), set(mesh) - {"algorithm", "distance_sampling"})
        self.assertEqual(len(MESH_OVERRIDE_FIELDS), len(set(MESH_OVERRIDE_FIELDS)))

    def test_original_cli_field_export_is_preserved(self):
        self.assertEqual(CLI_SIZE_FIELDS, SIZE_FIELDS)

    def test_old_saved_configuration_validates_without_rewriting(self):
        config = resolve("coarse")
        for key in VERTICAL_FIELDS:
            del config["mesh"][key]
        before = json.dumps(config, sort_keys=True)
        validate(config)
        self.assertEqual(json.dumps(config, sort_keys=True), before)

    def test_cli_help_handles_disabled_defaults(self):
        help_text = parser().format_help()
        self.assertIn("--vertical-size", help_text)
        self.assertIn("coarse=off", help_text)

    def test_cli_and_json_resolve_to_identical_configuration(self):
        arguments = []
        for name, value in CROSS_OVERRIDES.items():
            arguments.append("--" + name.replace("_", "-"))
            arguments.extend(str(item) for item in (value if isinstance(value, list) else [value]))
        parsed = parser().parse_args(arguments + ["--prepare-only"])
        overrides = {key: getattr(parsed, key) for key in ("domain_x", "domain_y", *MESH_OVERRIDE_FIELDS)}
        config = resolve(overrides=overrides)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "case_config.json"
            path.write_text(json.dumps(config))
            self.assertEqual(resolve(config_file=path), config)
        self.assertTrue(parsed.prepare_only)

    def test_json_null_disables_vertical_band(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "case_config.json"
            path.write_text(json.dumps({"mesh": {"vertical_size": None}}))
            self.assertIsNone(resolve(config_file=path)["mesh"]["vertical_size"])

    def test_cli_overrides_custom_json_settings(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "case_config.json"
            path.write_text(json.dumps({"mesh": {"vertical_size": 0.1}}))
            config = resolve(config_file=path, overrides={"vertical_size": 0.055})
            self.assertEqual(config["mesh"]["vertical_size"], 0.055)

    def test_unknown_mesh_setting_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "case_config.json"
            path.write_text(json.dumps({"mesh": {"vertical_szie": 0.055}}))
            with self.assertRaisesRegex(ValueError, "Unknown mesh"):
                resolve(config_file=path)

    def test_invalid_vertical_parameters_are_rejected(self):
        for name in VERTICAL_FIELDS:
            for value in (0, -1, float("inf"), float("nan"), True, "0.055"):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    MeshParameters(**{name: value})
        for name in ("vertical_half_width", "vertical_blend_width"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                MeshParameters(**{name: None})

    def test_vertical_size_respects_global_maximum(self):
        with self.assertRaisesRegex(ValueError, "far_size"):
            MeshParameters(vertical_size=0.8)

    def test_postprocess_arguments_remain_compatible(self):
        for arguments in (["--case-name", "existing", "--postprocess-only"],
                          ["--case-dir", "/tmp/existing", "--postprocess-only"]):
            self.assertTrue(parser().parse_args(arguments).postprocess_only)


class RefinementFieldTests(unittest.TestCase):
    def setUp(self):
        self.config = resolve(overrides=CROSS_OVERRIDES)
        self.parameters = MeshParameters(**self.config["mesh"])
        self.domain = Domain(-25, 45, -25, 25)

    def test_cross_targets_in_each_region(self):
        for point, expected in (((0, 10), 0.055), ((10, 0), 0.028),
                                ((10, 10), 0.175), ((0, 0), 0.028)):
            with self.subTest(point=point):
                self.assertAlmostEqual(field_value(self.parameters, *point), expected)

    def test_band_transitions_are_smooth_and_bounded(self):
        self.assertAlmostEqual(field_value(self.parameters, 3.5, 10), (0.055 + 0.175) / 2)
        self.assertAlmostEqual(field_value(self.parameters, 10, 3.5), (0.028 + 0.175) / 2)
        for coordinate in (3, 3.25, 3.5, 3.75, 4):
            self.assertGreaterEqual(field_value(self.parameters, coordinate, 10), 0.055)
            self.assertLessEqual(field_value(self.parameters, coordinate, 10), 0.175)
        epsilon = 1e-6
        for edge in (3, 4):
            self.assertAlmostEqual(field_value(self.parameters, edge - epsilon, 10),
                                   field_value(self.parameters, edge + epsilon, 10), places=10)

    def test_cross_field_is_symmetric(self):
        for x, y in ((3.5, 10), (10, 3.5), (2, 2)):
            expected = field_value(self.parameters, x, y)
            self.assertEqual(field_value(self.parameters, -x, y), expected)
            self.assertEqual(field_value(self.parameters, x, -y), expected)

    def test_graded_horizontal_wake_remains_available_with_vertical_band(self):
        parameters = MeshParameters(vertical_size=0.055)
        self.assertAlmostEqual(field_value(parameters, 10, 0), 0.038)
        self.assertAlmostEqual(field_value(parameters, 0, 10), 0.055)
        self.assertAlmostEqual(field_value(parameters, 10, 10), 0.7)

    def test_uniform_cross_metadata_records_all_refinement(self):
        bands = refinement_bands(self.parameters, self.domain)
        self.assertEqual(bands["horizontal"]["target_size"], 0.028)
        self.assertEqual(bands["horizontal"]["x"], [-25, 45])
        self.assertEqual(bands["vertical"], {
            "x": [-3.0, 3.0], "y": [-25, 25], "target_size": 0.055, "blend_width_x": 1.0})
        self.assertEqual(bands["far_field_target_size"], 0.175)

    def test_graded_metadata_does_not_claim_uniform_target_size(self):
        parameters = MeshParameters()
        bands = refinement_bands(parameters, self.domain)
        self.assertNotIn("target_size", bands["horizontal"])
        self.assertEqual(bands["horizontal"]["centerline_knots"], centerline_knots(parameters))
        self.assertNotIn("vertical", bands)

    def test_nondefault_vertical_widths_control_actual_field(self):
        parameters = MeshParameters(**{**asdict(self.parameters),
                                      "vertical_half_width": 2, "vertical_blend_width": 2})
        self.assertAlmostEqual(field_value(parameters, 2, 10), 0.055)
        self.assertAlmostEqual(field_value(parameters, 3, 10), 0.115)
        self.assertAlmostEqual(field_value(parameters, 4, 10), 0.175)


if __name__ == "__main__":
    unittest.main()
