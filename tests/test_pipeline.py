import tempfile
import unittest
from pathlib import Path

import riaf_pipeline as pipeline


class PipelineConfigurationTests(unittest.TestCase):
    def test_default_profile_and_thermal_switches_are_valid(self):
        pipeline.validate({"radiation": {
            "synchrotron": True,
            "bremsstrahlung": True,
            "comptonization": True,
        }})

    def test_unknown_profile_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "profile.source"):
            pipeline.validate({"profile": {"source": "mystery"}})

    def test_nested_set(self):
        value = {}
        pipeline.nested_set(value, "a.b.c", 4)
        self.assertEqual(value, {"a": {"b": {"c": 4}}})

    def write_profile(self, directory, rows):
        path = Path(directory) / "profile.dat"
        path.write_text("# header\n" + "\n".join(rows) + "\n", encoding="utf-8")
        return path

    def test_external_profile_validation_and_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_profile(directory, [
                "10 1e-10 0.01 0 0 2 0.5 1 1e10 0 0.1",
                "2.1 2e-10 0.1 0 0 5 0.3 2 2e10 0 0.2",
            ])
            result = pipeline.validate_external_profile(path, 10.0, 0.0)
            self.assertEqual(result["rows"], 2)
            self.assertEqual(result["radius_outer_rg"], 10.0)
            self.assertGreater(result["outer_accretion_power_erg_s"], 0.0)

    def test_external_profile_rejects_wrong_column_count(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_profile(directory, ["10 1 0.1"])
            with self.assertRaisesRegex(ValueError, "expected 11 columns"):
                pipeline.validate_external_profile(path, 10.0, 0.0)

    def test_external_profile_rejects_non_decreasing_radius(self):
        with tempfile.TemporaryDirectory() as directory:
            row = "10 1e-10 0.01 0 0 2 0.5 1 1e10 0 0.1"
            path = self.write_profile(directory, [row, row])
            with self.assertRaisesRegex(ValueError, "strictly decreasing"):
                pipeline.validate_external_profile(path, 10.0, 0.0)


if __name__ == "__main__":
    unittest.main()
