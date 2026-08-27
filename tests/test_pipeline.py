import tempfile
import unittest
import json
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

    def test_parallel_settings_are_validated(self):
        with self.assertRaisesRegex(ValueError, "omp_threads"):
            pipeline.validate({"run": {"omp_threads": 0}})
        with self.assertRaisesRegex(ValueError, "scattering_random_seed"):
            pipeline.validate({"radiation": {"scattering_random_seed": -1}})

    def test_nested_set(self):
        value = {}
        pipeline.nested_set(value, "a.b.c", 4)
        self.assertEqual(value, {"a": {"b": {"c": 4}}})

    def test_matching_compton_cache_prompts_for_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = Path(directory)
            (run_dir / "parameters.json").write_text(json.dumps({
                "calculateComptonScatt": "1",
                "readPrecomputedADAF": "0",
                "model": {"radius_samples": 2},
            }), encoding="utf-8")
            (run_dir / "adafFile.txt").write_text("profile", encoding="utf-8")
            (run_dir / "adafParameters.txt").write_text("parameters", encoding="utf-8")
            signature = pipeline.compton_cache_signature(run_dir)
            for name in pipeline.COMPTON_MATRIX_FILES:
                (run_dir / name).write_text("matrix", encoding="utf-8")
            pipeline.record_compton_cache(run_dir, signature)

            self.assertTrue(pipeline.ask_to_reuse_compton_cache(
                run_dir, signature, input_fn=lambda _prompt: "yes"))
            self.assertFalse(pipeline.ask_to_reuse_compton_cache(
                run_dir, signature, input_fn=lambda _prompt: "no"))

            (run_dir / "adafFile.txt").write_text("changed", encoding="utf-8")
            changed_signature = pipeline.compton_cache_signature(run_dir)
            self.assertFalse(pipeline.matching_compton_cache(run_dir, changed_signature))

    def test_matching_hydro_cache_prompts_for_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            hydro_dir = Path(directory)
            config = hydro_dir / "hydro-input.json"
            config.write_text('{"blackHoleMass": 10}', encoding="utf-8")
            signature = pipeline.hydro_cache_signature(config)
            for name in pipeline.HYDRO_OUTPUT_FILES:
                (hydro_dir / name).write_text("output", encoding="utf-8")
            pipeline.record_hydro_cache(hydro_dir, signature)

            self.assertTrue(pipeline.ask_to_reuse_hydro_cache(
                hydro_dir, signature, input_fn=lambda _prompt: ""))
            self.assertFalse(pipeline.ask_to_reuse_hydro_cache(
                hydro_dir, signature, input_fn=lambda _prompt: "n"))

            config.write_text('{"blackHoleMass": 11}', encoding="utf-8")
            changed_signature = pipeline.hydro_cache_signature(config)
            self.assertFalse(pipeline.matching_hydro_cache(hydro_dir, changed_signature))

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
