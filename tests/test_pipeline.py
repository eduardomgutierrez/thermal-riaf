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


if __name__ == "__main__":
    unittest.main()
