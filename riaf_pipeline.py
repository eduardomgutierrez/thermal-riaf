#!/usr/bin/env python3
"""Build and run a thermal RIAF spectrum from one TOML configuration."""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(command, *, cwd=None, env=None):
    print("+", " ".join(map(str, command)))
    subprocess.run(command, cwd=cwd, check=True, env=env)


def nested_set(data, dotted_key, value):
    node = data
    parts = dotted_key.split(".")
    for part in parts[:-1]:
        node = node.setdefault(part, {})
    node[parts[-1]] = value


def validate(cfg):
    source = cfg.get("profile", {}).get("source", "hydro")
    if source not in {"hydro", "external"}:
        raise ValueError("profile.source must be 'hydro' or 'external'")
    for name in ("synchrotron", "bremsstrahlung", "comptonization"):
        if not isinstance(cfg.get("radiation", {}).get(name, True), bool):
            raise ValueError(f"radiation.{name} must be true or false")


def prepare_radproc_config(cfg, run_dir, hydro_dir):
    template = ROOT / "radproc/src/adaf/parameters.json"
    with template.open(encoding="utf-8") as stream:
        result = json.load(stream)
    radiation = cfg.get("radiation", {})
    flags = [
        int(radiation.get("synchrotron", True)),
        int(radiation.get("bremsstrahlung", True)),
        0, 0,
        int(radiation.get("comptonization", True)),
    ]
    result["calculateThermal"] = "1"
    result["calculateNonThermal"] = "0"
    result["calculateComptonScatt"] = str(int(radiation.get("calculate_scattering_matrix", True)))
    result["thermal"]["processes"] = {
        "synchrotron": flags[0],
        "bremsstrahlung": flags[1],
        "proton_proton": 0,
        "cold_disk": 0,
        "comptonization": flags[4],
    }
    for index, enabled in enumerate(flags):
        result["thermal"]["processNumber"][str(index)] = enabled
    grid = cfg.get("grid", {})
    overrides = {
        "model.particle.default.dim.radius.samples": str(grid.get("radius_samples", 50)),
        "model.particle.photon.dim.energy.samples": str(grid.get("energy_samples", 100)),
        "model.particle.photon.dim.energy.min": float(grid.get("log10_energy_ev_min", -6.0)),
        "model.particle.photon.dim.energy.max": float(grid.get("log10_energy_ev_max", 8.0)),
    }
    for key, value in overrides.items():
        nested_set(result, key, value)

    # The thermal Compton implementation uses precomputed probability grids.
    # Stage them explicitly so a run never depends on files left in the source
    # or build directory by an earlier calculation.
    table_dir = ROOT / "radproc/src/adaf"
    for table_name in ("probs.bin", "nt.bin", "om.bin", "omp.bin"):
        table = table_dir / table_name
        if not table.is_file():
            raise FileNotFoundError(f"required Compton table not found: {table}")
        shutil.copy2(table, run_dir / table_name)

    if cfg.get("profile", {}).get("source", "hydro") == "external":
        external = Path(cfg["profile"]["file"]).expanduser().resolve()
        if not external.is_file():
            raise FileNotFoundError(external)
        shutil.copy2(external, run_dir / external.name)
        result["readPrecomputedADAF"] = "1"
        result["solSpinFile"] = external.name
        result["blackHoleMass"] = str(cfg["profile"]["black_hole_mass_msun"])
        result["spin"] = str(cfg["profile"].get("spin", 0.0))
        result["R_high"] = str(cfg["profile"].get("ion_electron_temperature_ratio", 1.0))
    else:
        result["readPrecomputedADAF"] = "0"
        shutil.copy2(hydro_dir / "adafFile.txt", run_dir / "adafFile.txt")
        shutil.copy2(hydro_dir / "adafParameters.txt", run_dir / "adafParameters.txt")
    with (run_dir / "parameters.json").open("w", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--hydro-only", action="store_true")
    parser.add_argument("--no-build", action="store_true")
    args = parser.parse_args()
    with args.config.open("rb") as stream:
        cfg = tomllib.load(stream)
    validate(cfg)
    output = (ROOT / cfg.get("run", {}).get("output_dir", "runs/example")).resolve()
    hydro_dir, run_dir = output / "hydro", output / "spectrum"
    hydro_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)

    if cfg.get("profile", {}).get("source", "hydro") == "hydro":
        hydro_cfg = hydro_dir / "hydro-input.json"
        hydro_cfg.write_text(json.dumps(cfg.get("hydro", {}), indent=2) + "\n", encoding="utf-8")
        run([sys.executable, ROOT / "RIAF/solve.py", "--config", hydro_cfg,
             "--output-dir", hydro_dir])
    prepare_radproc_config(cfg, run_dir, hydro_dir)
    if args.hydro_only:
        print(f"prepared run in {output}")
        return
    build_dir = ROOT / "radproc/build"
    executable = build_dir / "src/adaf/adaf"
    if not args.no_build:
        if not shutil.which("meson"):
            raise RuntimeError("Meson is not installed; see README.md")
        # A failed Meson setup can leave the directory behind without a valid
        # build definition, so test Meson's core data rather than the directory.
        if not (build_dir / "meson-private/coredata.dat").is_file():
            run(["meson", "setup", build_dir, ROOT / "radproc"])
        run(["meson", "compile", "-C", build_dir])
    if not executable.is_file():
        raise FileNotFoundError(f"radproc executable not found: {executable}")
    runtime_env = os.environ.copy()
    runtime_env["OMP_NUM_THREADS"] = str(cfg.get("run", {}).get("omp_threads", 1))
    run([executable], cwd=run_dir, env=runtime_env)
    plot_path = run_dir / "thermal-diagnostics.png"
    run([sys.executable, ROOT / "plot_diagnostics.py", args.config.resolve(),
         "--spectrum-dir", run_dir, "--output", plot_path])
    print(f"thermal spectrum: {run_dir / 'lumThermal.dat'}")
    print(f"diagnostic plot: {plot_path}")


if __name__ == "__main__":
    try:
        main()
    except (KeyError, ValueError, RuntimeError, FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
