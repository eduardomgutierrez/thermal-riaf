#!/usr/bin/env python3
"""Run a thermal RIAF spectrum from one TOML configuration."""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
C_LIGHT = 2.99792458e10
GRAVITATIONAL_CONSTANT = 6.67430e-8
SOLAR_MASS = 1.98847e33


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
    threads = cfg.get("run", {}).get("omp_threads", 4)
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("run.omp_threads must be a positive integer")
    seed = cfg.get("radiation", {}).get("scattering_random_seed", 5489)
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("radiation.scattering_random_seed must be a non-negative integer")


def validate_external_profile(path, black_hole_mass_msun, spin):
    """Validate an 11-column external profile and return useful metadata."""
    if black_hole_mass_msun <= 0.0:
        raise ValueError("profile.black_hole_mass_msun must be positive")
    if not -1.0 < spin < 1.0:
        raise ValueError("profile.spin must be strictly between -1 and 1")

    rows = []
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 11:
                raise ValueError(
                    f"{path}:{line_number}: expected 11 columns, found {len(fields)}")
            try:
                row = [float(field) for field in fields]
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: non-numeric value") from exc
            if not all(math.isfinite(value) for value in row):
                raise ValueError(f"{path}:{line_number}: values must be finite")
            radius, density, vr, _br, _bphi, bmag, h_over_r, beta, temp, _omega, vphi = row
            checks = (
                (radius > 0.0, "radius must be positive"),
                (density > 0.0, "density must be positive"),
                (0.0 < vr < 1.0, "radial velocity v/c must be between 0 and 1"),
                (bmag >= 0.0, "magnetic-field magnitude cannot be negative"),
                (h_over_r > 0.0, "H/r must be positive"),
                (beta > 0.0, "plasma beta must be positive"),
                (temp > 0.0, "ion temperature must be positive"),
                (abs(vphi) < 1.0, "azimuthal velocity magnitude must be below c"),
            )
            for valid, message in checks:
                if not valid:
                    raise ValueError(f"{path}:{line_number}: {message}")
            if rows and radius >= rows[-1][0]:
                raise ValueError(
                    f"{path}:{line_number}: radii must be strictly decreasing")
            rows.append(row)
    if len(rows) < 2:
        raise ValueError(f"{path}: expected at least two data rows")

    horizon_rg = 1.0 + math.sqrt(1.0 - spin * spin)
    if rows[-1][0] <= horizon_rg:
        raise ValueError(
            f"{path}: inner radius {rows[-1][0]:.8g} r_g is not outside "
            f"the Kerr horizon ({horizon_rg:.8g} r_g for spin {spin:g})")

    grav_radius = (GRAVITATIONAL_CONSTANT * black_hole_mass_msun * SOLAR_MASS
                   / C_LIGHT**2)

    def accretion_rate(row):
        radius_cm = row[0] * grav_radius
        return 4.0 * math.pi * radius_cm**2 * row[6] * row[1] * row[2] * C_LIGHT

    return {
        "source_file": str(path),
        "rows": len(rows),
        "radius_outer_rg": rows[0][0],
        "radius_inner_rg": rows[-1][0],
        "horizon_radius_rg": horizon_rg,
        "density_min_g_cm3": min(row[1] for row in rows),
        "density_max_g_cm3": max(row[1] for row in rows),
        "ion_temperature_min_k": min(row[8] for row in rows),
        "ion_temperature_max_k": max(row[8] for row in rows),
        "outer_accretion_rate_g_s": accretion_rate(rows[0]),
        "inner_accretion_rate_g_s": accretion_rate(rows[-1]),
        "outer_accretion_power_erg_s": accretion_rate(rows[0]) * C_LIGHT**2,
        "inner_accretion_power_erg_s": accretion_rate(rows[-1]) * C_LIGHT**2,
    }


def prepare_radproc_config(cfg, run_dir, hydro_dir, config_dir=Path.cwd()):
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
    result["scatt"]["randomSeed"] = str(radiation.get("scattering_random_seed", 5489))
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
        external = Path(cfg["profile"]["file"]).expanduser()
        if not external.is_absolute():
            external = config_dir / external
        external = external.resolve()
        if not external.is_file():
            raise FileNotFoundError(external)
        mass = float(cfg["profile"]["black_hole_mass_msun"])
        spin = float(cfg["profile"].get("spin", 0.0))
        r_high = float(cfg["profile"].get(
            "r_high", cfg["profile"].get("ion_electron_temperature_ratio", 1.0)))
        if r_high <= 0.0:
            raise ValueError("profile.r_high must be positive")
        summary = validate_external_profile(external, mass, spin)
        summary["r_high"] = r_high
        (run_dir / "profile-summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"external profile: {summary['rows']} rows, "
              f"{summary['radius_inner_rg']:.6g}--{summary['radius_outer_rg']:.6g} r_g")
        shutil.copy2(external, run_dir / external.name)
        result["readPrecomputedADAF"] = "1"
        result["solSpinFile"] = external.name
        result["blackHoleMass"] = str(mass)
        result["spin"] = str(spin)
        result["R_high"] = str(r_high)
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
    parser.add_argument("--build", action="store_true",
                        help="configure and compile radproc before running")
    parser.add_argument("--no-build", action="store_true", help=argparse.SUPPRESS)
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
    prepare_radproc_config(cfg, run_dir, hydro_dir, args.config.resolve().parent)
    if args.hydro_only:
        print(f"prepared run in {output}")
        return
    build_dir = ROOT / "radproc/build"
    executable = build_dir / "src/adaf/adaf"
    if args.build:
        if not shutil.which("cmake"):
            raise RuntimeError("CMake is not installed; see README.md")
        if not (build_dir / "CMakeCache.txt").is_file():
            run(["cmake", "-S", ROOT / "radproc", "-B", build_dir,
                 "-DCMAKE_BUILD_TYPE=Release"])
        run(["cmake", "--build", build_dir, "--parallel"])
    if not executable.is_file():
        raise FileNotFoundError(
            f"radproc executable not found: {executable}; run ./build.sh first "
            "or rerun with --build")
    runtime_env = os.environ.copy()
    runtime_env["OMP_NUM_THREADS"] = str(cfg.get("run", {}).get("omp_threads", 4))
    run([executable], cwd=run_dir, env=runtime_env)
    plot_path = run_dir / "thermal-diagnostics.png"
    plot_env = os.environ.copy()
    plot_env["MPLCONFIGDIR"] = str(run_dir / ".matplotlib")
    run([sys.executable, ROOT / "plot_diagnostics.py", args.config.resolve(),
         "--spectrum-dir", run_dir, "--output", plot_path], env=plot_env)
    print(f"thermal spectrum: {run_dir / 'lumThermal.dat'}")
    print(f"diagnostic plot: {plot_path}")


if __name__ == "__main__":
    try:
        main()
    except (KeyError, ValueError, RuntimeError, FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
