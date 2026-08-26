#!/usr/bin/env python3
"""Plot thermal RIAF spectra, radial emission, and energy-budget checks."""

import argparse
import json
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

C_LIGHT = 2.99792458e10
EDDINGTON_ACCRETION_RATE = 1.39e18  # g s^-1 per solar mass


def integrated_luminosity(frequency, nu_lnu):
    """Return integral L = integral (nu L_nu) d ln(nu)."""
    return float(np.trapezoid(nu_lnu, x=np.log(frequency)))


def positive(values, relative_floor=1.0e-12):
    """Mask numerical underflow tails on logarithmic plots."""
    maximum = np.nanmax(values)
    floor = maximum * relative_floor
    return np.where(values > floor, values, np.nan)


def energy_budget(cfg, spectrum, spectrum_dir=None):
    total = integrated_luminosity(spectrum[:, 0], spectrum[:, 8])
    if cfg.get("profile", {}).get("source", "hydro") != "hydro":
        summary_path = spectrum_dir / "profile-summary.json" if spectrum_dir else None
        if summary_path and summary_path.is_file():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            return (total, summary["outer_accretion_power_erg_s"],
                    summary["inner_accretion_power_erg_s"])
        return total, None, None
    hydro = cfg["hydro"]
    outer = (EDDINGTON_ACCRETION_RATE * hydro["blackHoleMass"]
             * hydro["accRateNorm"] * C_LIGHT**2)
    inner_radius = 1.1
    inner = outer * (inner_radius / hydro["rOut"])**hydro.get("s", 0.0)
    return total, outer, inner


def make_plot(config_path, spectrum_dir, output_path):
    with config_path.open("rb") as stream:
        cfg = tomllib.load(stream)
    spectrum = np.loadtxt(spectrum_dir / "lumThermal.dat")
    radial = np.loadtxt(spectrum_dir / "lumRadius.dat")
    if spectrum.ndim != 2 or spectrum.shape[1] < 9:
        raise ValueError("lumThermal.dat does not have the expected 9+ columns")
    if radial.ndim != 2 or radial.shape[1] < 7:
        raise ValueError("lumRadius.dat does not have the expected 7+ columns")

    frequency = spectrum[:, 0]
    total_lum, outer_power, inner_power = energy_budget(cfg, spectrum, spectrum_dir)
    seed_lum = integrated_luminosity(frequency, spectrum[:, 2] + spectrum[:, 3])
    compton_lum = integrated_luminosity(frequency, spectrum[:, 4])

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    ax = axes[0, 0]
    ax.loglog(frequency, positive(spectrum[:, 2]), label="Synchrotron")
    ax.loglog(frequency, positive(spectrum[:, 3]), label="Bremsstrahlung")
    ax.loglog(frequency, positive(spectrum[:, 4]), label="Comptonization")
    ax.loglog(frequency, positive(spectrum[:, 8]), color="black", linewidth=2,
              label="Total")
    ax.set(xlabel=r"Frequency $\nu$ [Hz]", ylabel=r"$\nu L_\nu$ [erg s$^{-1}$]",
           title="Thermal spectrum")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend()

    radius = radial[:, 0]  # C++ output uses Schwarzschild-radius units here.
    ax = axes[0, 1]
    ax.loglog(radius, positive(radial[:, 2]), label="Synchrotron")
    ax.loglog(radius, positive(radial[:, 3]), label="Bremsstrahlung")
    ax.loglog(radius, positive(radial[:, 4]), label="Comptonization")
    ax.loglog(radius, positive(radial[:, 6]), color="black", linewidth=2,
              label="Total")
    ax.set(xlabel=r"Radius [$R_\mathrm{S}$]", ylabel="Shell luminosity [erg s$^{-1}$]",
           title="Where the radiation is produced")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend()

    shell_total = np.maximum(radial[:, 6], 0.0)
    cumulative = np.cumsum(shell_total)
    cumulative = cumulative / cumulative[-1] if cumulative[-1] > 0.0 else cumulative
    ax = axes[1, 0]
    ax.semilogx(radius, cumulative, color="black", linewidth=2)
    ax.axhline(0.5, color="tab:gray", linestyle="--", linewidth=1)
    ax.axhline(0.9, color="tab:gray", linestyle=":", linewidth=1)
    ax.set(xlabel=r"Radius [$R_\mathrm{S}$]", ylabel="Cumulative luminosity fraction",
           title="Radial concentration", ylim=(0.0, 1.02))
    ax.grid(True, which="both", alpha=0.2)

    ax = axes[1, 1]
    labels = ["Seed", "Compton", "Total"]
    values = [seed_lum, compton_lum, total_lum]
    colors = ["tab:blue", "tab:orange", "black"]
    if outer_power is not None:
        labels.extend([r"$\dot M_\mathrm{in}c^2$", r"$\dot M_\mathrm{out}c^2$"])
        values.extend([inner_power, outer_power])
        colors.extend(["tab:green", "tab:red"])
    ax.bar(labels, values, color=colors)
    ax.set_yscale("log")
    ax.set_ylabel("Power [erg s$^{-1}$]")
    ax.set_title("Energy-budget sanity check")
    ax.tick_params(axis="x", rotation=18)
    ax.grid(True, axis="y", which="both", alpha=0.2)
    if outer_power is not None:
        efficiency = total_lum / outer_power
        ax.text(0.03, 0.97,
                (f"L / (Mdot_out c^2) = {efficiency:.3e}\n"
                 f"L / (Mdot_in c^2) = {total_lum/inner_power:.3e}"),
                transform=ax.transAxes, va="top",
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"})

    fig.suptitle("Thermal RIAF diagnostics", fontsize=15)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    summary = {
        "total_luminosity_erg_s": total_lum,
        "seed_luminosity_erg_s": seed_lum,
        "compton_luminosity_erg_s": compton_lum,
        "outer_accretion_power_erg_s": outer_power,
        "inner_accretion_power_erg_s": inner_power,
        "radiative_efficiency_outer": total_lum / outer_power if outer_power else None,
    }
    (spectrum_dir / "diagnostics.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"diagnostic plot: {output_path}")
    print(f"total thermal luminosity: {total_lum:.6e} erg s^-1")
    if outer_power:
        print(f"L_thermal / (Mdot_out c^2): {total_lum/outer_power:.6e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--spectrum-dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with args.config.open("rb") as stream:
        cfg = tomllib.load(stream)
    spectrum_dir = args.spectrum_dir or Path(
        cfg.get("run", {}).get("output_dir", "runs/example")) / "spectrum"
    output = args.output or spectrum_dir / "thermal-diagnostics.png"
    make_plot(args.config.resolve(), spectrum_dir.resolve(), output.resolve())


if __name__ == "__main__":
    main()
