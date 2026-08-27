#!/usr/bin/env python3
"""Non-interactive driver for the global transonic RIAF solution."""

import argparse
import contextlib
import io
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar


HYDRO_PLOT_FILENAMES = (
    "Temperatures.pdf",
    "MachNumber.pdf",
    "SurfaceDens.pdf",
    "accRate.pdf",
    "HR.pdf",
    "angularMom.pdf",
    "eDens.pdf",
    "magf.pdf",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def plot_hydro_diagnostics(eq, j, logr, fields, output_dir):
    """Reproduce the diagnostic PDFs from the historical hydro workflow."""
    import matplotlib.pyplot as plt

    log_ti, log_te, log_v = fields
    radius_norm = np.exp(logr)
    radius = radius_norm * eq.schwRadius
    log_radius = logr / np.log(10.0)
    ion_temperature = np.exp(log_ti) * eq.iMMW
    electron_temperature = np.exp(log_te) * eq.eMMW
    radial_velocity = -np.exp(log_v)
    sound_speed = np.sqrt(eq.sqrdSoundVel(ion_temperature, electron_temperature))
    density = eq.massDensity(radius, ion_temperature, electron_temperature,
                             radial_velocity)
    scale_height = eq.height(radius, ion_temperature, electron_temperature)
    angular_momentum = (j - eq.alpha * radius_norm * sound_speed**2
                        / eq.cLight / radial_velocity)
    gamma_i = 5.0 / 3.0
    mach_factor = np.sqrt(
        ((3.0 * gamma_i - 1.0) + 2.0 * (gamma_i - 1.0) * eq.alpha**2)
        / (gamma_i + 1.0)
    )

    def save(name, ylabel, *series, ylim=None):
        fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        for values, label, color in series:
            ax.plot(log_radius, values, label=label, color=color)
        ax.set(xlabel=r"$\log_{10}(r/R_S)$", ylabel=ylabel,
               xlim=(0.0, log_radius[-1]))
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.grid(True, alpha=0.25)
        if any(label for _, label, _ in series):
            ax.legend()
        fig.savefig(output_dir / name)
        plt.close(fig)

    save("Temperatures.pdf", r"$\log_{10}(T/\mathrm{K})$",
         (np.log10(ion_temperature), "Ion", "tab:blue"),
         (np.log10(electron_temperature), "Electron", "tab:red"))
    save("MachNumber.pdf", "Mach number",
         (-radial_velocity / sound_speed * mach_factor, None, "black"),
         ylim=(0.0, 2.2))
    save("SurfaceDens.pdf", r"$\log_{10}(\Sigma/[\mathrm{g\,cm^{-2}}])$",
         (np.log10(density * 2.0 * scale_height), None, "black"))

    accretion_radius = np.logspace(0.01, np.log10(eq.rOut), 1000) * eq.schwRadius
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.plot(np.log10(accretion_radius / eq.schwRadius),
            eq.accRateRIAF(accretion_radius) / eq.accRateOut,
            color="black", label=r"$\dot M/\dot M_\mathrm{out}$")
    ax.plot(np.log10(accretion_radius / eq.schwRadius), eq.gAux(accretion_radius),
            color="tab:red", label="g")
    ax.plot(np.log10(accretion_radius / eq.schwRadius), eq.fAux(accretion_radius),
            color="tab:blue", label="f")
    ax.set(xlabel=r"$\log_{10}(r/R_S)$", ylabel="Normalized accretion rate",
           xlim=(0.0, log_radius[-1]), ylim=(0.0, 1.0))
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(output_dir / "accRate.pdf")
    plt.close(fig)

    save("HR.pdf", r"$H/r$",
         (scale_height / radius, None, "black"), ylim=(0.0, 1.0))
    save("angularMom.pdf", r"$\log_{10}(l)$",
         (np.log10(angular_momentum), None, "black"))
    save("eDens.pdf", r"$\log_{10}(n_e/[\mathrm{cm^{-3}}])$",
         (np.log10(eq.eDens(radius, ion_temperature, electron_temperature,
                            radial_velocity)), None, "black"))
    magnetic_field = np.sqrt(
        8.0 * np.pi * (1.0 - eq.beta) * density * sound_speed**2)
    save("magf.pdf", r"$\log_{10}(B/\mathrm{G})$",
         (np.log10(magnetic_field), None, "black"))


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ["RIAF_CONFIG_JSON"] = str(args.config.resolve())

    import global_eqs as eq

    v_out = -eq.lamda*np.sqrt(eq.sqrdSoundVel(eq.temp_i_Out, eq.temp_e_Out))
    initial = (
        np.log(eq.temp_i_Out/eq.iMMW),
        np.log(eq.temp_e_Out/eq.eMMW),
        np.log(-v_out),
    )
    chatter = contextlib.nullcontext() if args.verbose else contextlib.redirect_stdout(io.StringIO())
    with chatter:
        eigen = root_scalar(
            eq.bounds_beta,
            bracket=[eq.log10j0, eq.log10j1],
            args=initial,
            method="toms748",
            maxiter=30,
        )
    if not eigen.converged:
        raise RuntimeError("angular-momentum eigenvalue search did not converge")
    j = 10.0**eigen.root
    y0 = np.array([*initial, j])
    with chatter:
        outer = solve_ivp(
            eq.rhs_beta, (np.log(eq.rOut), 0.23), y0,
            method="LSODA", events=eq.event_beta, dense_output=True,
        )
    if not outer.success or not outer.t_events[0].size:
        raise RuntimeError("hydrodynamic integration did not reach a sonic point")
    # Continue through the critical point by linearly extrapolating the last
    # three accepted steps, matching the established interactive workflow.
    back = min(3, len(outer.t) - 1)
    if back < 1:
        raise RuntimeError("too few integration points near the sonic radius")
    dt = outer.t[-1] - outer.t[-1-back]
    slope = (outer.y[:3, -1] - outer.y[:3, -1-back]) / dt
    critical_t = outer.t[-1] + dt
    critical_y = np.array([*(outer.y[:3, -1] + slope*dt), j])
    with chatter:
        inner = solve_ivp(
            eq.rhs_beta, (critical_t, np.log(1.1)), critical_y,
            method="LSODA", dense_output=True,
        )
    if not inner.success:
        raise RuntimeError(f"inner hydrodynamic integration failed: {inner.message}")

    logr = np.flip(np.concatenate((outer.t, inner.t)))
    fields = [np.flip(np.concatenate((outer.y[i], inner.y[i]))) for i in range(3)]
    np.savetxt(
        args.output_dir / "adafFile.txt",
        np.column_stack([logr, *fields]), fmt="%7.5f",
        header=str(logr.size), comments="",
    )
    np.savetxt(
        args.output_dir / "adafParameters.txt",
        (eq.blackHoleMass, eq.accRateNorm, eq.s, eq.beta, eq.alpha, j,
         eq.delta, eq.innerRadiusSSD, eq.pIndex),
    )
    plot_hydro_diagnostics(eq, j, logr, fields, args.output_dir)
    print(f"angular momentum l_in = {j:.10g}")
    print(f"sonic radius = {np.exp(outer.t[-1]):.8g} R_S")
    print(f"wrote {logr.size} radial samples to {args.output_dir}")
    print(f"wrote {len(HYDRO_PLOT_FILENAMES)} hydro diagnostic PDFs to {args.output_dir}")


if __name__ == "__main__":
    main()
