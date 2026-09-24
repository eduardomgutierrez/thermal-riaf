#!/usr/bin/env python3
"""Non-interactive driver for the global transonic RIAF solution."""

import argparse
import contextlib
import io
import json
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


class _TrialStuck(Exception):
    """Inner integration is crawling at the sonic singularity."""


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
    cfg = json.loads(args.config.read_text())
    if "eigenvalue" in cfg:
        # The eigenvalue function is not monotonic, so a bracketed root search
        # can land on a spurious root; a known eigenvalue can be given instead.
        j = float(cfg["eigenvalue"])
        print(f"using fixed angular-momentum eigenvalue j = {j:.10g}")
    else:
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

    def mach(y):
        temp_i, temp_e = np.exp(y[0])*eq.iMMW, np.exp(y[1])*eq.eMMW
        return np.exp(y[2])/np.sqrt(eq.sqrdSoundVel(temp_i, temp_e))

    log_r_in = np.log(1.1)
    with chatter, np.errstate(all="ignore"):
        full = solve_ivp(eq.rhs_beta, (np.log(eq.rOut), log_r_in), y0,
                         method="LSODA", dense_output=True)
    grid = np.linspace(np.log(eq.rOut), log_r_in, 4000)
    supersonic = np.nonzero(mach(full.sol(grid)) >= 1.0)[0]
    if not supersonic.size:
        raise RuntimeError("hydrodynamic integration did not reach a sonic point")
    log_r_sonic = grid[supersonic[0]]

    # Restart inward from the sonic point, extrapolating across it by a small
    # jump in ln r ("skipping a few steps"); accept the first restart whose
    # solution stays supersonic down to the inner boundary.
    slope_step = float(cfg.get("sonic_slope_step", 0.02))
    y_sonic = full.sol(log_r_sonic)
    slope = (full.sol(log_r_sonic + 3*slope_step)[:3] - y_sonic[:3]) / (3*slope_step)
    inner = None
    best_partial = None  # (t, y) of the deepest finite supersonic segment
    for skip in cfg.get("sonic_skips", (0.02, 0.04, 0.06, 0.1, 0.15, 0.2, 0.3)):
        start_t = log_r_sonic - skip
        start_y = np.array([*(y_sonic[:3] - slope*skip), j])
        calls = [0]

        def rhs_capped(logr, y):
            calls[0] += 1
            if calls[0] > 20000:
                raise _TrialStuck
            return eq.rhs_beta(logr, y)

        try:
            with chatter, np.errstate(all="ignore"):
                trial = solve_ivp(rhs_capped, (start_t, log_r_in), start_y,
                                  method="LSODA", dense_output=True)
        except _TrialStuck:
            print(f"sonic skip {skip:g}: stuck at the singularity")
            continue
        finite = np.isfinite(trial.y).all(axis=0)
        ok = np.cumprod(finite & (np.nan_to_num(mach(np.nan_to_num(trial.y))) > 1.0)).astype(bool)
        if ok.sum() >= 5 and (best_partial is None or trial.t[ok][-1] < best_partial[0][-1]):
            best_partial = (trial.t[ok], trial.y[:, ok])
        good = (trial.success and np.isfinite(trial.y).all()
                and trial.t[-1] <= log_r_in + 1e-6 and mach(trial.y[:, -1]) > 1.0)
        print(f"sonic skip {skip:g}: {'accepted' if good else 'rejected'} "
              f"(Mach at inner boundary {mach(trial.y[:, -1]):.2f})")
        if good:
            inner = trial
            break
    if inner is None:
        # Fallback: extrapolate ln T_i, ln T_e and ln|v| linearly in ln r (power
        # laws) from the deepest valid supersonic segment down to the inner
        # boundary. The velocity slope is capped at free fall, v ∝ r^(-1/2).
        if best_partial is not None:
            seg_t, seg_y = best_partial
        else:
            seg_t = np.linspace(log_r_sonic + 5*slope_step, log_r_sonic, 6)
            seg_y = full.sol(seg_t)
        n_fit = min(8, seg_t.size)
        slopes = [np.polyfit(seg_t[-n_fit:], seg_y[k, -n_fit:], 1)[0] for k in range(3)]
        slopes[2] = max(slopes[2], -0.5)
        ext_t = np.linspace(seg_t[-1], log_r_in, 30)[1:]
        ext_y = np.array([seg_y[k, -1] + slopes[k]*(ext_t - seg_t[-1]) for k in range(3)]
                         + [np.full_like(ext_t, j)])
        inner = argparse.Namespace(t=np.concatenate((seg_t, ext_t)),
                                   y=np.concatenate((seg_y, ext_y), axis=1))
        print(f"WARNING: no restart converged; profiles extrapolated as power laws from "
              f"r = {np.exp(seg_t[-1]):.3f} R_S to the inner boundary "
              f"(d ln T_i/d ln r = {slopes[0]:.2f}, d ln T_e/d ln r = {slopes[1]:.2f}, "
              f"d ln v/d ln r = {slopes[2]:.2f})")

    keep = full.t > log_r_sonic
    outer = argparse.Namespace(t=full.t[keep], y=full.y[:, keep])

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
    print(f"sonic radius = {np.exp(log_r_sonic):.8g} R_S")
    print(f"wrote {logr.size} radial samples to {args.output_dir}")
    print(f"wrote {len(HYDRO_PLOT_FILENAMES)} hydro diagnostic PDFs to {args.output_dir}")


if __name__ == "__main__":
    main()
