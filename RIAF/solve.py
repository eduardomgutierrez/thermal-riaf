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


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


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
    print(f"angular momentum l_in = {j:.10g}")
    print(f"sonic radius = {np.exp(outer.t[-1]):.8g} R_S")
    print(f"wrote {logr.size} radial samples to {args.output_dir}")


if __name__ == "__main__":
    main()
