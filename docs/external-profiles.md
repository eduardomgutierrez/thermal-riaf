# External radial profiles

The thermal pipeline can use radial quantities produced by another simulation.
Set the following section in the TOML configuration:

```toml
[profile]
source = "external"
file = "/path/to/profile.dat"
black_hole_mass_msun = 10.0
spin = 0.5
r_high = 3.0
```

The bundled, ready-to-run example is `examples/external-profile.toml`, using
`examples/external-profile.dat` (copied unchanged from `sol_spin_p50.dat`). Run:

```bash
python riaf_pipeline.py examples/external-profile.toml
```

The input may contain blank lines and comment lines beginning with `#`; every
data row must contain exactly 11 whitespace-separated values. Radius must be
strictly decreasing from the outer to the inner boundary. Columns are:

| Column | Quantity | Unit |
|---:|---|---|
| 1 | radius | `r_g = GM/c^2` |
| 2 | mass density | g cm^-3 |
| 3 | radial velocity | `c` |
| 4 | radial magnetic component | G |
| 5 | azimuthal magnetic component | G |
| 6 | magnetic-field magnitude | G |
| 7 | scale-height ratio `H/r` | dimensionless |
| 8 | plasma beta | dimensionless |
| 9 | ion temperature | K |
| 10 | angular velocity | code input; currently not used |
| 11 | azimuthal velocity | `c` |

The file temperature is interpreted as the ion temperature. Electron
temperature uses the beta-dependent prescription

`R(beta) = (1 + r_high beta^2) / (1 + beta^2)` and `T_e = T_i / R(beta)`.

Thus `T_e = T_i` when `r_high = 1`; for larger values, magnetically dominated
regions (`beta << 1`) remain near one temperature while gas-dominated regions
(`beta >> 1`) approach `T_e = T_i/r_high`. The older key
`ion_electron_temperature_ratio` remains accepted for compatibility, but
`r_high` is the clearer name.

`black_hole_mass_msun` must match the mass used to generate a dimensional
density profile: physical radii and emitting volumes scale with it. The bundled
`sol_spin_p50.dat` example was generated for `4.2e6` solar masses.

Before launching C++, the pipeline checks the column count, finite values,
positive density/temperature/scale height/beta, subluminal velocities,
strict radial order, and that the inner point lies outside the Kerr horizon.
It prints the accepted radial range and writes `profile-summary.json`, including
boundary accretion rates and powers. The profile is then copied into the run
directory and the same synchrotron, bremsstrahlung, and Compton calculation is
used as in Python-hydrodynamics mode.

The scattering Monte Carlo is safe to run with multiple OpenMP threads. Set
`run.omp_threads` for the desired parallelism and keep
`radiation.scattering_random_seed` fixed when comparing models. Each radial
cell receives its own deterministic random stream, so changing the thread count
does not change the resulting scattering matrix.
