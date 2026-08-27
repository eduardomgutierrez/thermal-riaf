# Configuration and parameter files

For normal runs, edit only a top-level TOML file such as
`examples/thermal-riaf.toml` or `examples/external-profile.toml`. The pipeline
creates the lower-level files consumed by the Python and C++ programs.

## `[run]`: execution and output

| Parameter | Meaning | Default |
|---|---|---:|
| `output_dir` | Run directory, relative to the repository root unless absolute. Hydro and spectrum products are placed in its `hydro/` and `spectrum/` subdirectories. | `runs/example` |
| `omp_threads` | Number of OpenMP threads used by radproc. It changes performance, not the result. | `4` |

## `[profile]`: source of the radial solution

`source = "hydro"` runs the bundled transonic solver and requires a `[hydro]`
section. `source = "external"` bypasses that solver and requires these entries:

| Parameter | Unit/range | Meaning |
|---|---|---|
| `file` | path | Eleven-column radial profile. A relative path is resolved from the TOML location. See `external-profiles.md`. |
| `black_hole_mass_msun` | solar masses, `> 0` | Mass used to convert the dimensionless profile radius to cm. It must match the mass for which a dimensional density profile was generated. |
| `spin` | dimensionless, `-1 < a < 1` | Kerr spin parameter. It sets the horizon and relativistic redshift/velocity corrections. |
| `r_high` | dimensionless, `> 0` | Ion-to-electron temperature control in the beta-dependent electron-heating prescription. It is not a constant `Ti/Te`; see `external-profiles.md`. |

## `[hydro]`: one-dimensional RIAF solution

These names retain their historical spelling because they are also understood
by `RIAF/global_initCond.py`.

| Parameter | Unit/range | Physical meaning |
|---|---|---|
| `blackHoleMass` | solar masses, `> 0` | Black-hole mass. |
| `accRateNorm` | Eddington accretion rate, `> 0` | Accretion rate at `rOut`, using `Mdot_Edd = 1.39e18 (M/Msun) g s^-1`. |
| `rOut` | Schwarzschild radii `R_S`, `> 1` | Outer boundary of the hydrodynamic integration. |
| `innerRadiusSSD` | `R_S` | Transition/inner thin-disk radius used by the optional historical SSD coupling. With the pipeline's thermal RIAF configuration (`SSDdisk = 0`), it does not activate a thin disk. |
| `beta` | dimensionless, normally `0 < beta < 1` | Gas-to-total pressure ratio. The magnetic pressure fraction is `1-beta`; values nearer one produce weaker magnetic fields. This is not the external profile's column-8 beta convention unless the supplying code uses the same definition. |
| `alpha` | dimensionless, normally `0 < alpha < 1` | Shakura–Sunyaev viscosity parameter. |
| `delta` | dimensionless, `0 <= delta <= 1` | Fraction of dissipative heating deposited directly into electrons; `1-delta` heats ions. |
| `s` | dimensionless | Wind/outflow index in `Mdot(r) = Mdot(rOut) (r/rOut)^s`. `s=0` conserves accretion rate; positive `s` reduces it inward. |
| `pIndex` | dimensionless | Shape index for the optional SSD/RIAF transition prescription; normally inactive when `SSDdisk=0`. |
| `log10j0`, `log10j1` | `log10(l_in/(R_S c))` | Lower and upper bracket for the angular-momentum eigenvalue. The bracket must straddle a solution smooth at the sonic point. These are search bounds, not the final angular momentum. |

Advanced outer-boundary quantities can also be supplied: `temp_i_Out` and
`temp_e_Out` in kelvin, and `lamda`, the magnitude of the outer radial velocity
in units of the local sound speed. If any of these three is omitted, the bundled
radius-dependent prescription sets all three. Other historical options accepted
by the hydro code are `coronaFactor`, `SSDdisk`, and `correctorAccRate`; they are
not part of the supported thermal-only interface and should normally be left at
their defaults.

## `[radiation]`: emitted components and scattering

| Parameter | Meaning | Default |
|---|---|---:|
| `synchrotron` | Include thermal synchrotron emission and self-absorption. | `true` |
| `bremsstrahlung` | Include thermal electron–ion/electron–electron bremsstrahlung. | `true` |
| `comptonization` | Include thermal Comptonization of the seed photons. | `true` |
| `calculate_scattering_matrix` | Recompute the Monte Carlo radial scattering/escape matrices. On a rerun with matching parameters and profile data, the pipeline asks whether to reuse the existing matrices. If false, compatible matrix files must already exist in the spectrum run directory. | `true` |
| `scattering_random_seed` | Non-negative master seed. Each radial cell derives an independent stream, making results reproducible across OpenMP thread counts. | `5489` |

The top-level pipeline always disables nonthermal particles, jets, hadronic
emission, and cold-disk emission, regardless of legacy options in the radproc
template.

## `[grid]`: numerical spectrum grid

| Parameter | Unit | Meaning | Default |
|---|---|---|---:|
| `radius_samples` | count | Logarithmically spaced radial cells used by radproc. More cells improve radial resolution and increase scattering-matrix cost. | `50` |
| `energy_samples` | count | Photon-energy bins. | `100` |
| `log10_energy_ev_min` | `log10(eV)` | Lower photon-energy boundary. | `-6` |
| `log10_energy_ev_max` | `log10(eV)` | Upper photon-energy boundary. | `8` |

## Generated hydro handoff files

For a hydro run, `runs/<name>/hydro/hydro-input.json` records the `[hydro]`
section exactly as supplied.

`adafFile.txt` begins with the number of radial samples. Each following row is:

| Column | Stored quantity | Recover with | Unit |
|---:|---|---|---|
| 1 | `ln(r/R_S)` | `r/R_S = exp(value)` | dimensionless |
| 2 | `ln(T_i/mu_i)` | `T_i = mu_i exp(value)` | K |
| 3 | `ln(T_e/mu_e)` | `T_e = mu_e exp(value)` | K |
| 4 | `ln(|v_r|)` | `v_r = -exp(value)` | cm s^-1 |

Here `mu_i=1.23` and `mu_e=1.14` are the ion and electron mean molecular
weights used by both codes.

`adafParameters.txt` is a legacy positional, single-row file. Its nine values
are, in order:

1. black-hole mass (`M_sun`);
2. outer Eddington-normalized accretion rate;
3. wind index `s`;
4. gas-to-total pressure ratio `beta`;
5. viscosity parameter `alpha`;
6. solved angular-momentum eigenvalue `l_in/(R_S c)`;
7. electron-heating fraction `delta`;
8. SSD transition radius (`R_S`);
9. SSD transition shape index `pIndex`.

These two legacy files should not be edited or copied manually when using the
pipeline.

## Generated radproc `parameters.json`

The pipeline writes `runs/<name>/spectrum/parameters.json`. This is the exact
low-level C++ configuration used for the run. Relevant mappings include:

| TOML setting | Generated JSON setting |
|---|---|
| `profile.source` | `readPrecomputedADAF` |
| external `file` | `solSpinFile` |
| external mass/spin/`r_high` | `blackHoleMass`, `spin`, `R_high` |
| radiation component switches | `thermal.processes` and legacy `thermal.processNumber` |
| `calculate_scattering_matrix` | `calculateComptonScatt` |
| `scattering_random_seed` | `scatt.randomSeed` |
| grid settings | `model.particle.default.dim.radius` and `model.particle.photon.dim.energy` |

The pipeline also forces `calculateThermal=1` and `calculateNonThermal=0`. Other
entries are inherited from `radproc/src/adaf/parameters.json` for compatibility
with the legacy C++ configuration. They are implementation defaults rather than
supported collaborator-facing controls. Specialized users may run radproc
directly, but then they are responsible for maintaining mutually consistent
flags, input files, Compton tables, and grid dimensions.
