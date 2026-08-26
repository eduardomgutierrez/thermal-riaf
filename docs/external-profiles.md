# External radial profiles

The thermal pipeline can use radial quantities produced by another simulation.
Set the following section in the TOML configuration:

```toml
[profile]
source = "external"
file = "/path/to/profile.dat"
black_hole_mass_msun = 10.0
spin = 0.5
ion_electron_temperature_ratio = 3.0
```

The input must have one header line followed by whitespace-separated rows. The
radius must decrease from the outer to the inner boundary. Columns are:

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

The present reader obtains the electron temperature from
`T_e = T_i / ion_electron_temperature_ratio`. The pipeline copies the profile
into the run directory, records all choices in `parameters.json`, and otherwise
uses the same synchrotron, bremsstrahlung, and Compton calculation as the
Python-hydrodynamics mode.

