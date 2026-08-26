# Thermal RIAF pipeline

This repository joins the one-dimensional transonic RIAF solver in `RIAF/`
to the thermal radiation code in `radproc/`. A collaborator edits one TOML
file and runs one command; the pipeline finds the angular-momentum eigenvalue,
integrates through the sonic point, transfers the radial solution, builds the
C++ program, and calculates synchrotron, bremsstrahlung, and Compton emission.

## Requirements

- Python 3.11 or newer
- NumPy, SciPy, Astropy, and Matplotlib
- A C++17 compiler, Meson, and Ninja
- Boost, GSL, HDF5 (including its C++ library), and OpenMP development files

On Ubuntu/Debian, install the native dependencies with:

```bash
sudo apt install build-essential meson ninja-build libboost-all-dev libgsl-dev libhdf5-dev
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

If Meson is not installed by the operating system, it may instead be installed
inside the active environment with `python -m pip install meson`.

## Quick start

```bash
python riaf_pipeline.py examples/thermal-riaf.toml
```

The example writes hydrodynamic files under `runs/thermal-riaf/hydro/` and the
radiative results under `runs/thermal-riaf/spectrum/`. The final spectrum is
`lumThermal.dat`. It also creates `thermal-diagnostics.png`, containing the
separate emission components, radial luminosity, cumulative emission, and an
accretion-power sanity check. Machine-readable values are saved in
`diagnostics.json`. The pipeline uses a separate `radproc/build/` directory and
reuses it on later runs. Required Compton probability tables are copied into
each run automatically.

The example defaults to one OpenMP thread. The legacy scattering-matrix routine
has shared mutable state and can fail nondeterministically with multiple threads;
do not increase `run.omp_threads` until that routine has been made thread-safe.

To test only the Python solution and parameter handoff:

```bash
python riaf_pipeline.py examples/thermal-riaf.toml --hydro-only
```

All normal model choices are documented directly in the example TOML. Copy it,
change the physical parameters, and keep `profile.source = "hydro"`. Thermal
processes have explicit Boolean names; nonthermal calculations are always
disabled by this interface.

## Radial profiles from another code

Set `profile.source = "external"` to bypass the Python hydrodynamics and feed a
radial profile directly to `radproc`. The expected columns and units are in
[docs/external-profiles.md](docs/external-profiles.md).

The repository includes a complete external-profile example:

```bash
python riaf_pipeline.py examples/external-profile.toml
```

## Existing low-level interfaces

The historical scripts and full `radproc/src/adaf/parameters.json` interface
remain available for specialized projects. New thermal RIAF runs should use the
top-level pipeline because it validates the input, sets the thermal flags
consistently, and keeps generated files out of the source directories.

## Reproducibility and troubleshooting

- Each run retains the exact generated `hydro-input.json` and
  `parameters.json` beside its outputs.
- If the angular-momentum bracket does not straddle a smooth solution, adjust
  `hydro.log10j0` and `hydro.log10j1` in the TOML.
- Delete `radproc/build/` after changing compilers or system dependencies, then
  rerun the command.
- `--no-build` runs an already compiled executable without invoking Meson.

Run the lightweight interface tests with `python -m unittest discover -s tests`.

To recreate a plot without rerunning either physical code:

```bash
python plot_diagnostics.py examples/thermal-riaf.toml
```
