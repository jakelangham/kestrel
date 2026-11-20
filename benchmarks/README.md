# Benchmarks

These measure Kestrel's accuracy against exact solutions, and how its errors
converge as the grid is refined. Unlike the tests in `../tests`, they have no
pass/fail tolerances, and they are not run by `make check` or CI.

## Running

Run them from this directory, which needs Julia and a built `../src/kestrel`:

```sh
julia runbench.jl                           # all cases, each at its default levels
julia runbench.jl 1d                        # the 1D cases (or 2d)
julia runbench.jl thacker_radial_2d         # one case
julia runbench.jl help                      # the options, and a summary of each case
julia runbench.jl --levels 2                # only the two coarsest grids
julia runbench.jl --markdown results.md     # also write the table as Markdown
julia runbench.jl --reference fb4bf4b      # compare with another commit
julia runbench.jl --save errors.txt         # also write the errors to a file
```

For each case, level 0 is the grid in `Input_<case>.txt`, and each further
level halves the grid spacing. Each run writes its output to `<case>_n<cells
per tile>/`, which is not tracked by git. Set `OMP_NUM_THREADS` for OpenMP
builds; the results do not depend on it. All levels of all cases take about 6
minutes on 8 threads, most of it the finest 2D grids.

The runner reports each run as it finishes, in the style of
`../tests/runall.jl`: the depth and discharge errors (defined below), their
observed orders of convergence from the previous grid, `log2(e_coarse /
e_fine)`, the reference's errors on the same grid and the ratios to them, the
simulation time,
and a status:

- `OK`: both ratios are within 10% of 1, or below;
- `WORSE`: either ratio exceeds 1.1;
- `BETTER`: both ratios are below 1/1.1;
- `NEW`: there is no reference for this case and grid;
- `FAIL`: the run failed or produced a negative depth.

It ends with a summary table of all runs and the number that are neither
`WORSE` nor `FAIL`.

The exit status is 1 if any run is `WORSE` or `FAIL`.

## The reference

By default the errors are compared with those of the head of `main`, so the
ratios show what a branch changes. `--reference` names another git ref, a
file written by `--save` (for instance from uncommitted code), or `none`.

The errors of a commit are kept in `references/<short hash>.txt`, which is not
tracked by git. If that file is missing, `runbench.jl` calls
`record_reference.sh`, which builds the commit, unchanged, in a temporary git
worktree, runs every case at its default levels there with this checkout's
harness and inputs, and writes the file. This takes several minutes, once per
commit; when `main` moves on, its new head is recorded on the next run.

The inputs are frictionless (`Drag = None`). A commit without `Drag = None`,
such as `main` before the `fix-timestep-rejection` branch is merged, runs them
with `Drag = Chezy` and `Chezy Co = 1e-8`, the smallest drag it accepts. That
changes the errors by up to about 7% from frictionless flow; the reference
file and the runner's output say so.

## Cases

### Thacker's oscillations

Three frictionless oscillations in a paraboloid bowl
`b = h0 (x^2 + y^2) / a^2` with `h0 = 10 m` and `a = 3000 m`, from the exact
solutions of Thacker (1981, *J. Fluid Mech.* 107, 499–508). Each run lasts two
periods, with wet/dry fronts moving throughout. They test smooth flow with
moving shorelines over curved topography.

| Case | Dimensions | Solution | Grids |
|---|---|---|---|
| `thacker_planar` | 1D | Planar free surface; the lake's centre oscillates by ±500 m (period 1346 s) | 50, 25, 12.5, 6.25 m |
| `thacker_planar_2d` | 2D | Planar free surface; the lake's centre goes round a circle of radius 500 m (period 1346 s) | 200, 100, 50 m |
| `thacker_radial_2d` | 2D | Radially symmetric, curved free surface; the shoreline's radius oscillates between 2280 and 3948 m (period 673 s) | 200, 100, 50 m |

### Dam breaks

A reservoir 10 m deep on `-1000 m < x < 0`, released at t = 0. They test a
sharp wet/dry front, a shock, and Coulomb friction and the geometric factors
on a slope. All are 1D, on grids of 10, 5, 2.5 and 1.25 m.

| Case | Solution |
|---|---|
| `dambreak_dry` | Onto a dry flat bed, without friction (Ritter 1892): a rarefaction, and a front at speed `2 sqrt(g h)`. Run to 40 s. |
| `dambreak_wet` | Onto a layer 2 m deep on `0 < x < 800 m`, flat and without friction (Stoker 1957): a rarefaction into the reservoir and a shock into the layer. Run to 40 s. |
| `dambreak_coulomb_slope` | Onto the dry slope `b = -0.5 x`, with Coulomb friction coefficient 0.3, geometric factors on (Mangeney, Heinrich & Roche 2000, *Pure Appl. Geophys.* 157, 1081–1096), released from rest. Run to 30 s, by when the flow has moved 700 m downhill. |
| `dambreak_coulomb_slope_moving` | The same, but with the reservoir already moving downhill at 10 m/s at t = 0. |

On a uniform slope, Kestrel's 1D equations are the shallow water equations in
the horizontal coordinate, with gravity `g cos^3 θ` and a uniform acceleration
`g (sin θ - μ cos θ) cos θ` where the flow moves downhill (the equations of
Savage and Hutter, written in the horizontal coordinate). A frame moving with
that acceleration removes it, leaving Ritter's solution (the derivation is in
`benchlib.jl`, above `dambreak_constants`). This needs the slope to be steeper
than the friction angle, so that the flow never moves uphill.

The two Coulomb cases differ in how the flow starts. Released from rest, the
direction of the friction is undefined at first, and Kestrel's implicit drag,
whose coefficient is proportional to `1 / |u|`, gets the first steps' friction
wrong while `u` changes by much of itself within a step. Since Coulomb
friction doesn't depend on speed, the resulting speed error persists: at
10 m, the whole reservoir moves 0.17 m/s too fast (1% at 30 s), which
dominates the error and makes it depend on the first few step sizes (a
rejected step changes it by tens of percent). It converges at about first
order (0.033 m/s at 1.25 m). Started at 10 m/s, the speed error is 0.026 m/s
at 10 m and falls at about third order, so `dambreak_coulomb_slope_moving`
measures the front and the rarefaction, and the difference between the two
cases the start from rest.

The reservoir and the layer are finite, so disturbances come in from their far
ends. The errors are measured only in the region they have not yet reached,
at most the speed `|v| + sqrt(g h)` in the moving frame, less a further
100 m. In the dry cases this keeps the whole rarefaction and the front; in the
wet case, also the shock.

Dam breaks converge at about first order, because the exact solutions are not
smooth: Ritter's has kinks at the front and at the head of the rarefaction,
and Stoker's a shock.

### Common to all

The exact solutions, initial conditions and errors are in `benchlib.jl`. It
includes `../tests/testlib.jl` for running Kestrel and reading its results,
and is the place for other non-critical checks: ones that measure accuracy or
behaviour without a firm pass/fail criterion, and so are kept out of the test
suite, which `make check` and CI run.

- **Initial condition.** Kestrel can't set up these states from an input
  block, so `make_exact_ic` writes one, and the input files load it with
  `Initial condition`. To get the tile layout, it first runs the input with a
  Cap or Cube block in place of the initial-condition file. It then replaces
  the depth and discharge in the t = 0 output with the exact solution,
  averaged over each cell on Kestrel's bed. That bed is the exact one at cell
  vertices, interpolated linearly in 1D and bilinearly in 2D.
- **Errors.** At each output time, the depth error is
  `sum |h - h_exact| / sum h_exact` over cells (over the region above, for dam
  breaks). The discharge error is `sum |q - q_exact| / (U sum h_exact)`, where
  `q` is the discharge vector and `U` a speed scale: the largest exact speed at
  the shoreline for Thacker's solutions, and the front's speed at the end of
  the run (behind the shock, for Stoker's) for dam breaks. The reported error
  is the largest over the output times. Depths are the bed-normal depth `Hn`,
  and velocities horizontal.

## Results at `fb4bf4b`

| Case | dx (m) | Depth error | Order | Discharge error | Order |
|---|---|---|---|---|---|
| `thacker_planar` | 50 | 2.79e-4 | | 1.00e-3 | |
| | 25 | 1.11e-4 | 1.33 | 3.38e-4 | 1.57 |
| | 12.5 | 3.34e-5 | 1.73 | 1.26e-4 | 1.42 |
| | 6.25 | 1.08e-5 | 1.63 | 3.69e-5 | 1.77 |
| `thacker_planar_2d` | 200 | 3.88e-3 | | 2.14e-2 | |
| | 100 | 1.22e-3 | 1.67 | 8.05e-3 | 1.41 |
| | 50 | 3.99e-4 | 1.62 | 2.80e-3 | 1.52 |
| `thacker_radial_2d` | 200 | 3.26e-2 | | 2.75e-2 | |
| | 100 | 9.57e-3 | 1.77 | 1.01e-2 | 1.45 |
| | 50 | 3.71e-3 | 1.37 | 3.44e-3 | 1.55 |

For comparison, the wet/dry correction of Chertock et al. (2015), which
preceded the Kurganov–Petrova correction now used, gives errors 1.5 to 3.5
times larger on every grid, with similar orders.
