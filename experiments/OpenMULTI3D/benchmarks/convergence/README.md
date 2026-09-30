# Convergence benchmarks

These namelists test how the NLTE iteration converges for `atom.mg86g` and `atom.h20` in
the five MARCS reference models of the API benchmarks. Each run iterates to
`convlim=1e-5` or 400 iterations, whichever comes first. Ng acceleration and the output
pass are switched off, so the iteration table shows the plain preconditioned iteration.

| Star | MARCS model |
|---|---|
| `mdw` | p4000 g4.5 [Fe/H]=0 |
| `sun` | p5777 g4.4 [Fe/H]=0 |
| `tos` | p6250 g3.5 [Fe/H]=-5 |
| `rsg` | s4000 g1.0 [Fe/H]=0, 15 solar masses |
| `rgb` | s5000 g2.0 [Fe/H]=-3 |

Other settings follow `departure_solve` in `benchmarks/api`:

- abundances from `abund_magg22`, scaled by the model's [Fe/H] (`relative_abundance=T`);
- the `absmet` background opacity;
- [Fe/H] and microturbulence taken from the model files.

Unlike `departure_solve`, the runs use the atoms' own frequency quadrature. They therefore
need no generated line grid.

Run them from the OpenMULTI3D directory. dispatch names the run folder after the namelist
argument, so the output goes to `data/benchmarks/convergence/<name>`:

```sh
OMP_NUM_THREADS=4 ./dispatch.x benchmarks/convergence/mg86g_rsg.nml
for f in benchmarks/convergence/*.nml; do OMP_NUM_THREADS=4 ./dispatch.x $f > ${f%.nml}.log; done
```

The iteration table in the standard output gives the largest relative population change
per iteration. It also gives the level and depth index of that change, with index 1 at the
bottom of the model. A healthy run falls steadily to `convlim`. A run that levels off and
wanders at a fixed change has a precision floor rather than slow convergence.

## Reference behaviour

Branch `rt-column-tasks`, 30 September 2026, 4 threads, Apple M1 on AC power. The changes
of the last iterations are all below 0.005%.

| Namelist | Iterations to 1e-5 | Wall time |
|---|---|---|
| `mg86g_mdw` | 37 | 2 s |
| `mg86g_sun` | 52 | 3 s |
| `mg86g_tos` | 69 | 4 s |
| `mg86g_rsg` | 79 | 5 s |
| `mg86g_rgb` | 45 | 3 s |
| `h20_mdw` | 9 | 6 s |
| `h20_sun` | 36 | 22 s |
| `h20_tos` | 44 | 27 s |
| `h20_rsg` | 16 | 10 s |
| `h20_rgb` | 57 | 28 s |

`mg86g_rsg` is the sensitive case. The core of the Mg I 2852 A resonance line in its
outer layers is so thick that its escape probability is close to the fp32 rounding of
one. Before the spontaneous factor of the `'largest'` rate formulation was formed without
cancellation, this run wandered at 0.2 to 0.5% per iteration for all 400 iterations (see
`docs/FAQ.md`).
