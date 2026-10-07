# API benchmarks

`intended_api.py` holds the benchmark calculations of `docs/API_benchmarks.txt` as the calls
a user writes: `single_lte`, `departure_solve`, `single_precomputed_nlte`,
`single_ondemand_nlte`, `nlte_correction`, `abundance_fit`, `ionisation_balance` and the
composed `stellar_analysis`. Each function takes public inputs and can be imported on its
own:

```python
from benchmarks.api.intended_api import single_lte
from tumult.synthesis import AtmosphereGrid, Star, SynthConfig

config = SynthConfig(root=".", linelist="input/linelists/sun_3000-9000.txt",
                     abundances="input/abund/abund_magg22", R_sampling=100_000)
result = single_lte(config, AtmosphereGrid("input/atmos/MARCS"),
                    Star(5777, 4.44, 0.0, 1.0), (5160, 5190))
```

`run_reference_benchmarks.py` binds the reference assets: the five supplied MARCS reference
atmospheres, `atom.mg86g`, the solar line list, the Fe masks and a default 4000–9000 Å
range. It preflights what the selected workflows need, times each calculation and writes a
`benchmark.json` beside its saved result, reporting computational completion separately
from the scientific outcome (termination, convergence, exclusions by flag). The target is
the KPNO atlas in `input/spectra/` (`python -m tumult.synthesis.input_data fetch kpno`):

```sh
python -m benchmarks.api.run_reference_benchmarks \
  --target-spectrum input/spectra/KPNO_4000_9000.txt --preflight

python -m benchmarks.api.run_reference_benchmarks \
  --target-spectrum input/spectra/KPNO_4000_9000.txt --only single_lte --wavelength 5160 5190
```

Repeat `--only` to select workflows. `single_precomputed_nlte` reads the departures saved by
`departure_solve` in the same output folder. The complete run is expensive: per reference
star an LTE synthesis, a departure solve, two NLTE syntheses and a correction ladder, then an
abundance fit, an ionisation balance and the composed analysis on the KPNO atlas. Each
calculation runs in its own cache, `data/synth/api_benchmarks/<benchmark>[/<star>]` for native
runs and `synth/api_benchmarks/<benchmark>[/<star>]` for preparations, so no calculation reuses
another's results. The one exception is `single_precomputed_nlte`, which reads the departures of
`departure_solve`.

## Standard timings

Cold run of every benchmark with the defaults above, on an Apple M1 (4 threads, AC power),
commit `2589b461` (wall seconds, summed over the five stars where applicable):

| Benchmark | Wall |
|---|---:|
| `single_lte` | 28 |
| `departure_solve` | 15 |
| `single_precomputed_nlte` | 32 |
| `single_ondemand_nlte` | 48 |
| `nlte_correction` | 204 |
| `abundance_fit` | 51 |
| `ionisation_balance` | 84 |
| `stellar_analysis` | 575 |
| **Suite** | **1037** |

Timings are single samples. On battery power the M1 takes about twice as long. Compare new
timings against these only from a cold run with the same settings.
