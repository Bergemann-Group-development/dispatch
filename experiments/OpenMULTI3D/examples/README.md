# Examples

Run each namelist from the experiment directory, where its relative input paths resolve:

```
$ OMP_STACKSIZE=256M ./dispatch.x examples/sun_ba_1d.nml
```

The run is written to `data/examples/<name>/` and read with `tumult.read`.

| Namelist | Content | Input |
|---|---|---|
| `sun_ba_1d.nml` | Ba II lines in LTE and NLTE, six-level Ba atom, 1D MARCS Sun | in the repository; under a second |
| `sun_ba_3d.nml` | the same atom on a 3D solar snapshot sampled to 30 × 30 × 96 cells | the M3DIS snapshot ("Optional: a 3D snapshot" in QUICKSTART.md); minutes on several cores |

[QUICKSTART.md](../QUICKSTART.md) walks through both.
