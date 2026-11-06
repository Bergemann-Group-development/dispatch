---
name: python-style
description: Naming and style rules for Python in OpenMULTI3D (tumult/, TumultStudio/, scripts, and code shown in docs/manual). Use whenever writing or editing Python in this experiment, including example code in the documentation.
---

# Python style

The policy is `docs/PythonNamingConventions.md`; read its "Variables" section before
writing code. `docs/GLOSSARY.md` says what names mean, and `CLAUDE.md` has the package
preferences (pathlib, dimensionless internals, astropy units at the boundaries).

## Checklist for each edit

1. **Locals use the abbreviation table** in `PythonNamingConventions.md` (`win`, `lvl`,
   `abnd`, `obs`, `syn`, `cnt`, `dep`, `err`, `cfg`, `tmp`, `tran`, `bb`/`bf`, `llist`),
   with `s` appended for a collection (`wins`, `abnds`, `trans`).
2. **`atom` and `atmos` are spelled out.** `star` is a set of stellar parameters.
3. **Signatures stay spelled out.** Parameters, keyword arguments, attributes and other
   public names keep full words (`abundance=`, `window=`); abbreviate only local and
   temporary variables.
4. **Comprehensions use the item's initial**: `[r.crop() for r in runs]`,
   `[a + 0.1 for a in abnds]`. Arrays outside comprehensions keep multi-character names
   (`xx`, `yy`).
5. **Touching old code**: rename long locals only inside functions already being edited
   for another reason; no rename-only sweeps.
6. **Manual code** follows the same rules, and code shown in the manual matches
   the script it is taken from (e.g. `docs/manual/figs/science_case.py`).

When a name is missing from the table and recurs, ask the user before adding it, then add
it to the table in `PythonNamingConventions.md`.
