"""Convert isotope/species labels between Turbospectrum and molecule_mod format.

Two conventions for tagging (isotopically resolved) molecular species meet in
this code base:

* **Turbospectrum** encodes the whole species in one token ``'CCDD.aaabbb'``.
  The part before the dot is the concatenation of the constituent atomic
  numbers, two digits each, in ascending-Z order; the part after the dot is the
  concatenation of the constituent mass numbers, three zero-padded digits each,
  in that same order.  CO(12,16) -> ``'0608.012016'``, H2O -> ``'010108.001001016'``,
  atomic Fe-56 -> ``'26.056'``.

* **molecule_mod.f90** (see ``generate_isotopic_combinations`` and the ``iso``
  index used throughout ``linelist.py``) tags a molecule only by its mass
  numbers: the constituent nuclei's integer masses, **sorted ascending**,
  written back to back with no separator and no leading zeros (Fortran ``I0``).
  CO(12,16) -> ``'1216'``, H2O -> ``'1116'``, C-13/O-16 CO -> ``'1316'``.

The molecule_mod label carries no element identity and no separators, so the
molecule -> Turbospectrum direction needs to be told which molecule it is (its
constituent elements).  ``species_elements`` derives those by parsing the species
name as a chemical formula against the element table (``_parse_formula``), so no
hardcoded molecule list is required; ``molecule_dict`` is now only a convenience
for recovering a canonical *name* on the reverse direction, which the Turbospectrum
token alone otherwise determines.

The canonical intermediate here is a list of ``(Z, A)`` pairs, one per
constituent nucleus (atomic number, mass number).  Everything else is a
parse/format around that.
"""

from itertools import permutations

from .auxfuncs import atomic_weights, molecule_dict

# atomic number <-> symbol, from the ordering of atomic_weights (H is 1)
_ELEMENTS = list(atomic_weights.keys())
_Z = {sym: i + 1 for i, sym in enumerate(_ELEMENTS)}
# main (most abundant) isotope mass number, good to the nearest integer weight
_MAIN_MASS = {sym: round(w) for sym, w in atomic_weights.items()}
# ascending-Z tuple of constituents -> molecule name, for TS -> name recovery
_MOL_BY_ZS = {
    tuple(sorted(_Z[e] for e in elems)): name
    for name, elems in molecule_dict.items()
}


def _symbol(z):
  return _ELEMENTS[z - 1]


def element_z(symbol):
  """Atomic number for an element symbol: 'Fe' -> 26."""
  return _Z[symbol]


def _parse_formula(formula):
  """Split a chemical formula into element symbols: 'H2O' -> ['H','H','O'].

  Case-sensitive, greedy longest-match against the element table, so 'CO'
  parses as C+O (the symbol 'Co' needs a lower-case 'o') and 'ZrO' as Zr+O.
  Digit multipliers are expanded ('C2' -> ['C','C'])."""
  elems = []
  i = 0
  while i < len(formula):
    if formula[i:i + 2] in _Z:      # two-letter symbol (e.g. 'Ti', 'Zr')
      sym, i = formula[i:i + 2], i + 2
    elif formula[i:i + 1] in _Z:    # one-letter symbol (e.g. 'C', 'O')
      sym, i = formula[i:i + 1], i + 1
    else:
      raise ValueError(f"cannot parse formula {formula!r} at {formula[i:]!r}")
    j = i
    while j < len(formula) and formula[j].isdigit():
      j += 1
    mult = int(formula[i:j]) if j > i else 1
    i = j
    elems.extend([sym] * mult)
  return elems


def species_elements(species):
  """Constituent element symbols of a species: 'CO' -> ['C', 'O']."""
  species = species.split()[0]  # drop any ionisation stage, e.g. 'Fe I'
  if species in atomic_weights:
    return [species]
  if species in molecule_dict:    # canonical listing for the common molecules
    return list(molecule_dict[species])
  return _parse_formula(species)  # general fallback (ZrO, VO, SiS, CrH, ...)


def species_name(constituents):
  """Molecule/atom name for a list of (Z, A) constituents, or None if unknown."""
  zs = tuple(sorted(z for z, _ in constituents))
  if len(zs) == 1:
    return _symbol(zs[0])
  return _MOL_BY_ZS.get(zs)


#=============================================================================
# Turbospectrum token <-> (Z, A) constituents
#=============================================================================
def parse_ts(code):
  """Parse a Turbospectrum species token into a list of (Z, A) constituents.

  ``'0608.012016'`` -> ``[(6, 12), (8, 16)]``; ``'26.056'`` -> ``[(26, 56)]``.
  A missing or all-zero mass part yields A = 0 (natural composition).
  """
  code = str(code).strip().strip("'").strip()
  head, _, tail = code.partition('.')

  head = head.zfill(len(head) + (len(head) % 2))  # pad to whole 2-digit atoms
  zs = [int(head[i:i + 2]) for i in range(0, len(head), 2)]
  n = len(zs)

  if tail == '':
    masses = [0] * n
  elif n == 1:
    masses = [int(tail)]
  else:
    tail = tail.zfill(3 * n)  # pad to whole 3-digit mass numbers
    masses = [int(tail[i:i + 3]) for i in range(0, len(tail), 3)]

  return list(zip(zs, masses))


def format_ts(constituents):
  """Format (Z, A) constituents as a Turbospectrum token (ascending Z)."""
  cons = sorted(constituents, key=lambda za: za[0])
  n = len(cons)
  if n == 1:
    z, a = cons[0]
    return f"{z}.{a:03d}" if a else str(z)
  head = ''.join(f"{z:02d}" for z, _ in cons)
  tail = ''.join(f"{a:03d}" for _, a in cons)
  return f"{head}.{tail}"


#=============================================================================
# molecule_mod label <-> (Z, A) constituents
#=============================================================================
def format_molec(constituents):
  """Format (Z, A) constituents as a molecule_mod label: masses ascending, I0."""
  return ''.join(str(a) for a in sorted(a for _, a in constituents))


def _partitions(digits, k):
  """Yield every split of ``digits`` into exactly ``k`` contiguous pieces, each
  1..3 digits with value >= 1 (i.e. each a plausible mass number)."""
  n = len(digits)

  def rec(pos, parts):
    if len(parts) == k:
      if pos == n:
        yield list(parts)
      return
    left = k - len(parts) - 1
    for L in (1, 2, 3):
      if pos + L > n:
        break
      rem = n - (pos + L)
      if rem < left or rem > 3 * left:  # remaining pieces must still fit
        continue
      val = int(digits[pos:pos + L])
      if val >= 1:
        parts.append(val)
        yield from rec(pos + L, parts)
        parts.pop()

  yield from rec(0, [])


def _assign_masses(digits, main):
  """Resolve an unpunctuated mass-number label against constituents ``main``.

  The molecule_mod label glues mass numbers together with no separators, and --
  as the master line lists show -- may list *fewer* masses than there are
  constituents: '50TiO' tags the isotopologue by the metal alone ('50'), the
  oxygen left unspecified.  So we consider every split into 1..n pieces, assign
  each piece to a distinct constituent, and keep the combination closest (in
  summed |mass - main|) to the main isotopes.  Constituents no piece is assigned
  to come back as 0 -- the Turbospectrum "natural / unspecified" mass, which it
  fills from the main isotope itself ('50TiO' -> O = 000, Ti = 050).  A full
  label ('1216') recovers the exact split at zero cost.
  """
  n = len(main)
  best = None  # (cost, -k, masses)

  for k in range(1, n + 1):
    for pieces in _partitions(digits, k):
      for combo in permutations(range(n), k):
        cost = sum(abs(pieces[t] - main[combo[t]]) for t in range(k))
        # tie-break toward the more explicit reading (more masses listed)
        key = (cost, -k)
        if best is None or key < best[:2]:
          masses = [0] * n  # unlisted constituents -> natural (0)
          for t in range(k):
            masses[combo[t]] = pieces[t]
          best = (cost, -k, masses)

  if best is None:
    raise ValueError(f"cannot parse isotope label {digits!r} for {n} constituents")
  return best[2]


def parse_molec(label, elements):
  """Turn a molecule_mod label + its elements into (Z, A) constituents.

  ``elements`` is a list of element symbols (or atomic numbers), e.g. from
  ``species_elements``.  A label of ``'0'`` (or empty) means natural
  composition: every constituent left unspecified (mass 0), which Turbospectrum
  fills from the main isotope.
  """
  syms = [_symbol(e) if isinstance(e, int) else e for e in elements]
  # the label lists masses ascending, so pair them with the atoms ordered the
  # same way (ascending weight) before splitting
  atoms = sorted(syms, key=lambda s: atomic_weights[s])
  main = [_MAIN_MASS[s] for s in atoms]

  label = str(label).strip().strip("'")
  if label.strip('0') == '':  # '', '0', '000', ... -> all natural (mass 0)
    masses = [0] * len(atoms)
  else:
    masses = _assign_masses(label, main)

  return [(_Z[s], m) for s, m in zip(atoms, masses)]


#=============================================================================
# High-level, name-aware conversions
#=============================================================================
def ts_to_molec(code):
  """Turbospectrum token -> (species_name_or_None, molecule_mod_label).

  >>> ts_to_molec('0608.012016')
  ('CO', '1216')
  >>> ts_to_molec('010108.001001016')
  ('H2O', '1116')
  """
  cons = parse_ts(code)
  return species_name(cons), format_molec(cons)


def molec_to_ts(species, label):
  """Species name + molecule_mod label -> Turbospectrum token.

  >>> molec_to_ts('CO', '1216')
  '0608.012016'
  >>> molec_to_ts('H2O', '1116')
  '010108.001001016'
  """
  cons = parse_molec(label, species_elements(species))
  return format_ts(cons)


if __name__ == '__main__':
  # round-trip check over the molecules we know, main + a swapped-isotope case
  cases = [
      ('CO', '1216'), ('CO', '1316'), ('CO', '1218'),
      ('H2O', '1116'), ('OH', '116'), ('CH', '112'),
      ('C2', '1212'), ('TiO', '1648'), ('SiO', '1628'),
      ('FeH', '156'), ('MgH', '124'), ('CN', '1214'),
  ]
  for name, label in cases:
    ts = molec_to_ts(name, label)
    back_name, back_label = ts_to_molec(ts)
    ok = (back_label == label) and (back_name == name)
    print(f"{name:4s} {label:6s} <-> {ts:20s} -> ({back_name}, {back_label})"
          f"  {'ok' if ok else 'MISMATCH'}")
