"""Species and isotope labels of line lists, in Turbospectrum's convention.

Turbospectrum encodes a species in one token ``'CCDD.aaabbb'``. The part before the
dot concatenates the constituent atomic numbers, two digits each, in ascending order;
the part after it concatenates their mass numbers, three digits each, in the same
order, with 000 for the natural mixture: CO(12,16) is ``'0608.012016'``, H2O is
``'010108.001001016'`` and atomic Fe-56 is ``'26.056'``.

The isotope label of a line-list block is that mass part, in both line-list formats.
A tumult-format header names the species and writes the label on its own (``'CN'``,
``'012014'``); HOP stores it as an integer and the number of constituents restores the
leading zeros. Hydrogen has no isotope table in HOP, so its mass reads as natural
(``001`` and ``000`` are the same label), and the masses of equal atoms sort ascending
(``012013`` for 12C13C). :func:`mass_label` gives that canonical form.

The canonical intermediate here is a list of ``(Z, A)`` pairs, one per constituent
nucleus (atomic number, mass number).
"""

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
_ROMAN = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI', 7: 'VII', 8: 'VIII'}


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


#=============================================================================
# isotope labels of line-list headers
#=============================================================================
def mass_label(constituents):
    """The isotope label of (Z, A) constituents, as an integer.

  Masses in ascending Z, ascending among equal Z, three digits each; hydrogen is
  always natural (0). ``[(6, 12), (1, 1)]`` -> 12 (written ``000012``).
  """
    cons = sorted((z, 0 if z == 1 else a) for z, a in constituents)
    return int(''.join(f"{a:03d}" for _, a in cons) or '0')


def iso_token(species, iso):
    """The isotope field of a tumult-format header: ``'000012'`` for CH, ``'56'`` for Fe."""
    n = len(species_elements(species))
    return f"{int(iso):0{3 * n}d}" if n > 1 else str(int(iso))


def ts_token(species, iso):
    """Species name and isotope label -> Turbospectrum token.

  >>> ts_token('CO', 12016)
  '0608.012016'
  >>> ts_token('Fe I', 56)
  '26.056'
  """
    zs = sorted(element_z(e) for e in species_elements(species))
    if len(zs) == 1:
        return f"{zs[0]}.{int(iso):03d}"
    return ''.join(f"{z:02d}" for z in zs) + '.' + f"{int(iso):0{3 * len(zs)}d}"


def ts_species(code, ion):
    """Turbospectrum token and ionisation stage -> (species name, isotope label).

  >>> ts_species('26.056', 1)
  ('Fe I', 56)
  >>> ts_species('0607.012014', 1)
  ('CN', 12014)
  """
    cons = parse_ts(code)
    name = species_name(cons)
    if name is None:
        raise ValueError(f"unknown molecule in Turbospectrum species {code!r}")
    if len(cons) == 1:
        return f"{name} {_ROMAN[int(ion)]}", cons[0][1]
    return name, mass_label(cons)
