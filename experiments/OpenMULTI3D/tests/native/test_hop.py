"""The standalone HOP library (libhop.so) through the tumult.hop wrapper.

Needs ``HOP/standalone/build/libhop.so`` (``make hop_standalone``). Chunking is only
neutral within one process: a spawned worker reads only the lines within dlam of its own
chunk, so a far Lorentz wing reaching in from outside that window is missing (as in a
Dispatch run). Hence nprocs=1 here.

References ``tests/golden/hop_eos.npz`` and ``hop_opacity.npz``
------------------------------------------------------------------
Purpose: catch unintended changes to HOP's equation of state, continuum and line-list
opacity on a 4 x 6 (T, rho) grid (3500-12000 K, 1e-10 to 1e-7 g/cm^3) over 5160-5170 AA
(vacuum, R = 20000), vmic = 1 km/s. Densities are in 1/cm^3, opacities in 1/cm. There is
no analytic answer for these states, so the reference records agreement with a reviewed
result, not correctness.

``hop_eos`` was generated at df4c20bf, ``hop_opacity`` at 6a71ff283, both with gfortran
15.1.0, ``OPTS=hoppe`` (-O3 -mcpu=apple-m1), from (md5 prefix): input/TS_absdat.dat
5b921bc86ca0, input/abund/abund_magg22 abcd0926adda, input/isotopes_blp d63a83de4230,
input/linelists/master_line.txt 2821ea404b92, and input/molec_data_BC.txt a4e2d3274cda
(``hop_eos``) or cfb2d3e18c35 (``hop_opacity``).

Tolerance rtol = 1e-5. The references are exact on the generating build. A change of the
optimisation flags alone (-O1 to -O3) moved line-opacity cells by float32 roundoff of at
most 5e-6 relative, so a tighter bound would fail on a different build without any change
in the physics. Changes to the physics or input data below 1e-5 go undetected.

Regenerate with ``python -m tests.support.golden hop_eos hop_opacity`` (review the
summary), then ``--write``, and give the reason in the commit message.
"""
import numpy as np
import pytest

from tests.native.hop_reference import (HOP_DLAM, HOP_LIB, HOP_LINELIST, HOP_VMIC,
                                        hop_eos_reference, hop_grid, hop_opacity_reference,
                                        hop_wavelengths)
from tests.support.golden import assert_matches_golden, load_golden

pytestmark = pytest.mark.skipif(
    not HOP_LIB.exists(),
    reason="libhop.so not built - run 'make hop_standalone'")

RTOL = 1e-5


def test_eos_matches_reference():
    assert_matches_golden(hop_eos_reference(), load_golden("hop_eos"), rtol=RTOL)


def test_opacity_matches_reference():
    assert_matches_golden(hop_opacity_reference(), load_golden("hop_opacity"), rtol=RTOL)


def test_chunking_is_exact_without_a_linelist():
    """Continuum and sampled-metal opacity do not depend on the wavelength grid.

    Each wavelength is computed independently there, so slicing the axis is exact.  This
    is the half of the chunking behaviour that a threaded hop_opacity may rely on.
    """
    from tumult.hop import HopEOS

    temp, rho = hop_grid()
    lam = hop_wavelengths()
    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho)
        for kw in ({}, dict(include_absmet=True)):
            chi1, scat1 = hop.opacity(lam, nchunks=1, **kw)
            chi3, scat3 = hop.opacity(lam, nchunks=3, **kw)
            np.testing.assert_array_equal(chi3, chi1)
            np.testing.assert_array_equal(scat3, scat1)


@pytest.mark.skipif(not HOP_LINELIST.exists(), reason="master line list not present")
def test_chunking_perturbs_line_opacity():
    """Line opacity DOES depend on the wavelength grid it is handed - pinned here.

    ``f0 = find_closest(sp%lam, lam0)`` (bbtrans_mod.f90:614) is the grid index nearest
    the line centre, and both the skip test and the analytic profile extent are keyed on
    the continuum *at that index* (linelist_mod.f90:373,389,419)::

        if (arg <= logref(f0,i,j,k)) cycle
        ref  = tol * cont(f0,i,j,k) / bi
        vm2  = max(real(arg - logref(f0,i,j,k)), a*Rij*over_pi/ref)

    A narrower chunk puts f0 on a different wavelength, so ``ref`` shifts and every line
    is truncated somewhere slightly different.  The effect sits at the ``tolerance``
    level, and widening dlam does not remove it (it is not a line-selection effect).

    Consequence for the threading work: the wavelength axis must NOT be split across
    threads, or results would depend on OMP_NUM_THREADS.  Thread the cell loop instead -
    cont/logref are indexed (f0,i,j,k), so per-cell work is untouched by that split.
    """
    from tumult.hop import HopEOS

    temp, rho = hop_grid()
    lam = hop_wavelengths()
    with HopEOS(feh=0.0) as hop:
        hop.load_state(temp, rho=rho)
        kw = dict(linelist=str(HOP_LINELIST), dlam=HOP_DLAM, vmic=HOP_VMIC,
                  line_scattering=True)
        chi1, _ = hop.opacity(lam, nchunks=1, **kw)
        chi3, _ = hop.opacity(lam, nchunks=3, **kw)
    assert np.any(chi3 != chi1), "chunking artifact vanished - re-read the docstring"
    rel = np.abs(chi3 - chi1) / np.abs(chi1)
    assert rel.max() < 1e-3, f"chunking artifact grew to {rel.max():.2e}"


@pytest.mark.skipif(not HOP_LINELIST.exists(), reason="master line list not present")
def test_linelist_rebuild_is_transparent():
    """Widening the line-list window silently rebuilds and re-solves - same answer.

    A line list can only ever be *added* to a HOP instance, so asking for a wider window
    than was read makes the wrapper free the instance, re-init it from its recipe, replay
    the kept species and re-solve the state (_rebuild/_replay_state).  Nothing about that
    may be visible in the result: reading a wide window directly must give the same
    numbers as narrowing to it after a narrow read.

    This is the contract to preserve if the rebuild ever moves into the C shim.
    """
    from tumult.hop import HopEOS

    temp, rho = hop_grid()
    lam = hop_wavelengths()
    kw = dict(linelist=str(HOP_LINELIST), vmic=HOP_VMIC, line_scattering=True)

    with HopEOS(feh=0.0) as hop:                      # wide window from the start
        hop.load_state(temp, rho=rho, species=["Fe I"])
        direct, _ = hop.opacity(lam, dlam=30.0, **kw)
        n_fe_direct = hop.species_density("Fe I")

    with HopEOS(feh=0.0) as hop:                      # narrow first, then widen
        hop.load_state(temp, rho=rho, species=["Fe I"])
        hop.opacity(lam, dlam=5.0, **kw)
        rebuilt, _ = hop.opacity(lam, dlam=30.0, **kw)
        n_fe_rebuilt = hop.species_density("Fe I")

    np.testing.assert_array_equal(rebuilt, direct)
    np.testing.assert_array_equal(n_fe_rebuilt, n_fe_direct)


def test_reopen_after_close():
    """HOP keeps process-global state; closing must leave it reusable.

    The reference generator and the two reference tests above all open an instance in
    the same process, so this is a precondition for the rest of the module.
    """
    from tumult.hop import HopEOS

    seen = []
    for _ in range(2):
        with HopEOS(feh=0.0) as hop:
            seen.append(hop.load_state([5777.0], rho=[2.5e-7]).ne[0])
    assert seen[0] == seen[1]


def test_option_errors_carry_the_librarys_message():
    """hop_set_option distinguishes an unknown name from an unparseable value, and
    hop_last_error explains which - neither was possible with the old status-1-for-both
    three-setter API."""
    from tumult.hop import HopEOS, HopError

    with HopEOS(feh=0.0) as hop:
        with pytest.raises(HopError, match="unknown option 'nonesuch'"):
            hop.set_option("nonesuch", 1)
        with pytest.raises(HopError, match="cannot read 'maybe' as a logical"):
            hop.set_option("mhd_eos", "maybe")
        # every listed name is settable, and the value crosses as a string either way
        names = hop.option_names()
        assert {"mhd_eos", "ref_lam", "abundances", "voigt_kernel"} <= set(names)
        hop.set_option("ref_lam", 5000)      # int for a real-valued option
        hop.set_option("ref_lam", 5000.0)
        hop.set_option("mhd_eos", False)
        hop.set_option("voigt_kernel", 1)  # inert since the QLG-bot switch; the ABI still carries it


def test_species_error_names_the_species():
    from tumult.hop import HopEOS, HopError

    with HopEOS(feh=0.0) as hop:
        with pytest.raises(HopError, match="unknown species 'Nonium XII'"):
            hop.keep_species("Nonium XII")


def test_second_instance_is_refused():
    """The single-instance guard is what the wrapper's rebuild logic relies on."""
    from tumult.hop import HopEOS, HopError

    with HopEOS(feh=0.0):
        with pytest.raises(HopError, match="already open"):
            HopEOS(feh=0.0)


def test_line_strengths_are_the_ranking_hop_itself_uses():
    """hop_get_linelist hands back the per-line strengths hop_opacity worked out.

    The strength is ``max over cells of (line-centre extinction / continuum extinction)`` --
    the one line ranking that is comparable *between* species, because it carries the
    abundance, the ionization fraction, the dissociation equilibrium and the partition
    function that a gf-and-Boltzmann estimate does not. It is the same quantity a full run
    writes as the trailing column of ``linelist_out.txt``; the point of the accessor is to
    get it from a solved state in memory instead of from a stored file.

    Pinned here by what it must *order*, not by its values: a solar-type photosphere over
    5000-5002 AA is headed by Fe I 5001.87, then Ni I 5000.34, then Ti I 5001.01, and the
    molecular bands in the window sit orders of magnitude below them. Ranked by the naive
    proxy those same bands come out on top, so this is the assertion that would catch the
    accessor silently returning nothing useful.
    """
    from tumult.hop import HopEOS
    from tumult.physics import air2vac

    lam = air2vac(np.arange(5000.0, 5002.0, 0.005))
    with HopEOS(feh=0.0) as hop:
        hop.load_state([5876.0], rho=[2.609e-07])   # tau_Ross ~ 2/3 of a solar MARCS model
        lines = hop.line_strengths(lam, str(HOP_LINELIST), dlam=2.0, vmic=HOP_VMIC)

    assert len(lines) > 0
    assert set(lines.columns) == {"ion", "lam0", "evi", "loggf", "strength"}
    assert (lines.strength >= 0).all()

    inside = lines[(lines.lam0 >= 5000.0) & (lines.lam0 <= 5002.0)]
    top = inside.sort_values("strength", ascending=False)
    assert list(top.ion.head(3)) == ["Fe I", "Ni I", "Ti I"]
    assert abs(top.lam0.iloc[0] - 5001.87) < 0.02

    # the molecules present in this window are far below the atomic lines, which is the
    # whole distinction the proxy cannot make
    molecules = inside[inside.ion.isin(["C2", "CN", "CH", "MgH"])]
    if len(molecules):
        assert molecules.strength.max() < 0.05 * top.strength.iloc[0]


def test_line_strengths_need_a_state():
    """The strengths are a by-product of an opacity call, so there has to be a state."""
    from tumult.hop import HopEOS, HopError

    lam = np.array([5000.0, 5001.0])
    with HopEOS(feh=0.0) as hop:
        with pytest.raises(HopError):
            hop.line_strengths(lam, str(HOP_LINELIST))


def test_linelist_numeric_reader_keeps_supported_and_fallback_syntax(tmp_path):
    """The fast reader covers normal decimals/exponents and preserves list-I/O fallback."""
    from tumult.hop import HopEOS
    from tumult.physics import air2vac

    lines_path = tmp_path / "numeric-forms.txt"
    lines_path.write_text(
        "'Fe I'\n"
        "'0' 3\n"
        "5000.000 0.0 2.0\n"
        "5000.500\t0.0\t2.0D+00 0 1.0E+08 0 1 2\n"
        "5001.000, 0.0, 2.0\n"
    )
    lam = air2vac(np.arange(4999.9, 5001.1, 0.005))
    with HopEOS(feh=0.0) as hop:
        hop.load_state([5876.0], rho=[2.609e-07])
        lines = hop.line_strengths(lam, lines_path, dlam=0.2, vmic=HOP_VMIC)

    inside = lines[(lines.lam0 >= 4999.9) & (lines.lam0 <= 5001.1)]
    np.testing.assert_allclose(inside.lam0, [5000.0, 5000.5, 5001.0], atol=2e-8)
    np.testing.assert_allclose(inside.loggf, 2.0)
