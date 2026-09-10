"""Pin the standalone HOP library (libhop.so) through the tumult.hop wrapper.

These are the gate for the FFI work: the ABI generation, the thicker C shim and the move
of wavelength chunking from the wrapper into OpenMP inside hop_opacity must all leave
these bit-identical.  ``test_chunking_is_neutral`` is the one that matters most for the
threading step - it is the invariant the in-Fortran split has to keep.

Note that chunking is only neutral *within one process*: a spawned worker reads only the
lines within dlam of its own chunk, so a far Lorentz wing reaching in from outside that
window is missing (exactly as in a dispatch run).  Hence nprocs=1 here.
"""
import numpy as np
import pytest

from cases import (HOP_DLAM, HOP_LIB, HOP_LINELIST, HOP_VMIC, hop_eos_cases, hop_grid,
                   hop_opacity_cases, hop_wavelengths)
from conftest import assert_matches_golden, load_golden

pytestmark = pytest.mark.skipif(
    not HOP_LIB.exists(),
    reason="libhop.so not built - run 'make -C HOP/standalone -j8'")


def test_eos_matches_golden():
    assert_matches_golden(hop_eos_cases(), load_golden("hop_eos"))


def test_opacity_matches_golden():
    """Continuum and line opacity, against the golden.

    REGENERATED 2026-09-10, because the previous golden described a stale
    ``HOP/standalone/build/libhop.so`` rather than the committed source. The shipped
    library predated the window line-selection commits (``facc5eb2`` "Keep only the lines
    the windows actually need", ``e49148bb`` "Confine each output patch to its own
    window's lines"), so the suite was green only as long as nobody rebuilt it; rebuilding
    from committed source moves the line terms. Verified by rebuilding the library from an
    untouched checkout of ``HOP/standalone/hop_c_api.f90``, which fails against the old
    golden identically.

    The selection is what moved, and only the line terms: chi_continuum, scat_continuum and
    lam came back bit-identical, chi_linelist shifted by 1.9e-4 max / 2.7e-5 rms over
    452 of 960 cells and scat_linelist by 3.0e-3 max / 2.2e-4 rms over 533. Both move
    *up* in the large majority (392 of 452, 440 of 533) and the median moved cell shifts
    3e-6, so this is a handful of marginal lines entering the set rather than a change to
    any line's physics. The shift also sits where line opacity is comparable to the
    continuum (median line/continuum 1.4 among the moved cells against 1.0 among the
    unmoved), which is where a selection change belongs -- a line that was dropped and is
    now kept matters near its own core, not in a far wing.

    Previously REGENERATED 2026-09-09, for the radiative-damping column of
    ``input/linelists/master_line.txt``. That file used to be written in a convention
    where small column-5 values were multiplied by 1e8 at run time; the code no longer
    does that (``init_linelist2`` reads the column straight into ``lin%ga_in``,
    composition_mod.f90:1119,1129), so the file was rewritten to carry the real s^-1
    values -- 1.160e+00 became 1.160e+08 and so on. The old pairing gave those lines a
    radiative damping ~1e8 too small, i.e. cores and wings that were far too narrow.

    Only the line entries moved: chi_continuum and scat_continuum came back
    bit-identical. chi_linelist shifted by 91% max / 15% rms and scat_linelist by 96%
    max / 25% rms, but the distribution is very skewed -- the median cell moves ~1%,
    and 133 of 960 move by more than 10%. Those are the lines whose column 5 actually
    changed; the shift is concentrated where the line opacity is comparable to the
    continuum (median line/continuum ~1 among them), not in the far wings, which is
    where a damping-constant change belongs. Line opacity goes *up* almost everywhere
    (97% of wing cells), as broader Lorentzian wings should.

    Previously REGENERATED 2026-09-04, when QLG-bot replaced Traving/Eriksson's
    voigt_function2 as the one Voigt kernel (HOP/qlg_mod.f90; the rivals it was chosen
    over are on the `voigt-kernel-zoo` branch): chi_linelist moved 2.9% max / 1e-3 rms
    and scat_linelist 18% max / 6e-3 rms, both in the far wings where line opacity is
    ~1e-13, with the continuum bit-identical. That shift was voigt_function2's error
    being removed, not introduced -- scored against the fp64 Humlicek-Weideman
    reference over a in [1e-4,10] x v in [1e-2,1e3], voigt_function2 reaches 3.0e-1
    relative and is worse than 1% at 722 of 24k sampled points, while qlg_bot_phi peaks
    at 6.1e-3 and is never worse than 1%.
    """
    assert_matches_golden(hop_opacity_cases(), load_golden("hop_opacity"))


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

    gen_goldens.py and the two golden tests above all open an instance in the same
    process, so this is a precondition for the rest of the module.
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
