"""The scripting API as used by the benchmark calculations of docs/API_benchmarks.txt.

Each function is a complete workflow in the calls a user writes. ``SynthRequest`` contains
the state needed to produce synthetic spectra: a request with departures produces an NLTE
spectrum, one without departures an LTE spectrum. Observed or synthetic targets are inputs
to analysis operations and are not request state.

The arguments are the public input types of :mod:`tumult.synthesis` (``SynthConfig``,
``Star``, ``Atmosphere``, ``AtmosphereGrid``), wavelength ranges, paths of model atoms,
spectra and line masks, or objects returned by an earlier call.
``run_reference_benchmarks.py`` binds the reference assets to them.
"""

from tumult import abundance, departures, sparams, synthesis
from tumult.synthesis import SynthRequest


def single_lte(config, atmosphere_grid, star, windows):
    """Synthesize an LTE spectrum with automatic transition pruning."""
    request = SynthRequest(
        config=config,
        atmosphere=atmosphere_grid,
        star=star,
        windows=windows,
    )
    return synthesis.synthesize(request)


def departure_solve(config, atmosphere, star, atom):
    """Optimize the atom wavelength grid and solve departures on one atmosphere."""
    request = SynthRequest(config=config, atmosphere=atmosphere, star=star)
    return departures.solve(request, atom=atom)


def single_precomputed_nlte(config, atmosphere, star, solved_departures, windows):
    """Synthesize an NLTE spectrum from previously solved departures."""
    request = SynthRequest(
        config=config,
        atmosphere=atmosphere,
        star=star,
        departures=solved_departures,
        windows=windows,
    )
    return synthesis.synthesize(request)


def single_ondemand_nlte(config, atmosphere, star, atom, windows):
    """Synthesize an NLTE spectrum, solving its departures on demand."""
    request = SynthRequest(
        config=config,
        atmosphere=atmosphere,
        star=star,
        windows=windows,
    )
    request = request.update(departures=departures.request(request, atom=atom))
    return synthesis.synthesize(request)


def nlte_correction(config, atmosphere_grid, star, atom, **kwargs):
    """Calculate equilibrium NLTE abundance corrections for the atom's transitions."""
    request = SynthRequest(
        config=config,
        atmosphere=atmosphere_grid,
        star=star,
    )
    request = request.update(departures=departures.request(request, atom=atom))
    return abundance.nlte_corrections(
        request,
        method="cog",
        background_lines=False,
        **kwargs,
    )


def abundance_fit(
    config,
    target,
    line_mask,
    atmosphere_grid,
    star,
    **kwargs,
):
    """Fit LTE abundances to an observed target with adaptive curves of growth."""
    request = SynthRequest(
        config=config,
        atmosphere=atmosphere_grid,
        star=star,
    )
    return abundance.fit(
        request,
        target=target,
        line_mask=line_mask,
        method="cog",
        background_lines=False,
        **kwargs,
    )


def ionisation_balance(
    config,
    target,
    line_mask,
    atmosphere_grid,
    initial_star,
    **kwargs,
):
    """Infer stellar parameters from a full four-diagnostic ionisation balance."""
    request = SynthRequest(
        config=config,
        atmosphere=atmosphere_grid,
        star=initial_star,
    )
    return sparams.ionbalance(
        request,
        target=target,
        line_mask=line_mask,
        method="cog",
        background_lines=False,
        **kwargs,
    )


def stellar_analysis(config, atmosphere_grid, initial_star, target, fe_lmask, atom, lmask,
                     windows):
    """Infer a star and its abundances, then synthesize the requested windows."""
    request = SynthRequest(config=config, atmosphere=atmosphere_grid, star=initial_star)

    star = sparams.ionbalance(
        request,
        target=target,
        line_mask=fe_lmask,
        method="cog",
        background_lines=False,
    )

    deps = departures.request(request, atom=atom)
    request = request.update(star=star, departures=deps)

    abnd = abundance.fit(
        request,
        target=target,
        line_mask=lmask,
        method="cog",
        background_lines=False,
    )

    request = request.update(windows=windows, abundances=abnd)
    synth = synthesis.synthesize(request)
    return star, abnd, synth
