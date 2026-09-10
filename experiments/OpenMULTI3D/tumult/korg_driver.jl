# Synthesise one window with Korg.jl, driven by the keyword file tumult/korg.py writes.
#
#     julia --project=<KORG_PATH> korg_driver.jl <input file>
#
# The input file holds one "KEY value" per line; a key may repeat, which is how lists
# (MU_VALUES) and per-element abundances (ABUNDANCE <symbol> <A(X)>) arrive. Keeping the
# format this plain avoids adding a JSON dependency to the Korg environment.
#
# Wavelengths are air on the tumult side and vacuum inside Korg. The window bounds and
# every line wavelength are converted once here, and the spectrum is written back on the
# air axis, so sample j is at air LAMBDA_MIN + j*LAMBDA_STEP. Do not reach for
# `air_wavelengths=true` instead: it converts the requested range and then reports
# `sol.wavelengths` converted a second time, putting line cores 1.4 A off at 5000 A.
#
# The output mirrors the Turbospectrum result file so tumult can read both with one
# parser: a "# mu-points" header, then columns
#
#     lambda_air  F/Fc  F  [I(mu1)/Ic  I(mu1)  I(mu2)/Ic  I(mu2)  ...]
#
# Intensity columns appear only when MU_VALUES is given, and cost a second synthesis:
# Korg returns a continuum flux but not a continuum intensity, so the normalisation is
# taken from a continuum-only pass with an empty line list.

using Printf, CSV, DataFrames, Korg

# --- input file ------------------------------------------------------------------

"Parse the keyword file into key => vector of raw value strings."
function read_input(path)
    input = Dict{String,Vector{String}}()
    for text in eachline(path)
        line = strip(text)
        (isempty(line) || startswith(line, '#')) && continue
        key, rest = split(line, limit=2)
        push!(get!(input, String(key), String[]), String(strip(rest)))
    end
    input
end

has(input, key) = haskey(input, key) && !isempty(input[key])
one_of(input, key) = input[key][1]
one_of(input, key, default) = has(input, key) ? one_of(input, key) : default
num(input, key, default::Float64) =
    has(input, key) ? parse(Float64, one_of(input, key)) : default
flag(input, key, default::Bool) = has(input, key) ? one_of(input, key) == "true" : default

"ABUNDANCE lines as the Dict format_A_X takes, keyed by element symbol."
function read_abundances(input)
    abundances = Dict{String,Float64}()
    for entry in get(input, "ABUNDANCE", String[])
        symbol, value = split(entry)
        abundances[String(symbol)] = parse(Float64, value)
    end
    abundances
end

"MU_VALUES as the mu grid to report intensities on, or nothing for flux only."
function read_mu_values(input)
    has(input, "MU_VALUES") || return nothing
    values = parse.(Float64, input["MU_VALUES"])
    length(values) == 1 && isinteger(values[1]) ? Int(values[1]) : values
end

# --- line list -------------------------------------------------------------------

"The CSV tumult/korg.py writes, as Korg Lines. Wavelengths are air in the file."
function read_linelist(path)
    df = CSV.read(path, DataFrame)
    [Korg.Line(Korg.air_to_vacuum(r.lam0), r.loggf, Korg.Species(r.ion), r.evi,
               r.gamma_rad, r.gamma_stark, r.vdw) for r in eachrow(df)]
end

# --- output ----------------------------------------------------------------------

function write_spectrum(path, lam_air, flux, cntm, intensity, cntm_intensity, mus)
    open(path, "w") do fh
        @printf(fh, "# mu-points %s\n", join((@sprintf("%.6f", m) for m in mus), " "))
        for j in eachindex(lam_air)
            @printf(fh, "%12.5f %14.6e %14.6e", lam_air[j], flux[j]/cntm[j], flux[j])
            for i in eachindex(mus)
                @printf(fh, " %14.6e %14.6e",
                        intensity[i, j]/cntm_intensity[i, j], intensity[i, j])
            end
            @printf(fh, "\n")
        end
    end
end

"One row per layer, prefixed so tumult can pick it out of stdout."
function write_atmosphere(atm)
    for layer in atm.layers
        @printf("tau,T,ne,nt,z %14.6e %10.2f %14.6e %14.6e %14.6e\n",
                layer.tau_ref, layer.temp, layer.electron_number_density,
                layer.number_density, layer.z)
    end
end

# --- main ------------------------------------------------------------------------

function main(args)
    length(args) == 1 || error("usage: korg_driver.jl <input file>")
    input = read_input(args[1])

    atm = Korg.read_model_atmosphere(one_of(input, "ATMOSPHERE"))
    linelist = has(input, "LINELIST") ? read_linelist(one_of(input, "LINELIST")) : Korg.Line[]
    A_X = Korg.format_A_X(num(input, "METALS_H", 0.0), num(input, "ALPHA_H", 0.0),
                          read_abundances(input); solar_relative=false)

    # NPOINTS wins over LAMBDA_STEP: a run that set a resolution rather than a step has
    # a constant-velocity grid, and matching its sample count over the same end points
    # is the closest a linearly spaced Korg grid can come to it.
    blue, red = num(input, "LAMBDA_MIN", 5000.0), num(input, "LAMBDA_MAX", 5002.0)
    npoints = if has(input, "NPOINTS")
        parse(Int, one_of(input, "NPOINTS"))
    else
        round(Int, (red - blue)/num(input, "LAMBDA_STEP", 0.01)) + 1
    end
    lam_air = range(blue, red; length=npoints)
    wls = Korg.Wavelengths(range(Korg.air_to_vacuum(first(lam_air)),
                                 Korg.air_to_vacuum(last(lam_air)); length=npoints))

    mu_values = read_mu_values(input)
    kwargs = (vmic=num(input, "VMIC", 1.0),
              line_buffer=num(input, "LINE_BUFFER", 10.0),
              hydrogen_lines=flag(input, "HYDROGEN_LINES", true),
              line_cutoff_threshold=num(input, "LINE_CUTOFF_THRESHOLD", 3e-4),
              tau_scheme=one_of(input, "TAU_SCHEME", "anchored"))

    if isnothing(mu_values)
        sol = Korg.synthesize(atm, linelist, A_X, wls; kwargs...)
        mus, intensity, cntm_intensity = Float64[], nothing, nothing
    else
        scheme = one_of(input, "I_SCHEME", "linear")
        sol = Korg.synthesize(atm, linelist, A_X, wls; mu_values=mu_values,
                              I_scheme=scheme, kwargs...)
        # Korg returns intensity for every layer; layer 1 is the surface, which is what
        # radiative_transfer itself integrates to get the flux
        intensity = sol.intensity[:, :, 1]
        mus = first.(sol.mu_grid)
        cont = Korg.synthesize(atm, Korg.Line[], A_X, wls; mu_values=mu_values,
                               I_scheme=scheme, kwargs...)
        cntm_intensity = cont.intensity[:, :, 1]
    end

    write_spectrum(one_of(input, "RESULTFILE"), lam_air, sol.flux, sol.cntm,
                   intensity, cntm_intensity, mus)
    write_atmosphere(atm)
    @printf("korg_driver: %d layers, %d lines, %d wavelengths, air %.3f - %.3f\n",
            length(atm.layers), length(linelist), npoints, blue, red)
end

main(ARGS)
