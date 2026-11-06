# Rising moist bubble with rain of `elixir_hoeck_bubble.jl` for the compressible rainy Euler
# equations in vector-invariant form with the virtual potential temperature, to compare
# with the `CompressibleRainyEulerEquations2D`.
#
# The setup (domain, background state, humidity bubble, and final time) and the physics are
# the same as in `elixir_hoeck_bubble.jl`, except for the constant-kappa thermodynamics:
# - instantaneous saturation adjustment after every Runge-Kutta stage with
#   `NonlinearSolveDG` instead of a condensation source term,
# - the rates of autoconversion, accretion, and rain evaporation of the
#   `CompressibleRainyEulerEquations2D`, with evaporation only (no condensation on rain),
# - the terminal velocity of rain of the `CompressibleRainyEulerEquations2D`
#   (`rain_fall_speed = :doppler_etal`),
# - no frictional heating by falling rain (`rain_frictional_heating = false`), which is not
#   contained in the energy equation of the `CompressibleRainyEulerEquations2D`.
# The hydrostatic background state and the bubble are computed as in
# `elixir_hoeck_bubble.jl`, but with the saturation vapor pressure of the constant-kappa
# system. The remaining difference is the discretization (`P4estMesh`, standard fluxes).
#
# References:
# - Wojciech W. Grabowski, Terry L. Clark (1991)
#   Cloud-environment interface instability: Rising thermal calculations in two spatial
#   dimensions
#   https://doi.org/10.1175/1520-0469(1991)048<0527:CIIRTC>2.0.CO;2
# - Sabine Doppler, Philip L. Lederer, Joachim Schöberl, Henry von Wahl (2024)
#   A discontinuous Galerkin approach for atmospheric flows with implicit condensation
#   https://doi.org/10.1016/j.jcp.2023.112713
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiAtmo
using TrixiAtmo: saturation_vapor_pressure, saturation_residual,
                 saturation_residual_jacobian
using NLsolve: nlsolve

# domain
coordinates_min = (0.0, 0.0)
coordinates_max = (2400.0, 2400.0)

equations = CompressibleRainyEulerVectorInvariantEquations2D(c_pd = 1004, c_vd = 717,
                                                             R_v = 461,
                                                             gravity = EARTH_GRAVITATIONAL_ACCELERATION,
                                                             rain_fall_speed = :doppler_etal,
                                                             rain_frictional_heating = false)

###############################################################################
# initial condition, computed as in `elixir_hoeck_bubble.jl`

# hydrostatic dry potential temperature of the background state
function theta_dry_background(z, equations)
    @unpack c_pd, R_d, p_0 = equations
    surface_temperature = 283
    surface_pressure = 85_000
    stratification = 1.3e-5

    theta_0 = surface_temperature * (p_0 / surface_pressure)^(R_d / c_pd)
    return theta_0 * exp(stratification * z)
end

# residual of the relative humidity as defined for the `CompressibleRainyEulerEquations2D`
function humidity_residual(rho_dry, rho_vapor, rho_vs, humidity, equations)
    eps_ = equations.R_d / equations.R_v
    return (rho_vs * (rho_dry + rho_vapor / eps_) * humidity -
            rho_vapor * (rho_dry + rho_vs / eps_))
end

# residual of the discrete hydrostatic balance between two layers with distance `dz`
function hydrostatic_residual!(residual, guess, pressure_lower, humidity, z, dz,
                               equations)
    @unpack c_pd, R_d, R_v, p_0, g = equations
    pressure, rho_dry, rho_vapor, temperature = guess
    rho_vs = saturation_vapor_pressure(temperature, equations) / (R_v * temperature)

    residual[1] = (pressure - pressure_lower) / dz + (rho_dry + rho_vapor) * g
    residual[2] = pressure - temperature * (rho_dry * R_d + rho_vapor * R_v)
    residual[3] = (theta_dry_background(z, equations) -
                   temperature * (p_0 / pressure)^(R_d / c_pd))
    residual[4] = 1000 * humidity_residual(rho_dry, rho_vapor, rho_vs, humidity, equations)
    return nothing
end

# layers of the hydrostatic background state with pressure, dry air density, vapor density,
# and temperature, as for `AtmosphereLayersRainyBubble`
function hydrostatic_layers(equations; total_height, dz = 1.0)
    humidity_background = 0.2
    surface_pressure = 85_000.0

    n = convert(Int, total_height / dz)
    layer_data = zeros(n + 1, 4)

    # (slightly above) surface layer
    z = 0.01
    function residual_surface!(residual, guess)
        return hydrostatic_residual!(residual, guess, surface_pressure,
                                     humidity_background, z, 0.01, equations)
    end
    layer_data[1, :] .= nlsolve(residual_surface!,
                                [surface_pressure, 1.4, 0.04, 300.0]).zero

    # iterate up the atmosphere
    for i in 1:n
        z += dz
        pressure_lower = layer_data[i, 1]
        function residual_layer!(residual, guess)
            return hydrostatic_residual!(residual, guess, pressure_lower,
                                         humidity_background, z, dz, equations)
        end
        layer_data[i + 1, :] .= nlsolve(residual_layer!, layer_data[i, :], ftol = 1e-10,
                                        iterations = 20).zero
    end

    return layer_data
end

total_height = coordinates_max[2] + 1
layer_dz = 1.0
layer_data = hydrostatic_layers(equations; total_height, dz = layer_dz)

function initial_condition_hoeck_bubble(x, t,
                                        equations::CompressibleRainyEulerVectorInvariantEquations2D)
    @unpack c_pd, R_d, R_v, p_0, g = equations

    # humidity bubble: saturated core and cos² transition to the background humidity
    humidity_background = 0.2
    humidity_max = 1
    radius_outer, radius_inner = 300, 200
    x_center, z_center = 1200, 800

    r = sqrt((x[1] - x_center)^2 + (x[2] - z_center)^2)
    if r > radius_outer
        humidity = humidity_background
    elseif r > radius_inner
        humidity = (humidity_background +
                    (humidity_max - humidity_background) *
                    cos(pi * (r - radius_inner) / (2 * (radius_outer - radius_inner)))^2)
    else
        humidity = humidity_max
    end

    # hydrostatic pressure interpolated linearly between the layers
    dz = layer_dz
    z = x[2]
    n = convert(Int, floor((z + eps()) / dz)) + 1
    z_lower = (n - 1) * dz
    z_upper = n * dz
    if z_lower == total_height
        z_upper = z_lower + dz
        n = n - 1
    end
    if n == 0
        n = 1
    end
    pressure = ((layer_data[n + 1, 1] * (z - z_lower) +
                 layer_data[n, 1] * (z_upper - z)) / dz)

    # dry air density, vapor density, and temperature at the hydrostatic pressure with the
    # dry potential temperature of the background state and the given humidity
    function perturbation_residual!(residual, guess)
        rho_dry, rho_vapor, temperature = guess
        rho_vs = saturation_vapor_pressure(temperature, equations) / (R_v * temperature)
        residual[1] = 30 * humidity_residual(rho_dry, rho_vapor, rho_vs, humidity,
                              equations)
        residual[2] = (theta_dry_background(z, equations) -
                       temperature * (p_0 / pressure)^(R_d / c_pd))
        residual[3] = pressure - (rho_dry * R_d + rho_vapor * R_v) * temperature
        return nothing
    end
    rho_dry, rho_vapor, _ = nlsolve(perturbation_residual!, layer_data[n, 2:4],
                                    ftol = 1e-9, iterations = 20).zero

    # no cloud water and rain
    return prim2cons(SVector(rho_dry + rho_vapor, 0, 0, pressure, rho_vapor, 0, 0, g * z),
                     equations)
end

initial_condition = initial_condition_hoeck_bubble

###############################################################################
# semidiscretization of the compressible rainy Euler equations in vector-invariant form

boundary_conditions = (; y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

polydeg = 3
surface_flux = (flux_lax_friedrichs, flux_nonconservative_isothermal_isentropic)
volume_flux = (flux_central, flux_nonconservative_isothermal_isentropic)
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

trees_per_dimension = (64, 64)
mesh = P4estMesh(trees_per_dimension; polydeg = polydeg,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = (true, false), initial_refinement_level = 0)

# Condensation by saturation adjustment (see `stage_limiter!` below) and rain evaporation
# only, as for the `CompressibleRainyEulerEquations2D`
source_terms = SourceTermsWarmRain(condensation_time_scale = nothing,
                                   condensation_on_rain = false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 600.0)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     output_directory = "out",
                                     solution_variables = cons2prim)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution)

# saturation adjustment of (rho_theta_v, rho_vapor, rho_cloud) after every stage
stage_limiter! = NonlinearSolveDG(saturation_residual, saturation_residual_jacobian,
                                  SVector(4, 5, 6))

###############################################################################
# run the simulation
sol = solve(ode, SSPRK43(; stage_limiter! = stage_limiter!); ode_default_options()...,
            maxiters = 1.0e7, save_everystep = false, callback = callbacks);
