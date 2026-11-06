# Rising moist bubble for the compressible Euler equations with warm-rain microphysics in
# vector-invariant form with the virtual potential temperature as prognostic variable.
#
# The background state has constant virtual potential temperature in hydrostatic balance
# and a relative humidity of 95%. A warm bubble in the virtual potential temperature at the
# same pressure rises, cools, and condenses; the cloud water converts to rain, which falls.
# Condensation is modeled as relaxation towards saturation, see `SourceTermsWarmRain`.
using OrdinaryDiffEqLowStorageRK
using Trixi, TrixiAtmo
using TrixiAtmo: saturation_vapor_pressure

equations = CompressibleRainyEulerVectorInvariantEquations2D(c_pd = 1004, c_vd = 717,
                                                             R_v = 461,
                                                             gravity = EARTH_GRAVITATIONAL_ACCELERATION)

function initial_condition_moist_bubble(x, t,
                                        equations::CompressibleRainyEulerVectorInvariantEquations2D)
    @unpack p_0, c_pd, R_d, R_v, g = equations
    z = x[2]

    # Hydrostatic background with constant virtual potential temperature
    theta_v_0 = 300
    exner = 1 - g * z / (c_pd * theta_v_0)
    p = p_0 * exner^(c_pd / R_d)

    # Warm bubble at the same pressure, which changes only the density
    x_c = 10_000
    z_c = 2_000
    radius = 2_000
    r = sqrt((x[1] - x_c)^2 + (z - z_c)^2) / radius
    delta_theta_v = r < 1 ? 2 * cospi(0.5 * r)^2 : zero(r)
    theta_v = theta_v_0 + delta_theta_v
    rho = p / (R_d * theta_v * exner)

    # Water vapor with prescribed relative humidity. The temperature
    # T = p / (R_d rho + (R_v - R_d) rho_vapor) depends weakly on the vapor density, so a
    # few fixed-point iterations suffice.
    relative_humidity = 0.95
    rho_vapor = zero(rho)
    for _ in 1:10
        T = p / (R_d * rho + (R_v - R_d) * rho_vapor)
        rho_vapor = relative_humidity * saturation_vapor_pressure(T, equations) /
                    (R_v * T)
    end

    return prim2cons(SVector(rho, 0, 0, p, rho_vapor, 0, 0, g * z), equations)
end

initial_condition = initial_condition_moist_bubble

###############################################################################
# semidiscretization

# Standard fluxes: central volume flux, local Lax-Friedrichs surface flux, and
# nonconservative terms preserving hydrostatic states with constant virtual temperature and
# constant virtual potential temperature.
polydeg = 3
surface_flux = (flux_lax_friedrichs, flux_nonconservative_isothermal)
volume_flux = (flux_central, flux_nonconservative_isothermal)
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (0.0, 0.0)
coordinates_max = (20_000.0, 10_000.0)
trees_per_dimension = (32, 16)
mesh = P4estMesh(trees_per_dimension; polydeg = polydeg,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = (true, false), initial_refinement_level = 0)

boundary_conditions = (; y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

source_terms = SourceTermsWarmRain(condensation_time_scale = 5.0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1000.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        save_solution, stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            # solve needs some value here, but it is overwritten by the stepsize_callback
            dt = 1,
            ode_default_options()..., callback = callbacks)
