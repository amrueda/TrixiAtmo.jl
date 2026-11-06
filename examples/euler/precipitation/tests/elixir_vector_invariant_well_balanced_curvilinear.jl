# Well-balanced test for the moist compressible Euler equations with warm-rain microphysics
# in vector-invariant form with the virtual potential temperature on a curvilinear mesh.
# The atmosphere is at rest in hydrostatic balance with constant virtual temperature
# (`initial_condition_isothermal`) or constant virtual potential temperature
# (`initial_condition_adiabatic`) and is subsaturated, so that the microphysics is inactive.
#
# `flux_nonconservative_isothermal` preserves both hydrostatic states up to round-off
# errors, whereas `flux_nonconservative_isentropic` preserves only the state with constant
# virtual potential temperature.
#
# References:
# - Marco Artiano, Kieran Ricardo, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2026)
#   Entropy-Stable and Well-Balanced Discontinuous Galerkin Methods for the Compressible
#   Euler Equations in Vector-Invariant Form
#   https://arxiv.org/abs/2610.05820 (pre-print)
using OrdinaryDiffEqSSPRK
using Trixi, TrixiAtmo

equations = CompressibleRainyEulerVectorInvariantEquations2D(c_pd = 1004, c_vd = 717,
                                                             R_v = 461,
                                                             gravity = EARTH_GRAVITATIONAL_ACCELERATION)

initial_condition = initial_condition_isothermal
nonconservative_flux = flux_nonconservative_isothermal

polydeg = 2
basis = LobattoLegendreBasis(polydeg)
surface_flux = (flux_lax_friedrichs, nonconservative_flux)
volume_flux = (flux_central, nonconservative_flux)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

function mapping(xi, eta)
    x = xi + 0.1 * sinpi(xi) * sinpi(eta)
    y = eta + 0.1 * sinpi(xi) * sinpi(eta)
    return SVector(1000 * 0.5 * (1 + x), 1000 * 0.5 * (1 + y))
end
trees_per_dimension = (16, 16)

mesh = P4estMesh(trees_per_dimension, polydeg = polydeg,
                 mapping = mapping, periodicity = (false, false),
                 initial_refinement_level = 0)

boundary_conditions = (; x_pos = boundary_condition_slip_wall,
                       x_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall,
                       y_neg = boundary_condition_slip_wall)

# The subsaturated atmosphere without cloud water and rain is not affected by the
# microphysics, which is included to check that it does not disturb the balance.
source_terms = SourceTermsWarmRain(condensation_time_scale = 5.0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)

dt = 0.01
tspan = (0.0, 2.0)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

ode = semidiscretize(semi, tspan)

sol = solve(ode, SSPRK43(); dt = dt, ode_default_options()..., callback = callbacks,
            adaptive = false)
