@testsnippet RainyEulerVectorInvariant2D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "euler", "precipitation")

    # Maximum absolute value of the velocity components and maximum relative change of the
    # density between the initial and the final solution
    function velocity_and_density_change(sol, semi)
        u_0 = Trixi.wrap_array(sol.u[1], semi)
        u_end = Trixi.wrap_array(sol.u[end], semi)
        max_velocity = maximum(abs, view(u_end, 2:3, :, :, :))
        max_density_change = maximum(abs,
                                     view(u_end, 1, :, :, :) .- view(u_0, 1, :, :, :)) /
                             maximum(view(u_0, 1, :, :, :))
        return max_velocity, max_density_change
    end

    # Maximum of ρ_v / ρ_vs(T) - 1 over all nodes
    function maximum_supersaturation(u, equations, solver, cache)
        max_supersaturation = -Inf
        for element in Trixi.eachelement(solver, cache), j in Trixi.eachnode(solver),
            i in Trixi.eachnode(solver)

            u_node = Trixi.get_node_vars(u, equations, solver, i, j, element)
            T = TrixiAtmo.temperature(u_node, equations)
            rho_vs = TrixiAtmo.saturation_vapor_pressure(T, equations) /
                     (equations.R_v * T)
            max_supersaturation = max(max_supersaturation, u_node[5] / rho_vs - 1)
        end
        return max_supersaturation
    end
end

@testitem "Rainy Euler vector invariant 2D: elixir_vector_invariant_well_balanced_curvilinear isothermal" setup=[
    Setup,
    RainyEulerVectorInvariant2D
] tags=[:rainy_euler_vector_invariant] begin
    # The hydrostatic state with constant virtual temperature is preserved up to round-off
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "tests",
                           "elixir_vector_invariant_well_balanced_curvilinear.jl"),
                  initial_condition = initial_condition_isothermal,
                  nonconservative_flux = flux_nonconservative_isothermal_isentropic,
                  tspan = (0.0, 1.0))
    max_velocity, max_density_change = velocity_and_density_change(sol, semi)
    @test max_velocity < 1.0e-9
    @test max_density_change < 1.0e-11
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Rainy Euler vector invariant 2D: elixir_vector_invariant_well_balanced_curvilinear adiabatic" setup=[
    Setup,
    RainyEulerVectorInvariant2D
] tags=[:rainy_euler_vector_invariant] begin
    # The hydrostatic state with constant virtual potential temperature is preserved up to
    # round-off by both nonconservative fluxes
    for nonconservative_flux in (flux_nonconservative_isothermal_isentropic,
                                 flux_nonconservative_isentropic)
        trixi_include(@__MODULE__,
                      joinpath(EXAMPLES_DIR, "tests",
                               "elixir_vector_invariant_well_balanced_curvilinear.jl"),
                      initial_condition = initial_condition_adiabatic,
                      nonconservative_flux = nonconservative_flux,
                      tspan = (0.0, 1.0))
        max_velocity, max_density_change = velocity_and_density_change(sol, semi)
        @test max_velocity < 1.0e-9
        @test max_density_change < 1.0e-11
    end
end

@testitem "Rainy Euler vector invariant 2D: elixir_vector_invariant_well_balanced_curvilinear isothermal with isentropic flux" setup=[
    Setup,
    RainyEulerVectorInvariant2D
] tags=[:rainy_euler_vector_invariant] begin
    # The arithmetic mean of the virtual potential temperature in
    # `flux_nonconservative_isentropic` does not preserve the hydrostatic state with
    # constant virtual temperature, which generates spurious velocities
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "tests",
                           "elixir_vector_invariant_well_balanced_curvilinear.jl"),
                  initial_condition = initial_condition_isothermal,
                  nonconservative_flux = flux_nonconservative_isentropic,
                  tspan = (0.0, 1.0))
    max_velocity, _ = velocity_and_density_change(sol, semi)
    @test max_velocity > 1.0e-7
end

@testitem "Rainy Euler vector invariant 2D: elixir_vector_invariant_moist_bubble" setup=[
    Setup,
    RainyEulerVectorInvariant2D
] tags=[:rainy_euler_vector_invariant] begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "buoyancy",
                           "elixir_vector_invariant_moist_bubble.jl"),
                  trees_per_dimension = (16, 8), tspan = (0.0, 20.0))
    u_end = Trixi.wrap_array(sol.u[end], semi)
    @test all(isfinite, u_end)

    # The bubble rises
    @test maximum(view(u_end, 3, :, :, :)) > 0

    # Total mass and total water are conserved, since no rain reaches the bottom
    total_mass(u, equations) = u[1]
    total_water(u, equations) = u[5] + u[6] + u[7]
    for func in (total_mass, total_water)
        @test isapprox(Trixi.integrate(func, sol.u[end], semi),
                       Trixi.integrate(func, sol.u[1], semi); rtol = 1.0e-12)
    end

    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end

@testitem "Rainy Euler vector invariant 2D: elixir_vector_invariant_hoeck_bubble" setup=[
    Setup,
    RainyEulerVectorInvariant2D
] tags=[:rainy_euler_vector_invariant] begin
    trixi_include(@__MODULE__,
                  joinpath(EXAMPLES_DIR, "buoyancy",
                           "elixir_vector_invariant_hoeck_bubble.jl"),
                  trees_per_dimension = (16, 16), tspan = (0.0, 10.0))
    u_end = Trixi.wrap_array(sol.u[end], semi)
    @test all(isfinite, u_end)

    # After the saturation adjustment, the air is not supersaturated, and the cloud water
    # and rain densities are nonnegative
    @test maximum_supersaturation(u_end, equations, solver, semi.cache) < 1.0e-8
    @test minimum(view(u_end, 6, :, :, :)) > -1.0e-12
    @test minimum(view(u_end, 7, :, :, :)) > -1.0e-12

    # The total mass is conserved, since no rain reaches the bottom
    total_mass(u, equations) = u[1]
    @test isapprox(Trixi.integrate(total_mass, sol.u[end], semi),
                   Trixi.integrate(total_mass, sol.u[1], semi); rtol = 1.0e-12)

    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
