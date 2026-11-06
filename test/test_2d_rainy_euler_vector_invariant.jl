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
end

@testitem "Rainy Euler vector invariant 2D: elixir_vector_invariant_well_balanced_curvilinear with flux_etec" setup=[
    Setup,
    RainyEulerVectorInvariant2D
] tags=[:rainy_euler_vector_invariant] begin
    # With the ETEC volume flux, both hydrostatic states (constant virtual temperature and
    # constant virtual potential temperature) are preserved up to round-off
    for initial_condition in (initial_condition_isothermal, initial_condition_adiabatic)
        trixi_include(@__MODULE__,
                      joinpath(EXAMPLES_DIR, "tests",
                               "elixir_vector_invariant_well_balanced_curvilinear.jl"),
                      initial_condition = initial_condition, volume_flux = flux_etec,
                      tspan = (0.0, 1.0))
        max_velocity, max_density_change = velocity_and_density_change(sol, semi)
        @test max_velocity < 1.0e-9
        @test max_density_change < 1.0e-11
    end
    # Ensure that we do not have excessive memory allocations
    # (e.g., from type instabilities)
    @test_allocations(TrixiAtmo.Trixi.rhs_hyperbolic!, semi, sol, 100)
end
