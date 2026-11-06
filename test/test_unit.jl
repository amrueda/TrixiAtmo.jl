@testitem "Unit: Consistency check for EC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 330.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_ec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for TEC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 330.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_tec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for ETEC flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 330.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3])
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_etec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for LMARS flux with Potential Temperature: CEPTE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquations2D(c_p = 1004.0, c_v = 717.0)
    flux_lmars = FluxLMARS(340)
    u = SVector(1.1, -0.5, 2.34, 330.0)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4])
    u_2d = SVector(u[1], u[2], 0, u[4])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquations1D(c_p = 1004.0,
                                                                    c_v = 717.0)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_lmars(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquations3D(c_p = 1004.0, c_v = 717.0)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], -u[2], 0.0, 0.0, u[5])
    u_1d = SVector(u[1], -u[2], u[5])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_lmars(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5]]
end

@testitem "Unit: Consistency check for EC flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], u[5])
    u_2d = SVector(u[1], u[2], 0, u[4], u[5])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3], u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_ec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_ec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5], u[6])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_ec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_ec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5, 6]]
end

@testitem "Unit: Consistency check for FluxTurbo EC flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # The `FluxTurbo` specialization of `(flux_ec, flux_nonconservative_isothermal)`
    # must return the same left and right contributions as the fused evaluation of the
    # conservative and the nonconservative flux.
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    volume_flux = (flux_ec, flux_nonconservative_isothermal)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5, 1700.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    @test Trixi.nturbovars(volume_flux, equations) == Val(7)
    turbo_ll = Trixi.cons2turbo(volume_flux, u_ll..., equations)
    turbo_rr = Trixi.cons2turbo(volume_flux, u_rr..., equations)
    @test length(turbo_ll) == 7
    # rho, v1, v2, v3, theta, p, phi
    @test SVector(turbo_ll[1], turbo_ll[2], turbo_ll[3], turbo_ll[4], turbo_ll[6],
                  turbo_ll[7]) ≈ cons2prim(u_ll, equations)
    @test turbo_ll[5] ≈ u_ll[5] / u_ll[1]

    for normal_direction in normal_directions
        flux_left, flux_right = Trixi.flux_turbo(volume_flux, turbo_ll..., turbo_rr...,
                                                 normal_direction..., equations)
        flux_cons = flux_ec(u_ll, u_rr, normal_direction, equations)
        flux_noncons = flux_nonconservative_isothermal(u_ll, u_rr,
                                                       normal_direction, equations)
        @test SVector(flux_left) ≈ flux_cons + 0.5 * flux_noncons
        @test SVector(flux_right) ≈ flux_cons - 0.5 * flux_noncons
    end
end

@testitem "Unit: Consistency check for TEC flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], u[5])
    u_2d = SVector(u[1], u[2], 0, u[4], u[5])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3], u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_tec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_tec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5], u[6])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_tec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_tec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5, 6]]
end

@testitem "Unit: Consistency check for ETEC flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], u[5])
    u_2d = SVector(u[1], u[2], 0, u[4], u[5])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # test when u_ll is not the same as u_rr
    u_rr_1d = SVector(2.1, 0.3, 280.5, 1700)
    u_rr_2d = SVector(u_rr_1d[1], u_rr_1d[2], 0.0, u_rr_1d[3], u_rr_1d[4])
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_rr_1d, 1, equations_1d)
    flux_2d = flux_etec(u_2d, u_rr_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_etec(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5], u[6])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_etec(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_etec(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5, 6]]
end

@testitem "Unit: Consistency check for LMARS flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    flux_lmars = FluxLMARS(340)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1700)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    #check consistency between 1D and 2D EC fluxes
    u_1d = SVector(u[1], u[2], u[4], u[5])
    u_2d = SVector(u[1], u[2], 0, u[4], u[5])
    normal_1d = SVector(-0.3)
    normal_2d = SVector(normal_1d[1], 0.0)
    equations_1d = CompressibleEulerPotentialTemperatureEquationsWithGravity1D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    equations_2d = equations
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    flux_2d = flux_lmars(u_2d, u_2d, normal_2d, equations_2d)
    @test flux_1d ≈ flux(u_1d, normal_1d, equations_1d)
    @test flux_1d ≈ flux_2d[[1, 2, 4, 5]]

    # check consistency for 3D EC flux
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1700)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_lmars(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    # check consistency between 1D and 3D EC fluxes
    u_3d = SVector(u[1], u[2], 0.0, 0.0, u[5], u[6])
    normal_3d = SVector(normal_1d[1], 0.0, 0.0)
    equations_3d = equations
    flux_3d = flux_lmars(u_3d, u_3d, normal_3d, equations_3d)
    flux_1d = normal_1d[1] * flux_lmars(u_1d, u_1d, 1, equations_1d)
    @test flux_1d ≈ flux_3d[[1, 2, 5, 6]]
end

@testitem "Unit: Consistency check for ES LMARS flux with Potential Temperature with gravity: CEPTEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleEulerPotentialTemperatureEquationsWithGravity2D(c_p = 1004.0,
                                                                            c_v = 717.0,
                                                                            gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u = SVector(1.1, -0.5, 2.34, 330.0, 1700)
    u_rr = SVector(0.9, 0.3, -1.2, 280.5, 1500)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        # Consistency with the physical flux
        for i in 1:2
            @test flux_lmars_es(u, u, normal_direction, equations)[i] ≈
                  flux(u, normal_direction, equations)
        end

        flux_left, flux_right = flux_lmars_es(u, u_rr, normal_direction, equations)
        @test flux_left == flux_right
        @test flux_lmars_es(u, u_rr, 2.5 * normal_direction, equations)[1] ≈ 2.5 * flux_left
    end

    # Check consistency for the 3D flux
    equations_3d = CompressibleEulerPotentialTemperatureEquationsWithGravity3D(c_p = 1004.0,
                                                                               c_v = 717.0,
                                                                               gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_3d = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1700)
    u_rr_3d = SVector(0.9, 0.3, -1.2, 1.7, 280.5, 1500)

    normal_directions_3d = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions_3d
        for i in 1:2
            @test flux_lmars_es(u_3d, u_3d, normal_direction, equations_3d)[i] ≈
                  flux(u_3d, normal_direction, equations_3d)
        end
        flux_left, flux_right = flux_lmars_es(u_3d, u_rr_3d, normal_direction, equations_3d)
        @test flux_left == flux_right
        @test flux_lmars_es(u_3d, u_rr_3d, 2.5 * normal_direction, equations_3d)[1] ≈
              2.5 * flux_left
    end

    # Check consistency between the 2D and 3D fluxes for different left and right states
    u_3d = SVector(u[1], u[2], u[3], 0.0, u[4], u[5])
    u_rr_3d = SVector(u_rr[1], u_rr[2], u_rr[3], 0.0, u_rr[4], u_rr[5])
    for normal_direction in normal_directions
        normal_3d = SVector(normal_direction[1], normal_direction[2], 0.0)
        flux_2d = flux_lmars_es(u, u_rr, normal_direction, equations)
        flux_3d = flux_lmars_es(u_3d, u_rr_3d, normal_3d, equations_3d)
        for i in 1:2
            @test flux_2d[i] ≈ flux_3d[i][[1, 2, 3, 5, 6]]
        end
    end
end

@testitem "Unit: Consistency check for gravity terms with energy with gravity: CEEEWG" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables states with the same temperature
    equations = CompressibleEulerEnergyEquationsWithGravity3D(c_p = 1004.0, c_v = 717.0,
                                                              gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = prim2cons(SVector(1.1, 2.0, -1.5, 0.5, 9.0e4, 1700.0), equations)
    u_rr = prim2cons(SVector(0.9, -1.0, 3.0, 1.2, 0.9 / 1.1 * 9.0e4, 2900.0), equations)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    # Both gravity terms coincide for an isothermal state
    for normal_direction in normal_directions
        @test flux_nonconservative_waruszewski_etal(u_ll, u_rr, normal_direction,
                                                    equations) ≈
              flux_nonconservative_isothermal(u_ll, u_rr, normal_direction, equations)
    end

    # Check consistency between the 2D and 3D gravity terms
    equations_2d = CompressibleEulerEnergyEquationsWithGravity2D(c_p = 1004.0, c_v = 717.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll_2d = prim2cons(SVector(1.1, 2.0, -1.5, 9.0e4, 1700.0), equations_2d)
    u_rr_2d = prim2cons(SVector(0.9, -1.0, 3.0, 7.5e4, 2900.0), equations_2d)
    u_ll_3d = prim2cons(SVector(1.1, 2.0, -1.5, 0.0, 9.0e4, 1700.0), equations)
    u_rr_3d = prim2cons(SVector(0.9, -1.0, 3.0, 0.0, 7.5e4, 2900.0), equations)
    for normal_direction in (SVector(1.0, 0.0), SVector(0.0, 1.0), SVector(0.5, -0.5),
                             SVector(-1.2, 0.3))
        normal_direction_3d = SVector(normal_direction[1], normal_direction[2], 0.0)
        for flux_gravity in (flux_nonconservative_waruszewski_etal,
                             flux_nonconservative_isothermal)
            @test flux_gravity(u_ll_2d, u_rr_2d, normal_direction, equations_2d) ≈
                  flux_gravity(u_ll_3d, u_rr_3d, normal_direction_3d,
                               equations)[[1, 2, 3, 5, 6]]
        end
    end
end

@testitem "Unit: Consistency check for 3D shallow water fluxes: SWE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = ShallowWaterEquations3D(gravity = 1.0)
    u = SVector(1.1, -0.5, 2.34, -3.5, 120.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        @test flux_wintermeyer_etal(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end

    for normal_direction in normal_directions
        @test flux_fjordholm_etal(u, u, normal_direction, equations) ≈
              flux(u, normal_direction, equations)
    end
end

@testitem "Unit: Consistency check for split covariant shallow water fluxes: SWE" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = SplitCovariantShallowWaterEquations2D(EARTH_GRAVITATIONAL_ACCELERATION,
                                                      EARTH_ROTATION_RATE)
    u = SVector(1.1, -0.5, 2.34)
    aux_vars = SVector{26}(ones(26))
    orientation = 1
    @test flux_ec(u, u, aux_vars, aux_vars, orientation, equations) ≈
          flux(u, aux_vars, orientation, equations)
end

@testitem "Unit: Consistency check for EC flux with Rainy Euler" setup=[Setup] tags=[:unit_fluxes] begin
    # Set up equations and dummy conservative variables state
    equations = CompressibleRainyEulerEquations2D()
    # Example state vector (ρ_d, ρ_m, ρ_r, ρu, ρv, ρe, ρq_v, ρq_c, T)
    u = SVector(1.0, 0.2, 0.1, 0.5, -0.4, 2.2, 0.1, 0.1, 300)

    normal_directions = [SVector(1.0, 0.0),
        SVector(0.0, 1.0),
        SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    for normal_direction in normal_directions
        equal = flux_ec_rain(u, u, normal_direction, equations) .≈
                flux(u, normal_direction, equations)
        # TODO
        expected = [true, true, true, true, true, false, true, true, true]
        @test equal == expected
    end
end

@testitem "Unit: check_axes for 2D manifolds in 3D" setup=[Setup] tags=[
    :unit_check_axes,
    :upstream
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "shallow_water", "cartesian",
                                 "elixir_unsteady_solid_body_rotation_ec_correction.jl"),
                        cells_per_dimension=(3, 3), maxiters=1)

    @testset "Cartesian form" begin
        mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
        u = Trixi.wrap_array(Trixi.compute_coefficients(0.0, semi), semi)

        @test Trixi.check_axes(u, mesh, equations, dg, cache) === nothing

        u_too_few = similar(u, size(u)[1:(end - 1)]..., size(u, ndims(u)) - 1)
        u_too_many = similar(u, size(u)[1:(end - 1)]..., size(u, ndims(u)) + 1)
        @test_throws DimensionMismatch Trixi.check_axes(u_too_few, mesh, equations, dg,
                                                        cache)
        @test_throws DimensionMismatch Trixi.check_axes(u_too_many, mesh, equations, dg,
                                                        cache)

        @test Trixi.ninterfaces(dg, cache) > 0
        for container in (cache.elements, cache.interfaces, cache.boundaries,
                          cache.mortars)
            @test Trixi.check_axes(container, equations, dg, cache) === nothing
        end

        # The element container of TrixiAtmo.jl is really checked, i.e., it does not use
        # the no-op fallback of Trixi.jl
        dg_wrong = DGSEM(polydeg = Trixi.polydeg(dg) + 1)
        @test_throws DimensionMismatch Trixi.check_axes(cache.elements, equations,
                                                        dg_wrong, cache)

        # Explicit bounds checks of TrixiAtmo.jl before assuming inbounds access
        @test_throws DimensionMismatch TrixiAtmo.calc_sources_2d_manifold_in_3d!(u_too_few,
                                                                                 u, 0.0,
                                                                                 semi.source_terms,
                                                                                 equations,
                                                                                 dg, cache)
    end

    @test_trixi_include(joinpath(EXAMPLES_DIR, "shallow_water", "covariant",
                                 "elixir_unsteady_solid_body_rotation_ec.jl"),
                        cells_per_dimension=(3, 3), maxiters=1)

    @testset "Covariant form" begin
        mesh, equations, dg, cache = Trixi.mesh_equations_solver_cache(semi)
        u = Trixi.wrap_array(Trixi.compute_coefficients(0.0, semi), semi)

        @test Trixi.check_axes(u, mesh, equations, dg, cache) === nothing

        u_too_few = similar(u, size(u)[1:(end - 1)]..., size(u, ndims(u)) - 1)
        u_too_many = similar(u, size(u)[1:(end - 1)]..., size(u, ndims(u)) + 1)
        @test_throws DimensionMismatch Trixi.check_axes(u_too_few, mesh, equations, dg,
                                                        cache)
        @test_throws DimensionMismatch Trixi.check_axes(u_too_many, mesh, equations, dg,
                                                        cache)

        @test Trixi.ninterfaces(dg, cache) > 0
        for container in (cache.elements, cache.interfaces, cache.boundaries,
                          cache.auxiliary_variables)
            @test Trixi.check_axes(container, equations, dg, cache) === nothing
        end

        # The containers of TrixiAtmo.jl are really checked, i.e., the element container
        # does not use the no-op fallback of Trixi.jl
        dg_wrong = DGSEM(polydeg = Trixi.polydeg(dg) + 1)
        @test_throws DimensionMismatch Trixi.check_axes(cache.elements, equations,
                                                        dg_wrong, cache)
        @test_throws DimensionMismatch Trixi.check_axes(cache.auxiliary_variables,
                                                        equations, dg_wrong, cache)

        # Explicit bounds checks of TrixiAtmo.jl before assuming inbounds access
        @test_throws DimensionMismatch Trixi.apply_jacobian!(nothing, u_too_few, mesh,
                                                             equations, dg, cache)
        @test_throws DimensionMismatch Trixi.calc_sources!(u_too_few, u, 0.0,
                                                           semi.source_terms, equations,
                                                           dg, cache)
    end
end

@testitem "Unit: Consistency check for ETEC flux with vector-invariant equations: CEVIEWG" setup=[Setup] tags=[:unit_fluxes] begin
    equations = CompressibleEulerVectorInvariantEquationsWithGravity3D(c_p = 1004.0,
                                                                       c_v = 717.0,
                                                                       gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5, 1700.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        # Consistency with the conservative part of the physical flux
        rho, v1, v2, v3, rho_theta = u_ll
        v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
                   v3 * normal_direction[3]
        kinetic_energy = 0.5 * (v1^2 + v2^2 + v3^2)
        flux_physical = SVector(rho * v_normal, kinetic_energy * normal_direction[1],
                                kinetic_energy * normal_direction[2],
                                kinetic_energy * normal_direction[3],
                                rho_theta * v_normal, 0.0)
        for i in 1:2
            @test flux_etec(u_ll, u_ll, normal_direction, equations)[i] ≈ flux_physical
        end

        # The nonconservative contribution is antisymmetric
        flux_left, flux_right = flux_etec(u_ll, u_rr, normal_direction, equations)
        flux_left_swapped, flux_right_swapped = flux_etec(u_rr, u_ll, normal_direction,
                                                          equations)
        @test flux_left ≈ flux_right_swapped
        @test flux_right ≈ flux_left_swapped

        # The `FluxTurbo` specialization returns the same contributions
        turbo_ll = Trixi.cons2turbo(flux_etec, u_ll..., equations)
        turbo_rr = Trixi.cons2turbo(flux_etec, u_rr..., equations)
        flux_left_turbo, flux_right_turbo = Trixi.flux_turbo(flux_etec, turbo_ll...,
                                                             turbo_rr...,
                                                             normal_direction...,
                                                             equations)
        @test SVector(flux_left_turbo) ≈ flux_left
        @test SVector(flux_right_turbo) ≈ flux_right
    end
    @test Trixi.nturbovars(flux_etec, equations) == Val(8)

    # The entropy variables are the derivatives of the entropy
    epsilon = 1.0e-6
    for v in (1, 5)
        e_v = SVector(ntuple(i -> i == v ? epsilon : 0.0, 6))
        @test isapprox(cons2entropy(u_ll, equations)[v],
                       (entropy(u_ll + e_v, equations) -
                        entropy(u_ll - e_v, equations)) / (2 * epsilon); rtol = 1.0e-6)
    end
end

@testitem "Unit: Consistency check for ES LMARS flux with vector-invariant equations: CEVIEWG" setup=[Setup] tags=[:unit_fluxes] begin
    equations = CompressibleEulerVectorInvariantEquationsWithGravity3D(c_p = 1004.0,
                                                                       c_v = 717.0,
                                                                       gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    u_ll = SVector(1.1, -0.5, 2.34, 2.4, 330.0, 1500.0)
    u_rr = SVector(0.9, 0.3, -1.2, 1.7, 280.5, 1500.0)

    normal_directions = [SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 1.0),
        SVector(0.5, -0.5, 0.2),
        SVector(-1.2, 0.3, 1.4)]

    for normal_direction in normal_directions
        # Consistency with the conservative part of the physical flux
        rho, v1, v2, v3, rho_theta = u_ll
        v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] +
                   v3 * normal_direction[3]
        kinetic_energy = 0.5 * (v1^2 + v2^2 + v3^2)
        flux_physical = SVector(rho * v_normal, kinetic_energy * normal_direction[1],
                                kinetic_energy * normal_direction[2],
                                kinetic_energy * normal_direction[3],
                                rho_theta * v_normal, 0.0)
        for i in 1:2
            @test flux_lmars_es(u_ll, u_ll, normal_direction, equations)[i] ≈ flux_physical
        end

        # The flux is homogeneous of degree one in the normal direction
        for i in 1:2
            @test flux_lmars_es(u_ll, u_rr, 2.5 * normal_direction, equations)[i] ≈
                  2.5 * flux_lmars_es(u_ll, u_rr, normal_direction, equations)[i]
        end
    end
end

@testitem "Unit: Consistency check for flux_etec of rainy vector-invariant equations: CREVIE" setup=[Setup] tags=[:unit_fluxes] begin
    equations = CompressibleRainyEulerVectorInvariantEquations2D(c_pd = 1004.0,
                                                                 c_vd = 717.0,
                                                                 R_v = 461.0,
                                                                 gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    # (rho, v1, v2, rho_theta_v, rho_vapor, rho_cloud, rho_rain, phi)
    u_ll = SVector(1.1, -0.5, 2.34, 330.0, 0.012, 0.002, 0.003, 1500.0)
    u_rr = SVector(0.9, 0.3, -1.2, 280.5, 0.008, 0.001, 0.004, 1700.0)
    normal_directions = [SVector(1.0, 0.0), SVector(0.0, 1.0), SVector(0.5, -0.5),
        SVector(-1.2, 0.3)]

    # Entropy variables of the total energy (up to constant water components, whose jumps
    # vanish) and the corresponding entropy potential ρ K v_n
    function energy_variables(u, equations)
        rho, v1, v2, _, _, _, _, phi = u
        kinetic_energy = 0.5 * (v1^2 + v2^2)
        return SVector(kinetic_energy + phi, rho * v1, rho * v2,
                       equations.c_pd * TrixiAtmo.exner_pressure(u, equations),
                       0.0, 0.0, 0.0, 0.0)
    end
    energy_potential(u, n) = 0.5 * u[1] * (u[2]^2 + u[3]^2) * (u[2] * n[1] + u[3] * n[2])

    # Residual of the discrete condition [[w]] ⋅ f - {{w}} ⋅ g = [[ψ]] for the conservative
    # part f = (flux_left + flux_right) / 2 and the nonconservative part
    # g = flux_left - flux_right
    function residual(w_ll, w_rr, u_ll, u_rr, n, psi_jump, equations)
        flux_left, flux_right = flux_etec(u_ll, u_rr, n, equations)
        f = 0.5 * (flux_left + flux_right)
        g = flux_left - flux_right
        return sum((w_rr - w_ll) .* f) - sum(0.5 * (w_ll + w_rr) .* g) - psi_jump
    end

    for normal_direction in normal_directions
        # Consistency and symmetry
        flux_left, flux_right = flux_etec(u_ll, u_ll, normal_direction, equations)
        @test flux_left ≈ flux(u_ll, normal_direction, equations)
        @test flux_right ≈ flux(u_ll, normal_direction, equations)
        flux_left, flux_right = flux_etec(u_ll, u_rr, normal_direction, equations)
        flux_left_swapped, flux_right_swapped = flux_etec(u_rr, u_ll, normal_direction,
                                                          equations)
        @test flux_left ≈ flux_right_swapped
        @test flux_right ≈ flux_left_swapped
        # Homogeneity of degree one in the normal direction
        flux_left_scaled, _ = flux_etec(u_ll, u_rr, 2.5 * normal_direction, equations)
        @test flux_left_scaled ≈ 2.5 * flux_left

        # Total energy conservation without rain, for jumping mass fractions
        u_a = Base.setindex(u_ll, 0.0, 7)
        u_b = Base.setindex(u_rr, 0.0, 7)
        psi_jump = energy_potential(u_b, normal_direction) -
                   energy_potential(u_a, normal_direction)
        @test isapprox(residual(energy_variables(u_a, equations),
                                energy_variables(u_b, equations),
                                u_a, u_b, normal_direction, psi_jump, equations),
                       0.0; atol = 1.0e-8)

        # Entropy conservation without rain, for continuous mass fractions
        rho_b = 0.9
        u_b = SVector(rho_b, 0.3, -1.2, 280.5, rho_b * u_a[5] / u_a[1],
                      rho_b * u_a[6] / u_a[1], 0.0, 1700.0)
        @test isapprox(residual(cons2entropy(u_a, equations), cons2entropy(u_b, equations),
                                u_a, u_b, normal_direction, 0.0, equations),
                       0.0; atol = 1.0e-8)

        # With rain, the total energy defect is m [[ρ]] |[[v]]|^2 / (4 {{ρ}}) n_2: the
        # frictional heating balances the potential energy released by the rain exactly
        rain_mass_flux_avg = 0.5 *
                             (u_ll[7] * TrixiAtmo.terminal_velocity_rain(u_ll, equations) +
                              u_rr[7] * TrixiAtmo.terminal_velocity_rain(u_rr, equations))
        defect = (rain_mass_flux_avg * normal_direction[2] * (u_rr[1] - u_ll[1]) *
                  ((u_rr[2] - u_ll[2])^2 + (u_rr[3] - u_ll[3])^2) /
                  (2 * (u_ll[1] + u_rr[1])))
        psi_jump = energy_potential(u_rr, normal_direction) -
                   energy_potential(u_ll, normal_direction)
        @test isapprox(residual(energy_variables(u_ll, equations),
                                energy_variables(u_rr, equations),
                                u_ll, u_rr, normal_direction, psi_jump, equations),
                       defect; atol = 1.0e-8)
    end

    # The pressure gradient and gravity terms vanish for both hydrostatic states
    x_ll = SVector(100.0, 250.0)
    x_rr = SVector(130.0, 410.0)
    normal_direction = SVector(0.3, 1.1)
    for initial_condition in (initial_condition_isothermal, initial_condition_adiabatic)
        u_bottom = initial_condition(x_ll, 0.0, equations)
        u_top = initial_condition(x_rr, 0.0, equations)
        flux_left, flux_right = flux_etec(u_bottom, u_top, normal_direction, equations)
        @test isapprox(flux_left[2] - flux_right[2], 0.0; atol = 1.0e-10)
        @test isapprox(flux_left[3] - flux_right[3], 0.0; atol = 1.0e-10)
    end
end
