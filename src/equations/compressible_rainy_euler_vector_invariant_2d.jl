@muladd begin
#! format: noindent

@doc raw"""
    CompressibleRainyEulerVectorInvariantEquations2D(; c_pd = 1004, c_vd = 717, R_v = 461,
                                                     gravity = EARTH_GRAVITATIONAL_ACCELERATION)

The compressible Euler equations with warm-rain microphysics in two space dimensions in
vector-invariant form, with the virtual potential temperature as thermodynamic prognostic
variable,
```math
\begin{aligned}
\partial_t \varrho + \nabla \cdot (\varrho \vec{v} - \varrho_r v_r \vec{e}_z) &= 0, \\
\partial_t \vec{v} + (\nabla \times \vec{v}) \times \vec{v}
+ \nabla \left( \tfrac{1}{2} \lVert \vec{v} \rVert^2 \right)
+ c_{pd} \theta_v \nabla \Pi + \nabla \phi
- \frac{\varrho_r v_r}{\varrho} \partial_z \vec{v} &= 0, \\
\partial_t (\varrho \theta_v) + \nabla \cdot (\varrho \theta_v \vec{v})
- \frac{\varrho_r v_r}{c_{pd} \Pi} \partial_z \phi
&= \frac{L_v^0 (\mathcal{C} - \mathcal{E})}{c_{pd} \Pi}, \\
\partial_t \varrho_v + \nabla \cdot (\varrho_v \vec{v}) &= -\mathcal{C} + \mathcal{E}, \\
\partial_t \varrho_c + \nabla \cdot (\varrho_c \vec{v}) &= \mathcal{C} - \mathcal{A}, \\
\partial_t \varrho_r + \nabla \cdot (\varrho_r (\vec{v} - v_r \vec{e}_z))
&= \mathcal{A} - \mathcal{E},
\end{aligned}
```
where the vertical direction ``z`` is the second coordinate direction. Here, ``\varrho`` is
the total density (dry air, water vapor ``\varrho_v``, cloud water ``\varrho_c``, and rain
``\varrho_r``), ``\vec{v} = (v_1, v_2)`` the velocity of the gas and the cloud water,
``v_r \ge 0`` the terminal fall speed of rain, ``\theta_v`` the virtual potential
temperature including the condensate loading, ``\phi = g z`` the geopotential (stored as
auxiliary variable), and
```math
\Pi = \left( \frac{p}{p_0} \right)^{R_d / c_{pd}}
    = \left( \frac{R_d \varrho \theta_v}{p_0} \right)^{R_d / c_{vd}},
\qquad
p = p_0 \left( \frac{R_d \varrho \theta_v}{p_0} \right)^{\gamma}
```
the Exner pressure and the pressure with ``\gamma = c_{pd} / c_{vd}``.

The thermodynamics is the constant-``\kappa`` approximation of Eldred, Taylor, and Guba
(2022): the heat capacities of water vapor are ``c_{pd} R_v / R_d`` and
``c_{vd} R_v / R_d``, and that of liquid water is zero. Thus, ``\kappa = R_d / c_{pd}`` and
``\gamma`` do not depend on the composition, the equation of state has the dry form above,
and the latent heat of vaporization is ``L_v(T) = c_{pd} R_v / R_d \, T + L_v^0``. The
saturation vapor pressure ``e_s(T)`` is computed consistently from the Clausius-Clapeyron
relation with this latent heat. The temperature is ``T = p / (R_d \varrho_d + R_v \varrho_v)``.

The conversion rates ``\mathcal{C}`` (condensation, see [`SourceTermsWarmRain`](@ref)),
``\mathcal{A}`` (autoconversion and accretion), and ``\mathcal{E}`` (rain evaporation), and
the terminal velocity of rain follow the Kessler-type warm-rain scheme of
[`CompressibleRainyEulerEquations2D`](@ref). The frictional heating by falling rain is
written as the nonconservative product
``-\varrho_r v_r / (c_{pd} \Pi) \, \partial_z \phi`` so that the hyperbolic part conserves
the total energy.

The variables are `(rho, v1, v2, rho_theta_v, rho_vapor, rho_cloud, rho_rain, phi)`.
Only `normal_direction` methods are implemented, i.e., the equations are intended for the
`P4estMesh`.

## References
- Christopher Eldred, Mark Taylor, Oksana Guba (2022)
  Thermodynamically consistent versions of approximations used in modelling moist air
  [DOI: 10.1002/qj.4318](https://doi.org/10.1002/qj.4318)
- Sabine Doppler, Philip L. Lederer, Joachim Schöberl, Henry von Wahl (2024)
  A discontinuous Galerkin approach for atmospheric flows with implicit condensation
  [DOI: 10.1016/j.jcp.2023.112713](https://doi.org/10.1016/j.jcp.2023.112713)
- Marco Artiano, Kieran Ricardo, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2026)
  Entropy-Stable and Well-Balanced Discontinuous Galerkin Methods for the Compressible
  Euler Equations in Vector-Invariant Form
  [arXiv: 2610.05820](https://arxiv.org/abs/2610.05820)
"""
struct CompressibleRainyEulerVectorInvariantEquations2D{RealT <: Real} <:
       AbstractCompressibleRainyEulerEquations{2, 8}
    p_0::RealT # reference pressure in Pa
    c_pd::RealT # specific heat of dry air at constant pressure in J/(kg K)
    c_vd::RealT # specific heat of dry air at constant volume in J/(kg K)
    R_d::RealT # gas constant of dry air in J/(kg K)
    R_v::RealT # gas constant of water vapor in J/(kg K)
    c_pv::RealT # = c_pd * R_v / R_d; constant-kappa heat capacity of vapor in J/(kg K)
    gamma::RealT # = c_pd / c_vd; ratio of specific heats of the mixture
    inv_gamma_minus_one::RealT # = inv(gamma - 1)
    K::RealT # = p_0 * (R_d / p_0)^gamma; scaling factor between p and (rho theta_v)^gamma
    g::RealT # gravitational acceleration in m/s²
    T_ref::RealT # reference temperature in K
    e_s_ref::RealT # saturation vapor pressure at T_ref in Pa
    L_ref::RealT # latent heat of vaporization at T_ref in J/kg
    L_0::RealT # = L_ref - c_pv * T_ref; constant part of the latent heat in J/kg
    rho_liquid_water::RealT # density of liquid water in kg/m³
    rain_water_distr::RealT # intercept parameter N_0 of the Marshall-Palmer distribution
    v_mean_rain::RealT # coefficient v_0 of the fall speed of raindrops v(D) = v_0 D^(1/2)
end

function CompressibleRainyEulerVectorInvariantEquations2D(; c_pd = 1004, c_vd = 717,
                                                          R_v = 461,
                                                          gravity = EARTH_GRAVITATIONAL_ACCELERATION)
    c_pd, c_vd, R_v, g = promote(float(c_pd), float(c_vd), float(R_v), float(gravity))
    RealT = typeof(c_pd)
    p_0 = convert(RealT, 100_000)
    R_d = c_pd - c_vd
    c_pv = c_pd * R_v / R_d
    gamma = c_pd / c_vd
    inv_gamma_minus_one = inv(gamma - 1)
    K = p_0 * (R_d / p_0)^gamma
    T_ref = convert(RealT, 273.15)
    e_s_ref = convert(RealT, 610.7)
    L_ref = convert(RealT, 2.5e6)
    L_0 = L_ref - c_pv * T_ref
    rho_liquid_water = convert(RealT, 1000)
    rain_water_distr = convert(RealT, 8.0e6)
    v_mean_rain = convert(RealT, 130)
    return CompressibleRainyEulerVectorInvariantEquations2D{RealT}(p_0, c_pd, c_vd, R_d,
                                                                   R_v, c_pv, gamma,
                                                                   inv_gamma_minus_one,
                                                                   K, g, T_ref, e_s_ref,
                                                                   L_ref, L_0,
                                                                   rho_liquid_water,
                                                                   rain_water_distr,
                                                                   v_mean_rain)
end

# Together with the specialization of `Adapt.adapt_structure` in Trixi.jl,
# this allows to adapt the floating point type of the equations, e.g., to `Float32`.
function Base.similar(equations::CompressibleRainyEulerVectorInvariantEquations2D,
                      ::Type{NewRealT}) where {NewRealT}
    c_pd, c_vd, R_v, gravity = convert.(NewRealT,
                                        (equations.c_pd, equations.c_vd, equations.R_v,
                                         equations.g))
    return CompressibleRainyEulerVectorInvariantEquations2D(; c_pd, c_vd, R_v, gravity)
end

function varnames(::typeof(cons2cons),
                  ::CompressibleRainyEulerVectorInvariantEquations2D)
    ("rho", "v1", "v2", "rho_theta_v", "rho_vapor", "rho_cloud", "rho_rain", "phi")
end

function varnames(::typeof(cons2prim),
                  ::CompressibleRainyEulerVectorInvariantEquations2D)
    ("rho", "v1", "v2", "p", "rho_vapor", "rho_cloud", "rho_rain", "phi")
end

function have_nonconservative_terms(::CompressibleRainyEulerVectorInvariantEquations2D)
    return Trixi.True()
end

###############################################################################
# Thermodynamics and microphysics

@inline function pressure(u, equations::CompressibleRainyEulerVectorInvariantEquations2D)
    return equations.K * u[4]^equations.gamma
end

# Exner pressure Π = (p / p_0)^(R_d / c_pd) = (R_d ρθ_v / p_0)^(R_d / c_vd)
@inline function exner_pressure(u,
                                equations::CompressibleRainyEulerVectorInvariantEquations2D)
    return (equations.R_d * u[4] / equations.p_0)^(equations.R_d / equations.c_vd)
end

# Temperature T = p / (R_d ρ_d + R_v ρ_v)
@inline function temperature(u, equations::CompressibleRainyEulerVectorInvariantEquations2D)
    rho, _, _, _, rho_vapor, rho_cloud, rho_rain = u
    rho_dry = rho - rho_vapor - rho_cloud - rho_rain
    return pressure(u, equations) /
           (equations.R_d * rho_dry + equations.R_v * rho_vapor)
end

# Virtual temperature T_v = θ_v Π = p / (R_d ρ)
@inline function virtual_temperature(u,
                                     equations::CompressibleRainyEulerVectorInvariantEquations2D)
    return pressure(u, equations) / (equations.R_d * u[1])
end

# Saturation vapor pressure from the Clausius-Clapeyron relation with the latent heat
# L_v(T) = c_pv T + L_0 of the constant-kappa system, i.e.,
# e_s(T) = e_s(T_ref) (T / T_ref)^(c_pd / R_d) exp(L_0 / R_v (1 / T_ref - 1 / T)).
# Here, c_pv / R_v = c_pd / R_d.
@inline function saturation_vapor_pressure(temperature,
                                           equations::CompressibleRainyEulerVectorInvariantEquations2D)
    @unpack c_pd, R_d, R_v, T_ref, e_s_ref, L_0 = equations
    return e_s_ref * (temperature / T_ref)^(c_pd / R_d) *
           exp(L_0 / R_v * (1 / T_ref - 1 / temperature))
end

# Mass-weighted terminal fall speed of rain for a Marshall-Palmer distribution with
# intercept parameter N_0 and fall speeds of raindrops v(D) = v_0 D^(1/2), i.e.,
# v_r = v_0 Γ(4.5) / 6 (ρ_r / (π ρ_w N_0))^(1/8), where ρ_w is the density of liquid
# water.
# Note that `terminal_velocity_rain` of `CompressibleRainyEulerEquations2D` uses the
# density of moist air and rain instead of ρ_w.
@inline function terminal_velocity_rain(rho_rain,
                                        equations::CompressibleRainyEulerVectorInvariantEquations2D)
    RealT = typeof(rho_rain)
    if rho_rain > 0
        # Γ(4.5) / 6 ≈ 1.9386213994279082
        return equations.v_mean_rain * convert(RealT, 1.9386213994279082) *
               (rho_rain /
                (convert(RealT, pi) * equations.rho_liquid_water *
                 equations.rain_water_distr))^(0.125f0)
    else
        return zero(RealT)
    end
end

###############################################################################
# Conversions

@inline function Trixi.prim2cons(prim,
                                 equations::CompressibleRainyEulerVectorInvariantEquations2D)
    rho, v1, v2, p, rho_vapor, rho_cloud, rho_rain, phi = prim
    rho_theta_v = (p / equations.K)^(1 / equations.gamma)
    return SVector(rho, v1, v2, rho_theta_v, rho_vapor, rho_cloud, rho_rain, phi)
end

@inline function cons2prim(u, equations::CompressibleRainyEulerVectorInvariantEquations2D)
    rho, v1, v2, _, rho_vapor, rho_cloud, rho_rain, phi = u
    return SVector(rho, v1, v2, pressure(u, equations), rho_vapor, rho_cloud, rho_rain,
                   phi)
end

# Specific entropies of dry air and water vapor of the constant-kappa system,
#   s_d = c_pd log(T / T_ref) - R_d log(p_d / p_0),
#   s_v = c_pv log(T / T_ref) - R_v log(p_v / e_s(T_ref)) + L_ref / T_ref,
# and s_l = 0 for liquid water (cloud and rain), where p_d = ρ_d R_d T and
# p_v = ρ_v R_v T.
# The constants are chosen such that the chemical potentials of vapor and liquid water are
# equal at saturation.
@inline function specific_entropies(u,
                                    equations::CompressibleRainyEulerVectorInvariantEquations2D)
    @unpack c_pd, c_pv, R_d, R_v, p_0, T_ref, e_s_ref, L_ref = equations
    rho, _, _, _, rho_vapor, rho_cloud, rho_rain = u
    rho_dry = rho - rho_vapor - rho_cloud - rho_rain
    T = temperature(u, equations)
    s_dry = c_pd * log(T / T_ref) - R_d * log(rho_dry * R_d * T / p_0)
    s_vapor = c_pv * log(T / T_ref) - R_v * log(rho_vapor * R_v * T / e_s_ref) +
              L_ref / T_ref
    return s_dry, s_vapor, T
end

"""
    entropy(u, equations::CompressibleRainyEulerVectorInvariantEquations2D)

Mathematical entropy ``U = -\\varrho s`` given by the negative thermodynamic entropy density
of the constant-``\\kappa`` system. It is convex, but not strictly convex, and does not
depend on the velocity. For vanishing vapor density, the limit
``\\varrho_v s_v \\to 0`` is used.
"""
@inline function entropy(u, equations::CompressibleRainyEulerVectorInvariantEquations2D)
    rho, _, _, _, rho_vapor, rho_cloud, rho_rain = u
    rho_dry = rho - rho_vapor - rho_cloud - rho_rain
    s_dry, s_vapor, _ = specific_entropies(u, equations)
    rho_s_vapor = rho_vapor > 0 ? rho_vapor * s_vapor : zero(rho_vapor)
    return -(rho_dry * s_dry + rho_s_vapor)
end

# Entropy variables w = ∂U/∂u of the entropy U = -ρ s, see `entropy`. With
# σ_n = s_n - c_pn (σ_l = 0 for liquid water) and ρ_d = ρ - ρ_v - ρ_c - ρ_r, they are
#   w = (-σ_d, 0, 0, -c_pd Π / T, σ_d - σ_v, σ_d, σ_d, 0).
# The velocity components vanish since U does not depend on the velocity. The vapor
# component is singular for vanishing vapor density.
@inline function cons2entropy(u, equations::CompressibleRainyEulerVectorInvariantEquations2D)
    s_dry, s_vapor, T = specific_entropies(u, equations)
    sigma_dry = s_dry - equations.c_pd
    sigma_vapor = s_vapor - equations.c_pv
    w4 = -equations.c_pd * exner_pressure(u, equations) / T

    zero_ = zero(eltype(u))
    return SVector(-sigma_dry, zero_, zero_, w4, sigma_dry - sigma_vapor, sigma_dry,
                   sigma_dry, zero_)
end

###############################################################################
# Fluxes and wave speeds

# Conservative part of the flux in the normal direction. The vertical direction is the
# second coordinate direction. The flux of the momentum equation is the gradient of the
# kinetic energy (Bernoulli term).
@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleRainyEulerVectorInvariantEquations2D)
    rho, v1, v2, rho_theta_v, rho_vapor, rho_cloud, rho_rain = u
    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    v_rain = terminal_velocity_rain(rho_rain, equations)
    rain_flux = rho_rain * v_rain * normal_direction[2]
    kinetic_energy = 0.5f0 * (v1^2 + v2^2)

    f1 = rho * v_normal - rain_flux
    f2 = kinetic_energy * normal_direction[1]
    f3 = kinetic_energy * normal_direction[2]
    f4 = rho_theta_v * v_normal
    f5 = rho_vapor * v_normal
    f6 = rho_cloud * v_normal
    f7 = rho_rain * v_normal - rain_flux

    return SVector(f1, f2, f3, f4, f5, f6, f7, zero(eltype(u)))
end

# Nonconservative terms in the standard (strong-form) way: the factors H of the
# nonconservative products H ∇g are evaluated at the own node and multiplied by the jumps
# of g, i.e., the vorticity terms, the transport of momentum by falling rain, and the
# frictional heating by falling rain. The pressure gradient and the gravity term are
# discretized as `c_pd * theta_v_num * jump(Π) + jump(φ)` with the two-point mean
# `theta_v_num` of the virtual potential temperature passed as argument, which determines
# the hydrostatic states that are preserved discretely.
@inline function _flux_nonconservative_vector_invariant(u_ll, u_rr,
                                                        normal_direction::AbstractVector,
                                                        theta_v_num,
                                                        equations::CompressibleRainyEulerVectorInvariantEquations2D)
    rho_ll, v1_ll, v2_ll, _, _, _, rho_rain_ll, phi_ll = u_ll
    _, v1_rr, v2_rr, _, _, _, _, phi_rr = u_rr
    exner_ll = exner_pressure(u_ll, equations)
    exner_rr = exner_pressure(u_rr, equations)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_phi = phi_rr - phi_ll

    # Vorticity terms (∇ × v) × v = ω (-v2, v1) with ω = ∂_x v2 - ∂_y v1
    vorticity = normal_direction[1] * jump_v2 - normal_direction[2] * jump_v1
    g2 = -v2_ll * vorticity
    g3 = v1_ll * vorticity

    # Transport of momentum by falling rain, -(ρ_r v_r / ρ) ∂_z v
    rain_mass_flux_ll = rho_rain_ll * terminal_velocity_rain(rho_rain_ll, equations)
    rain_advection = rain_mass_flux_ll / rho_ll * normal_direction[2]
    g2 = g2 - rain_advection * jump_v1
    g3 = g3 - rain_advection * jump_v2

    # Pressure gradient and gravity, c_pd θ_v ∇Π + ∇φ
    gradient = equations.c_pd * theta_v_num * (exner_rr - exner_ll) + jump_phi
    g2 = g2 + gradient * normal_direction[1]
    g3 = g3 + gradient * normal_direction[2]

    # Frictional heating by falling rain, -(ρ_r v_r / (c_pd Π)) ∂_z φ
    g4 = -rain_mass_flux_ll / (equations.c_pd * exner_ll) * normal_direction[2] *
         jump_phi

    zero_ = zero(eltype(u_ll))
    return SVector(zero_, g2, g3, g4, zero_, zero_, zero_, zero_)
end

"""
    flux_nonconservative_isothermal(u_ll, u_rr, normal_direction::AbstractVector,
                                    equations::CompressibleRainyEulerVectorInvariantEquations2D)

Nonconservative terms of the [`CompressibleRainyEulerVectorInvariantEquations2D`](@ref).
The pressure gradient and gravity terms are discretized as
``c_{pd} \\theta_v^{num} [\\![\\Pi]\\!] + [\\![\\phi]\\!]`` with the inverse logarithmic
mean ``\\theta_v^{num} = 1 / \\{\\!\\{ 1 / \\theta_v \\}\\!\\}_{\\log}``. This preserves
hydrostatic states at rest with constant virtual temperature ``T_v = \\theta_v \\Pi`` and
with constant virtual potential temperature ``\\theta_v``, analogously to the isothermal and
isentropic
balances of the dry equations in
- Marco Artiano, Kieran Ricardo, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2026)
  Entropy-Stable and Well-Balanced Discontinuous Galerkin Methods for the Compressible
  Euler Equations in Vector-Invariant Form
  [arXiv: 2610.05820](https://arxiv.org/abs/2610.05820)

The remaining nonconservative terms (vorticity, transport of momentum by falling rain, and
frictional heating by falling rain) are discretized in the standard way with the local
factors of the own state multiplied by the jumps.
"""
@inline function flux_nonconservative_isothermal(u_ll, u_rr,
                                                 normal_direction::AbstractVector,
                                                 equations::CompressibleRainyEulerVectorInvariantEquations2D)
    theta_v_num = inv_ln_mean(u_ll[1] / u_ll[4], u_rr[1] / u_rr[4])
    return _flux_nonconservative_vector_invariant(u_ll, u_rr, normal_direction,
                                                  theta_v_num, equations)
end

"""
    flux_nonconservative_isentropic(u_ll, u_rr, normal_direction::AbstractVector,
                                    equations::CompressibleRainyEulerVectorInvariantEquations2D)

Nonconservative terms of the [`CompressibleRainyEulerVectorInvariantEquations2D`](@ref).
The pressure gradient and gravity terms are discretized as
``c_{pd} \\theta_v^{num} [\\![\\Pi]\\!] + [\\![\\phi]\\!]`` with the arithmetic mean
``\\theta_v^{num} = \\{\\!\\{ \\theta_v \\}\\!\\}``. This preserves hydrostatic states at
rest with constant virtual potential temperature ``\\theta_v``, analogously to the
isentropic balance of the dry equations in
- Marco Artiano, Kieran Ricardo, Oswald Knoth, Peter Spichtinger, Hendrik Ranocha (2026)
  Entropy-Stable and Well-Balanced Discontinuous Galerkin Methods for the Compressible
  Euler Equations in Vector-Invariant Form
  [arXiv: 2610.05820](https://arxiv.org/abs/2610.05820)

The remaining nonconservative terms are discretized as in
[`flux_nonconservative_isothermal`](@ref).
"""
@inline function flux_nonconservative_isentropic(u_ll, u_rr,
                                                 normal_direction::AbstractVector,
                                                 equations::CompressibleRainyEulerVectorInvariantEquations2D)
    theta_v_num = 0.5f0 * (u_ll[4] / u_ll[1] + u_rr[4] / u_rr[1])
    return _flux_nonconservative_vector_invariant(u_ll, u_rr, normal_direction,
                                                  theta_v_num, equations)
end

# The speed of sound c = sqrt(γ p / ρ) of the constant-kappa system does not depend on the
# composition. Rain moves with the additional vertical velocity -v_r.
@inline function max_abs_speeds(u,
                                equations::CompressibleRainyEulerVectorInvariantEquations2D)
    rho, v1, v2 = u
    c = sqrt(equations.gamma * pressure(u, equations) / rho)
    v_rain = terminal_velocity_rain(u[7], equations)
    return abs(v1) + c, abs(v2) + c + v_rain
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleRainyEulerVectorInvariantEquations2D)
    rho_ll, v1_ll, v2_ll = u_ll
    rho_rr, v1_rr, v2_rr = u_rr

    v_ll = abs(v1_ll * normal_direction[1] + v2_ll * normal_direction[2]) +
           terminal_velocity_rain(u_ll[7], equations) * abs(normal_direction[2])
    v_rr = abs(v1_rr * normal_direction[1] + v2_rr * normal_direction[2]) +
           terminal_velocity_rain(u_rr[7], equations) * abs(normal_direction[2])

    c_ll = sqrt(equations.gamma * pressure(u_ll, equations) / rho_ll)
    c_rr = sqrt(equations.gamma * pressure(u_rr, equations) / rho_rr)

    return max(v_ll, v_rr) + max(c_ll, c_rr) * norm(normal_direction)
end

# Slip-wall boundary condition: the velocity is mirrored at the wall. The surface fluxes
# are a tuple of the conservative and the nonconservative surface flux. Falling rain leaves
# the domain through the bottom wall since the rain flux is not mirrored.
@inline function boundary_condition_slip_wall(u_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_functions,
                                              equations::CompressibleRainyEulerVectorInvariantEquations2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # compute the normal velocity
    v_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * v_normal * normal[1],
                         u_inner[3] - 2 * v_normal * normal[2],
                         u_inner[4], u_inner[5], u_inner[6], u_inner[7], u_inner[8])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)
    return flux, noncons_flux
end

###############################################################################
# Source terms

"""
    SourceTermsWarmRain(; condensation_time_scale)

Source terms of the warm-rain microphysics for the
[`CompressibleRainyEulerVectorInvariantEquations2D`](@ref). They contain
- condensation of vapor and evaporation of cloud water, modeled as relaxation towards
  saturation, ``\\mathcal{C} = \\max(\\varrho_v - \\varrho_{vs}(T), -\\varrho_c) / \\tau``,
  with the time scale ``\\tau`` = `condensation_time_scale`,
- autoconversion ``10^{-3}\\,\\mathrm{s}^{-1}\\,\\varrho_c`` and accretion
  ``1.72\\, \\varrho_c \\varrho_r^{7/8}`` of cloud water to rain,
- evaporation of rain
  ``\\mathcal{E} = (3.86 \\cdot 10^{-3} - 9.41 \\cdot 10^{-5} (T - T_{ref}))
  (1 + 9.1\\, \\varrho_r^{3/16}) (\\varrho_{vs}(T) - \\varrho_v) \\varrho_r^{1/2}``,
- the corresponding latent heating ``L_v^0 (\\mathcal{C} - \\mathcal{E}) / (c_{pd} \\Pi)``
  in the equation of the virtual potential temperature.

The rates of autoconversion, accretion, and rain evaporation are those of
[`CompressibleRainyEulerEquations2D`](@ref). Since the condensation rate has the sign of
``\\varrho_v - \\varrho_{vs}``, the phase changes do not decrease the thermodynamic entropy.
Gravity is not part of the source terms, but of the nonconservative terms.
"""
struct SourceTermsWarmRain{RealT <: Real}
    condensation_time_scale::RealT
end

function SourceTermsWarmRain(; condensation_time_scale)
    return SourceTermsWarmRain(condensation_time_scale)
end

@inline function (source::SourceTermsWarmRain)(u, x, t,
                                               equations::CompressibleRainyEulerVectorInvariantEquations2D)
    RealT = eltype(u)
    _, _, _, _, rho_vapor, rho_cloud, rho_rain = u
    T = temperature(u, equations)
    exner = exner_pressure(u, equations)
    rho_vs = saturation_vapor_pressure(T, equations) / (equations.R_v * T)
    rho_rain_pos = max(rho_rain, zero(RealT))

    # Condensation (positive) or evaporation of cloud water (negative)
    condensation = max(rho_vapor - rho_vs, -rho_cloud) /
                   convert(RealT, source.condensation_time_scale)

    # Autoconversion and accretion of cloud water
    conversion = convert(RealT, 1.0e-3) * rho_cloud +
                 convert(RealT, 1.72) * rho_cloud * rho_rain_pos^(0.875f0)

    # Evaporation of rain
    evaporation = (convert(RealT, 3.86e-3) -
                   convert(RealT, 9.41e-5) * (T - equations.T_ref)) *
                  (1 + convert(RealT, 9.1) * rho_rain_pos^(0.1875f0)) *
                  (rho_vs - rho_vapor) * sqrt(rho_rain_pos)

    heating = equations.L_0 * (condensation - evaporation) / (equations.c_pd * exner)

    zero_ = zero(RealT)
    return SVector(zero_, zero_, zero_, heating, -condensation + evaporation,
                   condensation - conversion, conversion - evaporation, zero_)
end

###############################################################################
# Hydrostatic initial conditions at rest

# Vapor density of the hydrostatic initial conditions below: a specific humidity decaying
# exponentially with height, which is subsaturated in the default setups. The water species
# enter the hydrostatic balance only through the total density.
@inline function _hydrostatic_vapor_density(rho, z)
    RealT = typeof(rho)
    return rho * convert(RealT, 2.0e-3) * exp(-z / convert(RealT, 5000))
end

@doc raw"""
    initial_condition_isothermal(x, t,
                                 equations::CompressibleRainyEulerVectorInvariantEquations2D)

Moist atmosphere at rest in hydrostatic balance with constant virtual temperature
``T_v = \theta_v \Pi = 280\,\mathrm{K}``, i.e., with
``p = p_0 \exp(-g z / (R_d T_v))`` and ``\varrho = p / (R_d T_v)``. The specific
humidity decays exponentially with height, and there is no cloud water or rain.
It is preserved by [`flux_nonconservative_isothermal`](@ref).
"""
function initial_condition_isothermal(x, t,
                                      equations::CompressibleRainyEulerVectorInvariantEquations2D)
    RealT = eltype(x)
    @unpack p_0, R_d, g = equations
    z = x[2]
    T_v = convert(RealT, 280)

    p = p_0 * exp(-g * z / (R_d * T_v))
    rho = p / (R_d * T_v)
    rho_vapor = _hydrostatic_vapor_density(rho, z)

    return prim2cons(SVector(rho, 0, 0, p, rho_vapor, 0, 0, g * z), equations)
end

@doc raw"""
    initial_condition_adiabatic(x, t,
                                equations::CompressibleRainyEulerVectorInvariantEquations2D)

Moist atmosphere at rest in hydrostatic balance with constant virtual potential temperature
``\theta_v = 300\,\mathrm{K}``, i.e., with ``\Pi = 1 - g z / (c_{pd} \theta_v)`` and
``\varrho = p / (R_d \theta_v \Pi)``. The specific humidity decays exponentially with
height, and there is no cloud water or rain. It is preserved by
[`flux_nonconservative_isothermal`](@ref) and [`flux_nonconservative_isentropic`](@ref).
"""
function initial_condition_adiabatic(x, t,
                                     equations::CompressibleRainyEulerVectorInvariantEquations2D)
    RealT = eltype(x)
    @unpack p_0, c_pd, R_d, g = equations
    z = x[2]
    theta_v = convert(RealT, 300)

    exner = 1 - g * z / (c_pd * theta_v)
    p = p_0 * exner^(c_pd / R_d)
    rho = p / (R_d * theta_v * exner)
    rho_vapor = _hydrostatic_vapor_density(rho, z)

    return prim2cons(SVector(rho, 0, 0, p, rho_vapor, 0, 0, g * z), equations)
end
end # @muladd
