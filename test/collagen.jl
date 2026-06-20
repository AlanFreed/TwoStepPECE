"""
-------------------------------------------------------------------------------

    author:  Alan Freed
    date:    26 May 2026
    updated: 27 May 2026

    This example describes the response of a biologic fiber under load, where a
    change in strain dε and a change in entropy dη are caused by changes in 
    temperature dθ and stress dσ via
    
        dε = α dθ + (1/E₁ + ((1-½γε₂²) / (E₂ + (β+γε₂)(σ-σ₀)))) dσ
        dη = C dθ/θ - α/ρ dσ
    
    with hidden or internal variables of 
    
        ε₁ = α(θ-θ₀) + (σ-σ₀)/E
        ε₂ = ε - ε₁
    
    where ε₁ represents a strain caused by a stretching of molecular bonds
    while ε₂ represents a strain caused by a molecular reconfiguration, 
    whose sum describes the observed strain ε. 
    
    θ₀ and θ are the initial and current temperatures, in an absolute sense.
    σ₀ and σ are the initial and current stresses carried by the fiber.
    
    The model is described in terms of the following material constants:
    
        ρ   mass density
        C   specific heat
        E₁  elastic modulus due to molecular stretching (Young's modulus)
        E₂  elastic modulus due to molecular reconfiguration
        α   thermal strain coefficient
        β   exponential growth in σ with increasing ε₂  (Fung's parameter)
        γ   exponential growth in σ with increasing ε₂² (strain limits ε₂)
    
"""
module collagen

using
    Measures,        # needed to pad white space (margins) around a plot
    PhysicalFields,
    Plots,
    TwoStepPECE
    
import
    PhysicalFields as PF
    
export
    run
    
function model(par::Vector{<:Real}, x::Real, y::Vector{<:Real})::Tuple{Vector{<:Real}, Vector{<:Real}}

    # model
    mod = par[1]
    
    # control variables
    ctrl = par[2]
    if ctrl ≈ 1.0
        # stress control
        σ_prev = y[1]
        σ_curr = y[2]
        σ_next = x
        ε_prev = y[3]
        ε_curr = y[4]
        θ_prev = y[5]
        θ_curr = y[6]
    else # ctrl ≈ 2.0
        # strain control
        σ_prev = y[1]
        σ_curr = y[2]
        ε_prev = y[3]
        ε_curr = y[4]
        ε_next = x
        θ_prev = y[5]
        θ_curr = y[6]
    end
    
    # physical constants applicable to all models
    ρ = par[2] # mass density
    C = par[3] # specific heat
    α = par[4] # thermal strain coefficient

    # model specific material constants
    if mod ≈ 1.0
        # Hooke's model
        E = par[5] # Young's modulus
    elseif mod ≈ 2.0
        # Fung's model
        E = par[5] # elastic tangent modulus at σ = 0
        β = par[6] # exponential growth in σ with increasing ε
    else # mod ≈ 3.0
        # Freed's model
        E₁ = par[5] # elastic modulus due to molecular stretching
        E₂ = par[6] # elastic modulus due to molecular reconfiguration
        β  = par[8] # exponential growth in σ with increasing ε₂
        γ  = par[9] # exponential growth in σ with increasing ε₂²
    end
    
    # initial conditions
    θ₀ = c[1]
    σ₀ = c[2]
    
    # independent variables
    θ  = x[1]
    σ  = x[2]
    dθ = x[3]
    dσ = x[4]

    # differential equations for dependent variables
    dη = C*dθ/θ - (α/ρ)*dσ
    dε = α*dθ + (1/E₁ + ((1-0.5γ*ε₂*ε₂) / (E₂ + (β+γ*ε₂)*(σ-σ₀))))*dσ
    
    ode = Vector{Float64}(undef, 2)
    ode[1] = dε
    ode[2] = dη
    
    # internal variables
    ε₁ = α*(θ-θ₀) + (σ-σ₀)/E
    ε₂ = ε - ε₁

    z = Vector{Float64}(undef, 2)
    z[1] = ε₁
    z[2] = ε₂
    
    return ode, z
end # model

function run()

end # run
end # collagen