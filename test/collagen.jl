"""
-------------------------------------------------------------------------------

    author:  Alan Freed
    date:    26 May 2026
    updated: 26 May 2026

    This example describes the response of a biologic fiber under load, where a
    change in strain dε and a change in entropy dη is caused by changes in 
    temperature dθ and stress dσ via
    
        dε = α dθ + (1/E₁ + ((1-½γε₂²) / (E₂ + (β+γε₂)(σ-σ₀)))) dσ
        dη = C dθ/θ - α/ρ dσ
    
    with hidden or internal variables of 
    
        ε₁ = α(θ-θ₀) + (σ-σ₀)/E
        ε₂ = ε - ε₁
    
    where ε₁ represents a strain caused by a stretching of molecular bonds
    while ε₂ represents a strain caused by a molecular reconfiguration, 
    whose sum describes the observed strain ε. 
    
    Here θ₀ and θ are the initial and current temperatures, in an absolute sense,
    with σ₀ and σ being the initial and current stresses carried by the fiber.
    
    The model is described in terms of the following material constants:
    
        ρ   mass density
        C   specific heat
        E₁  elastic modulus due to molecular stretching (Young's modulus)
        E₂  elastic modulus due to molecular reconfiguration
        α   thermal strain coefficient
        β   exponential growth in σ with increasing ε₂  (Fung's parameter)
        γ   exponential growth in σ with increasing ε₂²
    
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
    
end # collagen