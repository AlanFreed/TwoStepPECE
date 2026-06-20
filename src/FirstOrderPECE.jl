#=
-------------------------------------------------------------------------------

This file solves a pair of first-order ODEs:
    x′ = fₓ(model, t, x, y)        controls
    y′ = fᵧ(model, t, x, y)        responses
subject to initial conditions of x₀ and y₀ associated with time t₀ so that
    x′₀ = fₓ(model, t₀, x₀, y₀)    controls
    y′₀ = fᵧ(model, t₀, x₀, y₀)    responses
where model is an implementation of abstract type Model.  It is used to pass
model parameters or constants present in the ODEs fₓ and fᵧ.  Scalar t is the
independent variable, typically denoting time.  While x is a vector-valued
collection of control variables, and y is a vector-valued collection of their
responses, whose derivatives x′ and y′ denote dx/dt and dy/dt, respectively.  
The control function fₓ is often simpler than the response function fᵧ, but
this need not be the case.  Time t₀ is usually taken to be 0.

A local solution advances along a sub-grid with a local step size h that is
finer than the global step size dt in which h embeds.  The solution spans an
interval [t₀, T] with dt = (T - t₀) / N such that there will be N nodes of
integration whereat solutions are sought.  The global step dt is taken to be
uniform over the entire span of integration, while the local step h dynamically
adjusts to ensure that truncation error remains less than an user specified 
error tolerance denoted as tol.  Global nodes increment as: n = 0, 1, 2, ⋯, N;
while local nodes decrement as: s = S, S-1, S-2, ⋯, 0.

Usually, the error will be driven by the response y, but this is not certain 
nor necessary, so the errors arising from solving control x and response y are
both computed, with the larger of the two being sent to the PI controller.

Local truncation error comes from taking a difference between corrected and
predicted values, i.e.,
    errorₓ = ∥x_corr - x_pred∥
    errorᵧ = ∥y_corr - y_pred∥
where the local truncation error for Heun's PECE method is
    errorₓ = (h/2)∥x′_pred - x′₀∥
    errorᵧ = (h/2)∥y′_pred - y′₀∥
while for Freed's PECE method, which is based upon Gear's BDF2 method, 
the local truncation error is determined to be
    errorₓ = (2h/3)∥x′_pred - 2x′_curr + x′_prev∥
    errorᵧ = (2h/3)∥y′_pred - 2y′_curr + y′_prev∥
with the controlling error sent to the PI controller being
    εₓnext = errorₓ / max(1, ∥x_next∥)
    εᵧnext = errorᵧ / max(1, ∥y_next∥)
    error = max(εₓnext, εᵧnext)
where, e.g., ∥x∥ is a norm for vector x. 

To provide an estimate for the initial step size h to be used when taking 
the first integration step, begin with the initial condition and assign
    h₀ = ∥y₀∥ / ∥y′₀∥
where, e.g., ∥y∥ is a norm for y. To help avoid a potential wind-down or 
wind-up instability from occurring, constrain this interval so that 
dt/100 < h₀ < dt/10, and then integrate
    t₁ = t₀ + h₀
    x₁ = x₀ + h₀x′₀
    y₁ = y₀ + h₀y′₀
   x′₁ = fₓ(model, t₁, x₁, y₁)
   y′₁ = fᵧ(model, t₁, x₁, y₁)
    x₁ ← x₀ + (h₀/2)(x′₁ + x′₀)
    y₁ ← y₀ + (h₀/2)(y′₁ + y′₀)
   x′₁ ← fₓ(model, t₁, x₁, y₁)
   y′₁ ← fᵧ(model, t₁, x₁, y₁)
which is the PECE method of Heun applied to both the controls and responses. 
Afterwords, refine this estimate for an initial h  according to the formula
    h₁ = 2|(∥y₁∥-∥y₀∥) / (∥y′₁∥+∥y′₀∥)|
now constrained by dt/1000 < h₁ to help avoid potential wind-down instability.
The local step required to advance towards the first global node comes from
    S = max(2, round(dt/h₁))
from which the initial, local, step size is determined to be
    h ← dt / S
where dt advances the independent variable from global step to global step.

A predict/evaluate/correct/evaluate (PECE) scheme has been implemented (cf., 
Freed, 2017). Because this method is a two-step integrator, it is not self 
starting, so a one-step method is needed to take the first integration step.

n = 0
s = S
repeat
    if n == 0 and s == S then
        Advance the independent variable
            t₁ = t₀ + h
        after which the forward Euler method is integrated as a predictor
            x₁ = x₀ + hx′₀
            y₁ = y₀ + hy′₀
        followed by a first approximation for their rates
            x′₁ = fₓ(model, t₁, x₁, y₁)
            y′₁ = fᵧ(model, t₁, x₁, y₁)
        saving
            x′_pred = x′₁
            y′_pred = y′₁
        for use when computing error. The trapezoidal rule is its corrector
            x₁ ← x₀ + (h/2)(x′₁ + x′₀)  
            y₁ ← y₀ + (h/2)(y′₁ + y′₀)  
        after which a refined approximation for its rate is re-evaluated
            x′₁ ← fₓ(model, t₁, x₁, y₁)
            y′₁ ← fᵧ(model, t₁, x₁, y₁)
        whose local truncation error advances as
            ε_curr ← 1  
            errorₓ = (h/2)∥x′_pred - x′₀∥
            errorᵧ = (h/2)∥y′_pred - y′₀∥
            εₓnext = errorₓ / max(1, ∥x₁∥)
            εᵧnext = errorᵧ / max(1, ∥y₁∥)
            ε_next ← max(εₓnext, εᵧnext)
        Upon completion of a first step, assign
            t_prev  ← t₀
            x_prev  ← x₀     
            y_prev  ← y₀   
            x′_prev ← x′₀ 
            y′_prev ← y′₀     
            t_curr  ← t₁ 
            x_curr  ← x₁ 
            y_curr  ← y₁ 
            x′_curr ← x′₁
            y′_curr ← y′₁
        This is Heun's method, which is second-order accurate.
    else 
        Advance the independent variable
            t_next = t_curr + h
        The main PECE solver begins with the predictor
            x_pred = (1/3)(4x_curr - x_prev) + (2h/3)(2x′_curr - x′_prev)
            y_pred = (1/3)(4y_curr - y_prev) + (2h/3)(2y′_curr - y′_prev) 
        with a first approximation for their rates then being evaluated as
            x′_pred = fₓ(model, t_next, x_pred, y_pred)
            y′_pred = fᵧ(model, t_next, x_pred, y_pred)
        saving x′_pred and y′_pred for use when computing error. Gear's BDF2
        method is the corrector
            x_next ← (1/3)(4x_curr - x_prev) + (2h/3)x′_pred   
            y_next ← (1/3)(4y_curr - y_prev) + (2h/3)y′_pred   
        after which refined approximations for their ratee are re-evaluated via
            x′_next ← fₓ(model, t_next, x_next, y_next)
            y′_next ← fᵧ(model, t_next, x_next, y_next)
        whose local truncation error advances as
            errorₓ = (2h/3)∥x′_pred - 2x′_curr + x′_prev∥
            errorᵧ = (2h/3)∥y′_pred - 2y′_curr + y′_prev∥
            εₓnext = errorₓ / max(1, ∥x_next∥)
            εᵧnext = errorᵧ / max(1, ∥y_next∥)
            ε_next ← max(εₓnext, εᵧnext)
        This integrator is second-order accurate. Most importantly, the 
        corrector, is Gear's BDF2 method which is A stable.
    end
    
    A PI controller adjusts the local step size h according to the scheme: 
    if ε_curr < tol and ε_next < tol then  
        use the PI controller:
            C = (tol/ε_next)^(0.3/3) (ε_curr/ε_next)^(0.4/3)
    else
        use the I controller:
            C = √(tol/ε_next)
    end
    
    Manage the history by advancing the counters and variables:
    if C > 2 and s > 4 with s mod 2 == 1 then
        t_curr  ← t_next   
        x_curr  ← x_next   
        y_curr  ← y_next  
        x′_curr ← x′_next 
        y′_curr ← y′_next  
        ε_curr  ← ε_next
        h ← 2h 
        s ← (s - 1) ÷ 2
    else if C > 1 then
        t_prev  ← t_curr 
        x_prev  ← x_curr  
        y_prev  ← y_curr      
        x′_prev ← x′_curr    
        y′_prev ← y′_curr   
        t_curr  ← t_next   
        x_curr  ← x_next      
        y_curr  ← y_next    
        x′_curr ← x′_next
        y′_curr ← y′_next
        ε_curr  ← ε_next
        s ← s - 1
    else if ε_next < tol then
        t_prev  ← (1/2)(t_next + t_curr)
        x_prev  ← (1/2)(x_next + x_curr) - (h/8)(x′_next - x′_curr)
        y_prev  ← (1/2)(y_next + y_curr) - (h/8)(y′_next - y′_curr)
        x′_prev ← fₓ(model, t_prev, x_prev, y_prev)
        y′_prev ← fᵧ(model, t_prev, x_prev, y_prev)
        t_curr  ← t_next     
        x_curr  ← x_next     
        y_curr  ← y_next  
        x′_curr ← x′_next
        y′_curr ← y′_next
        ε_curr  ← (1/2)*(ε_next + ε_curr)
        h ← h/2
        s ← 2(s - 1)
    else  
        t_prev  ← (1/2)(t_curr + t_prev)  
        x_prev  ← (1/2)(x_curr + x_prev) - (h/8)(x′_curr - x′_prev)
        y_prev  ← (1/2)(y_curr + y_prev) - (h/8)(y′_curr - y′_prev)
        x′_prev ← fₓ(model, t_prev, x_prev, y_prev)
        y′_prev ← fᵧ(model, t_prev, x_prev, y_prev)
        ε_curr  ← 1
        h ← h/2
        s ← 2s
        repeat the integration step with half the step size
    end
    if s == 0 then
        n ← n + 1
        s ← round(dx/h)
    end
until n > N

References:
    G. Sőderlind, "Automatic control and adaptive time-stepping", Numerical
    Algorithms, Vol. 31, 2002, 281-310.

    A. D. Freed, "A Technical Note: Two-Step PECE Methods for Approximating
    Solutions to First- and Second-Order ODEs", arXiv:1707.02125 [cs.{NA}], 
    2017.
    
    A. D. Freed and I. Iskovitz, "Development and Applications of a Rosenbrock
    Integrator," NASA TM 4709, 1996.
-------------------------------------------------------------------------------   
=#

struct FirstOrderPECE <: PECE
    fₓ::Function            # Differential equation governing the controls.
    fᵧ::Function            # Differential equation governing the responses.
    model::Model            # Model constants for model being integrated.
    N::Integer              # Number of global steps to be integrated over.
    dt::Real                # Step size for for the global integrator.
    h::PF.MReal             # Current step size for the local integrator.
    n::PF.MInteger          # Global step counter, n increments to N.
    s::PF.MInteger          # Local step counter, s decrements to 0.
    t₀::Real                # Initial value for the independent variable.
    x₀::Vector{<:Real}      # Initial condition for the control variables.
    y₀::Vector{<:Real}      # Initial condition for the response variables.
    t_prev::PF.MReal        # Previous value for the independent variable.
    t_curr::PF.MReal        # Current value for the independent variable.
    t_next::PF.MReal        # Next value for the independent variable.
    x_prev::PF.MVector      # Previous values for the control variables.
    x_curr::PF.MVector      # Current values for the control variables.
    x_next::PF.MVector      # Next values for the control variables.
    y_prev::PF.MVector      # Previous values for the response variables.
    y_curr::PF.MVector      # Current values for the response variables.
    y_next::PF.MVector      # Next values for the response variables.
    x′_prev::PF.MVector     # Previous values for the derivatives dx/dt.
    x′_curr::PF.MVector     # Current values for the derivatives dx/dt.
    x′_next::PF.MVector     # Next values for the derivatives dx/dt.
    y′_prev::PF.MVector     # Previous values for the derivatives dy/dt.
    y′_curr::PF.MVector     # Current values for the derivatives dy/dt.
    y′_next::PF.MVector     # Next values for the derivatives dy/dt.
    tol::Real               # Error tolerance targetted by the PI controller.
    ε_curr::PF.MReal        # Truncation error at the current step.
    ε_next::PF.MReal        # Truncation error at the next step.
    steps::PF.MInteger      # Counter for successful steps taken.
    doubled::PF.MInteger    # Counter for times where step size was doubled.
    halved::PF.MInteger     # Counter for times where step size was halved.
    repeats::PF.MInteger    # Counter for times where a step was repeated.
    atNode::PF.MBoolean     # True if local step coincides with a global step.

    # internal constructors 
    
    function FirstOrderPECE(my_fₓ::Function,       # control ODE
                            my_fᵧ::Function,       # response ODE
                            my_model::Model,       # constants for these ODEs
                            N::Integer,            # global steps to take
                            t₀::Real,              # solution begins at
                            t_N::Real,             # solution ends at
                            x₀::Vector{<:Real},    # initial condition for x
                            y₀::Vector{<:Real},    # initial condition for y
                            tol::Real)             # error tolerance
        # verify inputs
        if N < 2
            N = convert(UInt32, 2)
        else
            N = convert(UInt32, N)
        end
        if t_N > t₀
            dx = convert(Float64, (t_N-t₀)/N)
        else
            error("Final value t_N must be greater than initial value t₀.")
        end
        if tol < 1.0e-8
            tol = 1.0e-8
        end
        if tol > 1.0e-1
            tol = 1.0e-1
        end
        
        # The initial conditions.
        if !(t₀ isa Float64)
            t₀ = convert(Float64, t₀)
        end
        if eltype(x₀) ≠ Float64
            len = length(x₀)
            x_0 = Vector{Float64}(undef, len)
            for i in 1:len
                x_0[i] = convert(Float64, x₀[i])
            end
            x₀ = x_0
        end
        if eltype(y₀) ≠ Float64
            len = length(y₀)
            y_0 = Vector{Float64}(undef, len)
            for i in 1:len
                y_0[i] = convert(Float64, y₀[i])
            end
            y₀ = y_0
        end
        
        if 4 ≠ (only(methods(my_fₓ)).nargs - 1)
            error("Function my_fₓ has four arguments: model, t, x and y.")
        end
        if 4 ≠ (only(methods(my_fᵧ)).nargs - 1)
            error("Function my_fᵧ has four arguments: model, t, x and y.")
        end
        
        x′₀ = my_fₓ(my_model, t₀, x₀, y₀)
        if x′₀ isa Vector{<:Real}
            if length(x′₀) ≠ length(x₀)
                msg = "Vector lengths for x and x′ = fₓ(model,t,x,y) differ."
                throw(DimensionMismatch(msg))
            end
        else
            error("The assigned ODE does not return a vector x′ of reals.")
        end
        
        y′₀ = my_fᵧ(my_model, t₀, x₀, y₀)
        if y′₀ isa Vector{<:Real}
            if length(y′₀) ≠ length(y₀)
                msg = "Vector lengths for y and y′ = fᵧ(model,t,x,y) differ."
                throw(DimensionMismatch(msg))
            end
        else
            error("The assigned ODE does not return a vector y′ of reals.")
        end
        
        # Establish the global step-size dt.
        dt = (t_N - t₀) / N
        
        # Determine an initial approximate for the local step-size h₀.
        # Use only the response variables to determine this estimate.
        norm_y₀  = LA.norm(y₀)
        norm_y′₀ = LA.norm(y′₀)
        if norm_y′₀ ≈ 0.0
            h₀ = dt
        else
            h₀ = norm_y₀ / norm_y′₀
        end
        # To assist in avoiding a wind-down or a wind-up instability.
        if h₀ < dt/100
            h₀ = dt / 100
        end
        if h₀ > dt/10
            h₀ = dt / 10
        end
        
        # Take a first step with this approximated, local, step-size h₀.
        t₁  = t₀ + h₀
        x₁  = x₀ + h₀*x′₀
        y₁  = y₀ + h₀*y′₀
        x′₁ = my_fₓ(my_model, t₁, x₁, y₁)
        y′₁ = my_fᵧ(my_model, t₁, x₁, y₁)
        
        x₁  = x₀ + (h₀/2)*(x′₁ + x′₀)
        y₁  = y₀ + (h₀/2)*(y′₁ + y′₀)
        x′₁ = my_fₓ(my_model, t₁, x₁, y₁)
        y′₁ = my_fᵧ(my_model, t₁, x₁, y₁)
        
        norm_y₁  = LA.norm(y₁)
        norm_y′₁ = LA.norm(y′₁)
        if norm_y′₀ + norm_y′₁ ≈ 0.0
            h₁ = h₀
        else
            h₁ = 2abs((norm_y₁ - norm_y₀) / (norm_y′₁ + norm_y′₀))
        end
        # To assist in avoiding a wind-down instability.
        if h₁ < dt/1000
            h₁ = dt / 1000
        end
        
        # Establish the initial, local, step-size h.
        S = Int(max(2, round(dt/h₁)))
        h = dt / S
       
        # Take a first step using the actual step size h.
        t₁  = t₀ + h
        x₁  = x₀ + h*x′₀
        y₁  = y₀ + h*y′₀
        x′₁ = my_fₓ(my_model, t₁, x₁, y₁)
        y′₁ = my_fᵧ(my_model, t₁, x₁, y₁)
        
        # Determine the local truncation errors.
        εₓ  = (h/2)*LA.norm(x′₁ - x′₀)
        εᵧ  = (h/2)*LA.norm(y′₁ - y′₀)
        
        # Finish integration with the corrector.
        x₁  = x₀ + (h/2)*(x′₁ + x′₀)
        y₁  = y₀ + (h/2)*(y′₁ + y′₀)
        x′₁ = my_fₓ(my_model, t₁, x₁, y₁)
        y′₁ = my_fᵧ(my_model, t₁, x₁, y₁)
        
        # Assign the history variables.
        t_prev  = PF.MReal(t₀)
        t_curr  = PF.MReal(t₀+h)
        t_next  = PF.MReal(0)
        x_prev  = PF.MVector(x₀)
        x_curr  = PF.MVector(x₁)
        x_next  = PF.MVector(length(x₀))
        y_prev  = PF.MVector(y₀)
        y_curr  = PF.MVector(y₁)
        y_next  = PF.MVector(length(y₀))
        x′_prev = PF.MVector(x′₀)
        x′_curr = PF.MVector(x′₁)
        x′_next = PF.MVector(length(x₀))
        y′_prev = PF.MVector(y′₀)
        y′_curr = PF.MVector(y′₁)
        y′_next = PF.MVector(length(y₀))
        
        # Assign the counters and step size.
        n = PF.MInteger(0)
        s = PF.MInteger(S-1)
        h = PF.MReal(h)
        
        # Assign truncation errors.
        norm_x₁ = LA.norm(x₁)
        norm_y₁ = LA.norm(y₁)
        ε_curr  = PF.MReal(1)
        ε_next  = PF.MReal(max(εₓ/max(1,norm_x₁), εᵧ/max(1,norm_y₁)))
        
        # Create integration counters.
        steps   = PF.MInteger(1)
        doubled = PF.MInteger(0)
        halved  = PF.MInteger(0)
        repeats = PF.MInteger(0)
        atNode  = PF.MBoolean(false)
        
        print("\n.")
        new(my_fₓ, my_fᵧ, my_model, N, dt, h, n, s, t₀, x₀, y₀, 
            t_prev, t_curr, t_next, x_prev, x_curr, x_next, 
            y_prev, y_curr, y_next, x′_prev, x′_curr, x′_next, 
            y′_prev, y′_curr, y′_next,
            tol, ε_curr, ε_next, steps, doubled, halved, repeats, atNode)
    end 
    
    # basic constructor
    
    function FirstOrderPECE(fₓ::Function,
                            fᵧ::Function,
                            model::Model,
                            N::Integer,
                            dt::Real,
                            h::PF.MReal,
                            n::PF.MInteger,
                            s::PF.MInteger,
                            t₀::Real,
                            x₀::Vector{<:Real},
                            y₀::Vector{<:Real},
                            t_prev::PF.MReal,
                            t_curr::PF.MReal,
                            t_next::PF.MReal,
                            x_prev::PF.MVector,
                            x_curr::PF.MVector,
                            x_next::PF.MVector,
                            y_prev::PF.MVector,
                            y_curr::PF.MVector,
                            y_next::PF.MVector,
                            x′_prev::PF.MVector,
                            x′_curr::PF.MVector,
                            x′_next::PF.MVector,
                            y′_prev::PF.MVector,
                            y′_curr::PF.MVector,
                            y′_next::PF.MVector,
                            tol::Real,
                            ε_curr::PF.MReal,
                            ε_next::PF.MReal,
                            steps::PF.MInteger,
                            doubled::PF.MInteger,
                            halved::PF.MInteger,
                            repeats::PF.MInteger,
                            atNode::PF.MBoolean)
        
        new(fₓ, fᵧ, model, N, dt, h, n, s, t₀, x₀, y₀, 
            t_prev, t_curr, t_next, x_prev, x_curr, x_next, 
            y_prev, y_curr, y_next, x′_prev, x′_curr, x′_next, 
            y′_prev, y′_curr, y′_next,
            tol, ε_curr, ε_next, steps, doubled, halved, repeats, atNode)
    end
end # FirstOrderPECE

function advance!(pece::FirstOrderPECE)
    if pece.n > pece.N
        print("\nThe ODE has been solved.\n")
        return nothing
    end
    PF.set!(pece.atNode, false)
    
    # get history variables
    h       = PF.get(pece.h)
    t_prev  = PF.get(pece.t_prev)
    x_prev  = PF.Vector(pece.x_prev)
    y_prev  = PF.Vector(pece.y_prev)
    x′_prev = PF.Vector(pece.x′_prev)
    y′_prev = PF.Vector(pece.y′_prev)
    t_curr  = PF.get(pece.t_curr)
    x_curr  = PF.Vector(pece.x_curr)
    y_curr  = PF.Vector(pece.y_curr)
    x′_curr = PF.Vector(pece.x′_curr)
    y′_curr = PF.Vector(pece.y′_curr)
    ε_curr  = PF.get(pece.ε_curr)
    
    # advance the independent variable
    t_next  = t_curr + h
    
    # P: apply Freed's predictor that pairs with Gear's BDF2 formula
    x_pred  = (1/3)*(4x_curr - x_prev) + (2h/3)*(2x′_curr - x′_prev)
    y_pred  = (1/3)*(4y_curr - y_prev) + (2h/3)*(2y′_curr - y′_prev)
    # E: evaluate the ODEs
    x′_pred = pece.fₓ(pece.model, t_next, x_pred, y_pred)
    y′_pred = pece.fᵧ(pece.model, t_next, x_pred, y_pred)
    # C: apply Gear's BDF2 formula as the corrector
    x_next  = (1/3)*(4x_curr - x_prev) + (2h/3)*x′_pred
    y_next  = (1/3)*(4y_curr - y_prev) + (2h/3)*y′_pred
    # E: re-evaluate the ODEs
    x′_next = pece.fₓ(pece.model, t_next, x_next, y_next)
    y′_next = pece.fᵧ(pece.model, t_next, x_next, y_next)
    
    # determine the local truncation error
    errorₓ = (2h/3)*LA.norm(x′_pred - 2x′_curr + x′_prev)
    εₓnext = errorₓ / max(1, LA.norm(x_next))
    errorᵧ = (2h/3)*LA.norm(y′_pred - 2y′_curr + y′_prev)
    εᵧnext = errorᵧ / max(1, LA.norm(y_next))
    ε_next = max(εₓnext, εᵧnext)
    
    # apply the step controller
    if (ε_curr < pece.tol) && (ε_next < pece.tol) 
        # use a PI controller:
        C = (pece.tol/ε_next)^(0.1) * (ε_curr/ε_next)^(0.4/3)
    else
        # use an I controller:
        C = sqrt(pece.tol/ε_next)
    end
    
    # Manage the history by advancing its counters and variables
    if (C > 2) && (pece.s > 4) && (pece.s%2 == 1)
        PF.set!(pece.t_curr, t_next)
        for i in 1:pece.x_curr.len
            pece.x_curr[i]  = x_next[i]
            pece.x′_curr[i] = x′_next[i]
        end
        for i in 1:pece.y_curr.len
            pece.y_curr[i]  = y_next[i]
            pece.y′_curr[i] = y′_next[i]
        end
        PF.set!(pece.ε_curr, ε_next)
        PF.set!(pece.steps, PF.get(pece.steps)+1)
        PF.set!(pece.doubled, PF.get(pece.doubled)+1)
        PF.set!(pece.h, 2PF.get(pece.h))
        PF.set!(pece.s, (PF.get(pece.s)-1)÷2)
    elseif C > 1
        PF.set!(pece.t_prev, t_curr)
        for i in 1:pece.x_curr.len
            pece.x_prev[i]  = x_curr[i]
            pece.x′_prev[i] = x′_curr[i]
        end
        for i in 1:pece.y_curr.len
            pece.y_prev[i]  = y_curr[i]
            pece.y′_prev[i] = y′_curr[i]
        end
        PF.set!(pece.t_curr, t_next)
        for i in 1:pece.x_curr.len
            pece.x_curr[i]  = x_next[i]
            pece.x′_curr[i] = x′_next[i]
        end
        for i in 1:pece.y_curr.len
            pece.y_curr[i]  = y_next[i]
            pece.y′_curr[i] = y′_next[i]
        end
        PF.set!(pece.ε_curr, ε_next)
        PF.set!(pece.steps, PF.get(pece.steps)+1)
        PF.set!(pece.s, PF.get(pece.s)-1)
    elseif ε_next < pece.tol    # with C ≤ 1
        PF.set!(pece.t_prev, (1/2)*(t_next+t_curr))
        for i in 1:pece.x_curr.len
            pece.x_prev[i] = ((1/2)*(x_next[i] + x_curr[i])
                              - (h/8)*(x′_next[i] - x′_curr[i]))
        end
        for i in 1:pece.y_curr.len
            pece.y_prev[i] = ((1/2)*(y_next[i] + y_curr[i])
                              - (h/8)*(y′_next[i] - y′_curr[i]))
        end
        t_prev  = PF.get(pece.t_prev)
        x_prev  = PF.Vector(pece.x_prev)
        y_prev  = PF.Vector(pece.y_prev)
        x′_prev = pece.fₓ(pece.model, t_prev, x_prev, y_prev)
        for i in 1:pece.x′_prev.len
            pece.x′_prev[i] = x′_prev[i]
        end
        y′_prev = pece.fᵧ(pece.model, t_prev, x_prev, y_prev)
        for i in 1:pece.y′_prev.len
            pece.y′_prev[i] = y′_prev[i]
        end
        PF.set!(pece.t_curr, t_next)
        for i in 1:pece.x_curr.len
            pece.x_curr[i]  = x_next[i]
            pece.x′_curr[i] = x′_next[i]
        end
        for i in 1:pece.y_curr.len
            pece.y_curr[i]  = y_next[i]
            pece.y′_curr[i] = y′_next[i]
        end
        PF.set!(pece.ε_curr, (1/2)*(ε_next+ε_curr))
        PF.set!(pece.steps, PF.get(pece.steps)+1)
        PF.set!(pece.halved, PF.get(pece.halved)+1)
        PF.set!(pece.h, PF.get(pece.h)/2)
        PF.set!(pece.s, 2(PF.get(pece.s)-1))
    else    # ε_next ≥ tol and C ≤ 1
        PF.set!(pece.t_prev, (1/2)*(t_curr+t_prev))
        for i in 1:pece.x_curr.len
            pece.x_prev[i] = ((1/2)*(x_curr[i] + x_prev[i])
                              - (h/8)*(x′_curr[i] - x′_prev[i]))
        end
        for i in 1:pece.y_curr.len
            pece.y_prev[i] = ((1/2)*(y_curr[i] + y_prev[i])
                              - (h/8)*(y′_curr[i] - y′_prev[i]))
        end
        t_prev  = PF.get(pece.t_prev)
        x_prev  = PF.Vector(pece.x_prev)
        y_prev  = PF.Vector(pece.y_prev)
        x′_prev = pece.fₓ(pece.model, t_prev, x_prev, y_prev)
        for i in 1:pece.x′_prev.len
            pece.x′_prev[i] = x′_prev[i]
        end
        y′_prev = pece.fᵧ(pece.model, t_prev, x_prev, y_prev)
        for i in 1:pece.y′_prev.len
            pece.y′_prev[i] = y′_prev[i]
        end
        PF.set!(pece.ε_curr, 1)    # forces the I controller to be used
        PF.set!(pece.repeats, PF.get(pece.repeats)+1)
        PF.set!(pece.h, PF.get(pece.h)/2)
        PF.set!(pece.s, 2PF.get(pece.s))
        advance!(pece)             # repeat this integration step
    end
    if pece.s == 0
        PF.set!(pece.n, PF.get(pece.n)+1)
        PF.set!(pece.s, Int(PF.round(pece.dt/pece.h)))
        counts = Int(max(1, pece.N÷50))
        if pece.n % counts == 0
            print(".")    # prints out about 50  .  at completion
        end
        if pece.n > pece.N
            print("\nThe ODE has been solved.\n")
        end
        PF.set!(pece.atNode, true)
    end
    return nothing
end # advance!
    