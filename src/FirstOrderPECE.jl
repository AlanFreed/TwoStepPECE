#=
-------------------------------------------------------------------------------

This file solves the following system of first-order ODEs:
    y′  = ode(x, y)
subject to an initial condition of y₀ associated with x₀ so that
    y′₀ = ode(x₀, y₀)
where x is a scalar-valued independent variable, and where y is a vector-valued
dependent variable. Notation y′ denotes dy/dx.

A local solution advances along a sub-grid with a local step size h that is
finer than the global step size dx in which h embeds. The solution spans an
interval [x₀, X] with dx = (X - x₀) / N such that there will be N nodes of
integration whereat solutions are sought. Global step size dx is taken to be
uniform over the entire span of integration, while local step h dynamically
adjusts to ensure that the truncation error remains less than an user specified 
error tolerance denoted as tol. Global nodes increment as: n = 0, 1, 2, ⋯, N;
while local nodes decrement as: s = S, S-1, S-2, ⋯, 0.

A local truncation error comes from taking a difference between corrected and
predicted values, i.e.,
    error = ∥y_corr - y_pred∥
where for Heun's PECE method
    error = (h/2)*∥y′_pred - y′₀∥
and for Freed's PECE method based upon Gear's BDF2 method
    error = (2h/3)*∥y′_pred - 2y′_curr + y′_prev∥
where ∥y∥ is a norm for y. The local truncation error is then given by
    if ∥y_next∥ < 1  then 
        An absolute measure of error is used.
        ε_next = error
    else  
        A relative measure for error is used.
        ε_next = error / ∥y_next∥
    end

To provide an estimate for the initial step size h to be used when taking the
the first integration step then, from the initial condition, assign
    h₀ = ∥y₀∥ / ∥y′₀∥
where ∥y∥ is a norm for y. To help avoid a potential wind-down or a wind-up
instability, constrain this interval so that dx/100 < h₀ < dx/10 and then 
integrate
    x₁  = x₀ + h₀
    y₁  = y₀ + h₀*y′₀
    y′₁ = ode(x₁, y₁)
    y₁  ← y₀ + (h₀/2)*(y′₁ + y′₀)
    y′₁ ← ode(x₁, y₁)
which is the PECE method of Heun. Afterwords, refine this estimate for h  
according to the formula
    h = 2|[(∥y₁∥ - ∥y₀∥) / (∥y′₁∥ + ∥y′₀∥)]|
now constrained by dx/1000 < h to help avoid a potential wind-down instability.
The local steps required to advance toward the first global node comes from
    S = max(2, round(dx/h))
from which the initial, local, step size is determined to be
    h ← dx / S
where dx advances the independent variable from global step to global step.

A predict/evaluate/correct/evaluate (PECE) scheme has been implemented (see 
Freed, 2017). Because this method is a two-step integrator, it is not self 
starting, so a one-step method is needed to take the first integration step.

n = 0
s = S
repeat
    if n == 0 and s == S then
        Advance the independent variable
            x₁ = x₀ + h
        after which the forward Euler method is integrated as a predictor
            y₁ = y₀ + h*y′₀
        followed by a first approximation for its rate
            y′₁ = ode(x₁, y₁)
        saving
            y′_pred = y′₁
        for use when computing error. The trapezoidal rule is its corrector
            y₁ ← y₀ + (h/2)*(y′₁ + y′₀)  
        after which a refined approximation for its rate is re-evaluated
            y′₁ ← ode(x₁, y₁)
        whose local truncation error advances as
            error = (h/2)*∥y′_pred - y′₀∥
            ε_curr ← 1    
            ε_next ← error / max(1, ∥y₁∥)
        Upon completion of a first step, assign
            x_prev  ← x₀     
            y_prev  ← y₀    
            y′_prev ← y′₀     
            x_curr  ← x₀ + h 
            y_curr  ← y₁ 
            y′_curr ← y′₁
        This is Heun's method, which is second-order accurate.
    else 
        Advance the independent variable
            x_next = x_curr + h
        The main PECE solver begins with the predictor
            y_next  = (1/3)*(4y_curr - y_prev) + (2h/3)*(2y′_curr - y′_prev) 
        with a first approximation for its rate then being evaluated as
            y′_next = ode(x_next, y_next)
        saving
            y′_pred = y′_next
        for use when computing error. Gear's BDF2 method is the corrector
            y_next  ← (1/3)*(4y_curr - y_prev) + (2h/3)*y′_pred   
        after which a refined approximation for its rate is re-evaluated via
            y′_next ← ode(x_next, y_next)
        whose local truncation error advances as
            error  = (2h/3)*∥y′_pred - 2y′_curr + y′_prev∥
            ε_next ← error / max(1, ∥y_next∥)
        This corrector is second-order accurate. Most importantly, Gear's
        BDF2 method is A stable.
    end
    
    A PI controller adjusts the local step size h according to the scheme 
    if ε_curr < tol and ε_next < tol then  
        use the PI controller:
            C = (tol/ε_next)^(0.3/3) * (ε_curr/ε_next)^(0.4/3)
    else
        use the I controller:
            C = (tol/ε_next)^(1/2)
    end
    
    Manage the history by advancing the counters and variables:
    if C > 2 and s > 4 and s mod 2 == 1 then
        x_curr  ← x_next   
        y_curr  ← y_next    
        y′_curr ← y′_next  
        ε_curr  ← ε_next
        h ← 2h 
        s ← (s - 1) ÷ 2
    else if C > 1 then
        x_prev  ← x_curr  
        y_prev  ← y_curr    
        y′_prev ← y′_curr   
        x_curr  ← x_next     
        y_curr  ← y_next   
        y′_curr ← y′_next
        ε_curr  ← ε_next
        s ← s - 1
    else if ε_next < tol then
        x_prev  ← (1/2)*(x_next + x_curr)
        y_prev  ← (1/2)*(y_next + y_curr) - (h/8)*(y′_next - y′_curr)
        y′_prev ← (1/8)*(3y′_next + 6y′_curr - y′_prev)
        x_curr  ← x_next     
        y_curr  ← y_next    
        y′_curr ← y′_next
        ε_curr  ← (1/2)*(ε_next + ε_curr)
        h ← h/2
        s ← 2(s - 1)
    else  
        x_prev  ← (1/2)*(x_curr + x_prev)  
        y_prev  ← (1/2)*(y_curr + y_prev) - (h/8)*(y′_curr - y′_prev)
        y′_prev ← (1/2)*(y′_curr + y′_prev) 
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
    ode::Function           # The differential equation that is to be solved.
    N::UInt32               # Number of global steps to be integrated over.
    dx::Float64             # Step size for for the global integrator.
    h::PF.MReal             # Current step size for the local integrator.
    n::PF.MInteger          # Global step counter, n increments to N.
    s::PF.MInteger          # Local step counter, s decrements to 0.
    x₀::Float64             # Initial value for the independent variable.
    y₀::Vector{<:Real}      # Initial condition for the dependent variable.
    x_prev::PF.MReal        # Previous value for the independent variable.
    x_curr::PF.MReal        # Current value for the independent variable.
    x_next::PF.MReal        # Next value for the independent variable.
    y_prev::PF.MVector      # Previous value for the dependent variable.
    y_curr::PF.MVector      # Current value for the dependent variable.
    y_next::PF.MVector      # Next value for the dependent variable.
    y′_prev::PF.MVector     # Previous value for the derivative dy/dx.
    y′_curr::PF.MVector     # Current value for the derivative dy/dx.
    y′_next::PF.MVector     # Next value for the derivative dy/dx.
    tol::Float64            # Error tolerance targetted by the PI controller.
    ε_curr::PF.MReal        # Truncation error at the current step.
    ε_next::PF.MReal        # Truncation error at the next step.
    steps::PF.MInteger      # Counter for successful steps taken.
    doubled::PF.MInteger    # Counter for times where step size was doubled.
    halved::PF.MInteger     # Counter for times where step size was halved.
    repeats::PF.MInteger    # Counter for times where a step was repeated.
    atNode::PF.MBoolean     # True if local step coincides with a global step.

    # internal constructors
    
    function FirstOrderPECE(my_ode::Function,    # differential equation
                            N::Integer,          # global steps to take
                            x₀::Real,            # solution begins at
                            x_N::Real,           # solution ends at
                            y₀::Vector{<:Real},  # initial condition
                            tol::Real)           # error tolerance
        # verify inputs
        if N < 2
            N = convert(UInt32, 2)
        else
            N = convert(UInt32, N)
        end
        if x_N > x₀
            dx = convert(Float64, (x_N-x₀)/N)
        else
            error("Final value x_N must be greater than initial value x₀.")
        end
        if tol < 1.0e-8
            tol = 1.0e-8
        end
        if tol > 1.0e-1
            tol = 1.0e-1
        end
        
        # The initial conditions.
        if !(x₀ isa Float64)
            x₀ = convert(Float64, x₀)
        end
        if eltype(y₀) ≠ Float64
            len = length(y₀)
            y_0 = Vector{Float64}(undef, len)
            for i in 1:len
                y_0[i] = convert(Float64, y₀[i])
            end
            y₀ = y_0
        end
        if 2 ≠ (only(methods(my_ode)).nargs - 1)
            error("Function my_ode must have two arguments, viz., x and y.")
        end
        y′₀ = my_ode(x₀, y₀)
        if y′₀ isa Vector{<:Real}
            if length(y′₀) ≠ length(y₀)
                msg = "The length of vectors y and y′ = f(x,y) differ."
                throw(DimensionMismatch(msg))
            end
            if eltype(y′₀) ≠ Float64
                len = length(y′₀)
                y′0 = Vector{Float64}(undef, len)
                for i in 1:len
                    y′0[i] = convert(Float64, y′₀[i])
                end
                y′₀ = y′0
            end
        else
            error("The assigned ODE does not return a vector of reals.")
        end
        
        # Determine the initial, local, step size.
        norm_y₀  = LA.norm(y₀)
        norm_y′₀ = LA.norm(y′₀)
        if norm_y′₀ ≈ 0.0
            h₀ = dx / 100
        else
            h₀ = norm_y₀ / norm_y′₀
        end
        x₁  = x₀ + h₀
        y₁  = y₀ + h₀*y′₀
        y′₁ = my_ode(x₁, y₁)
        y₁  = y₀ + (h₀/2)*(y′₁ + y′₀)
        y′₁ = my_ode(x₁, y₁)
        h   = 2abs((LA.norm(y₁) - norm_y₀) / (LA.norm(y′₁) + norm_y′₀))
        if h < dx/1000
            h = dx / 1000
        end
        if h > dx/10
            h = dx / 10
        end
        S = Int(max(2, round(dx/h)))
        h = dx / S
        
        # Determine the truncation error.
        x₁  = x₀ + h
        y₁  = y₀ + h*y′₀
        y′₁ = my_ode(x₁, y₁)
        err = (h/2)*LA.norm(y′₁ - y′₀)
        
        # Finish integration with the corrector.
        y₁  = y₀ + (h/2)*(y′₁ + y′₀)
        y′₁ = my_ode(x₁, y₁)
        
        # Assign the history variables.
        x_prev  = PF.MReal(x₀)
        x_curr  = PF.MReal(x₀+h)
        x_next  = PF.MReal()
        y_prev  = PF.MVector(y₀)
        y_curr  = PF.MVector(y₁)
        y_next  = PF.MVector(length(y₀))
        y′_prev = PF.MVector(y′₀)
        y′_curr = PF.MVector(y′₁)
        y′_next = PF.MVector(length(y₀))
        
        # Assign the counters and step size.
        n = PF.MInteger(0)
        s = PF.MInteger(S-1)
        h = PF.MReal(h)
        
        # Assign truncation errors.
        ε_curr = PF.MReal(1)
        ε_next = PF.MReal(err/max(1, LA.norm(y₁)))
        
        # Create integration counters.
        steps   = PF.MInteger(1)
        doubled = PF.MInteger(0)
        halved  = PF.MInteger(0)
        repeats = PF.MInteger(0)
        atNode  = PF.MBoolean(false)
        
        print("\n.")
        new(my_ode, N, dx, h, n, s, x₀, y₀, x_prev, x_curr, x_next, 
            y_prev, y_curr, y_next, y′_prev, y′_curr, y′_next, tol, 
            ε_curr, ε_next, steps, doubled, halved, repeats, atNode)
    end 
    
    # constructor called by JSON3
    
    function FirstOrderPECE(ode::Function,
                            N::UInt32,
                            dx::Float64,
                            h::PF.MReal,
                            n::PF.MInteger,
                            s::PF.MInteger,
                            x₀::Float64,
                            y₀::Float64,
                            x_prev::PF.MReal,
                            x_curr::PF.MReal,
                            x_next::PF.MReal,
                            y_prev::PF.MVector,
                            y_curr::PF.MVector,
                            y_next::PF.MVector,
                            y′_prev::PF.MVector,
                            y′_curr::PF.MVector,
                            y′_next::PF.MVector,
                            tol::Float64,
                            ε_curr::PF.MReal,
                            ε_next::PF.MReal,
                            steps::PF.MInteger,
                            doubled::PF.MInteger,
                            halved::PF.MInteger,
                            repeats::PF.MInteger,
                            atNode::PF.MBoolean)
        
        new(my_ode, N, dx, h, n, s, x₀, y₀, x_prev, x_curr, x_next, 
            y_prev, y_curr, y_next, y′_prev, y′_curr, y′_next, tol, 
            ε_curr, ε_next, steps, doubled, halved, repeats, atNode)
    end
end # FirstOrderPECE
 
# Functions that exist for instances of FirstOrderPECE.

StructTypes.StructType(::Type{FirstOrderPECE}) = StructTypes.Struct()

function toFile(y::FirstOrderPECE, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, y)
        write(json_stream, '\n')
    else
        error("The supplied JSON stream is not open.")
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{FirstOrderPECE}, json_stream::IOStream)::FirstOrderPECE
    if isopen(json_stream)
        y = JSON3.read(readline(json_stream), FirstOrderPECE)
    else
        error("The supplied JSON stream is not open.")
    end
    return y
end

function advance!(pece::FirstOrderPECE)
    if pece.n > pece.N
        print("\nThe ODE has been solved.")
        return nothing
    end
    
    pece.atNode = false
    # advance the independent variable
    x_next = pece.x_curr + pece.h
    # apply a predictor that pairs with Gear's BDF2 formula
    y_next = (1/3)*((4pece.y_curr - pece.y_prev) 
                    + (2pece.h/3)*(2pece.y′_curr - pece.y′_prev))
    # evaluate the ODE
    y′_next = pece.ode(PF.Real(x_next), PF.Vector(y_next))
    # save this value for evaluating truncation error
    y′_pred = PF.MVector(y′_next)
    # apply Gear's BDF2 formula as the corrector:
    y_next = (1/3)*(4pece.y_curr - pece.y_prev) + (2pece.h/3)*y′_pred
    # re-evaluate the ODE:
    y′_next = pece.ode(PF.Real(x_next), PF.Vector(y_next))
    
    # determine the local truncation error:
    err = (2pece.h/3)*PF.norm(y′_pred - 2pece.y′_curr + pece.y′_prev)
    ε_next = err / max(1, PF.norm(y_next))
    
    # apply the step controller:
    if (pece.ε_curr < pece.tol) && (ε_next < pece.tol) 
        # use a PI controller:
        C = (pece.tol/ε_next)^(0.1) * (pece.ε_curr/ε_next)^(0.4/3)
    else
        # use an I controller:
        C = PF.sqrt(pece.tol/ε_next)
    end
    
    # Manage the history by advancing its counters and variables:
    if (C > 2) && (pece.s > 4) && (pece.s%2 == 1)
        pece.x_curr  = x_next   
        pece.y_curr  = y_next    
        pece.y′_curr = y′_next
        pece.ε_curr  = ε_next
        pece.steps   = pece.steps + 1
        pece.doubled = pece.doubled + 1
        pece.h = 2pece.h 
        pece.s = (pece.s - 1) ÷ 2
    elseif C > 1
        pece.x_prev  = pece.x_curr  
        pece.y_prev  = pece.y_curr    
        pece.y′_prev = pece.y′_curr   
        pece.x_curr  = x_next     
        pece.y_curr  = y_next   
        pece.y′_curr = y′_next 
        pece.ε_curr  = ε_next
        pece.steps   = pece.steps + 1
        pece.s = pece.s - 1
    elseif ε_next < pece.tol    # with C ≤ 1
        pece.x_prev  = (1/2)*(x_next + pece.x_curr)
        pece.y_prev  = (1/2)*((y_next + pece.y_curr) 
                       - (pece.h/8)*(y′_next - pece.y′_curr))
        pece.y′_prev = (1/8)*(3y′_next + 6pece.y′_curr - pece.y′_prev)
        pece.x_curr  = x_next     
        pece.y_curr  = y_next    
        pece.y′_curr = y′_next
        pece.ε_curr  = (1/2)*(ε_next + pece.ε_curr)
        pece.steps   = pece.steps + 1
        pece.halved  = pece.halved + 1
        pece.h = pece.h/2
        pece.s = 2(pece.s - 1)
    else    # ε_next ≥ tol and C ≤ 1
        pece.x_prev  = (1/2)*(pece.x_curr + pece.x_prev)  
        pece.y_prev  = (1/2)*((pece.y_curr + pece.y_prev) 
                       - (pece.h/8)*(pece.y′_curr - pece.y′_prev))
        pece.y′_prev = (1/2)*(pece.y′_curr + pece.y′_prev) 
        pece.ε_curr  = 1    # forces the I controller to be used
        pece.repeats = pece.repeats + 1
        pece.h = pece.h/2
        pece.s = 2pece.s
        advance!(pece)      # repeat this integration step
    end
    if pece.s == 0 then
        pece.n = pece.n + 1
        pece.s = Int(PF.round(pece.dx/pece.h))
        counts = Int(max(1, pece.N÷50))
        if pece.n % counts == 0
            print(".")    # prints out about 50  .  at completion
        end
        if pece.n > pece.N
            print("\nThe ODE has been solved.")
        end
        pece.atNode = true
    end
    return nothing
end # advance!
    