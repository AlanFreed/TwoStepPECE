#=
-------------------------------------------------------------------------------

This file solves the following system of second-order ODEs
    y″  = ode(x, y, y′)
subject to an initial condition of y₀ and y′₀ associated with x₀ so that
    y″₀ = ode(x₀, y₀, y′₀)
where x is a scalar-valued independent variable, and y is a vector-valued
dependent variable, with y′ being its first-order derivative. Notation y′
denotes dy/dx, while notation y″ denotes d²y/dx².

A local solution advances along a sub-grid with a local step size h that is
finer than the global step size dx in which h embeds. The solution spans an
interval [x₀, X] with dx = (X - x₀) / N such that there will be N nodes of
integration whereat solutions are sought. Global step size dx is taken to be
uniform over the entire span of integration, while local step h dynamically
adjusts to ensure that the truncation error remains less than an user specified 
error tolerance denoted as tol. Global nodes increment as: n = 0, 1, 2, ⋯, N;
while local nodes decrement as: s = S, S-1, S-2, ⋯, 0.

For a second-order ODE, local truncation error comes from taking a difference 
between corrected and predicted values for y, viz.
    error = ∥y_corr - y_pred∥
where, for the first solution step, the local truncation error is
    error = (h²/6)*∥y″_pred - y″₀∥
with the remaining integration steps having a truncation error of
    error = (h/24)*∥y′_pred + 2y′_curr - 3y′_prev
                    + (h/3)*(10y″_pred - 11y″_curr + y″_prev)∥
where ∥y∥ is a norm for y. The local truncation error is then given by
    if ∥y_next∥ < 1  then 
        An absolute measure of error is used.
        ε_next = error
    else  
        A relative measure for error is used.
        ε_next = error / ∥y_next∥
    end

To provide an estimate for the initial step size h to be used when taking
the first integration step, from initial conditions, assign
    h₀ = ∥y₀∥ / ∥y′₀∥
where ∥y∥ is a norm for y. To help avoid a potential wind-down or a wind-up
instability, constrain this interval so that dx/100 < h₀ < dx/10 and then 
integrate
    x₁  = x₀ + h₀
    y₁  = y₀ + h₀*y′₀ + (h₀²/2)*y″₀
    y′₁ = y′₀ + h₀*y″₀
    y″₁ = ode(x₁, y₁, y′₁)
    y₁  ← y₀ + (h₀/2)*(y′₁ + y′₀) - (h₀²/12)*(y″₁ - y″₀)
    y′₁ ← y′₀ + (h₀/2)*(y″₁ + y″₀)
    y″₁ = ode(x₁, y₁, y′₁)
which is a one-step PECE integrator. Afterwords, refine this estimate for
h according to the formula
    h ← 2|[(∥y₁∥ - ∥y₀∥) / (∥y′₁∥ + ∥y′₀∥)]|
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
        after which a predictor integrates as
            y₁  = y₀ + h*y′₀ + (h²/2)*y″₀
            y′₁ = y′₀ + h*y″₀
        followed by a first approximation for its second-order derivative
            y″₁ = ode(x₁, y₁, y′₁)
        saving
            y″_pred = y″₁
        for use when computing error. A corrector then integrates as
            y₁  ← y₀ + (h/2)*(y′₁ + y′₀) - (h²/12)*(y″₁ - y″₀)
            y′₁ ← y′₀ + (h/2)*(y″₁ + y″₀)
        after which a refined approximation for its derivative is re-evaluated
            y″₁ ← ode(x₁, y₁, y′₁)
        whose local truncation error advances as
            error = (h²/6)*∥y″_pred - y″₀∥
            ε_curr ← 1
            ε_next ← error / max(1, ∥y₁∥)
        Upon completion of a first step, assign
            x_prev  ← x₀     
            y_prev  ← y₀    
            y′_prev ← y′₀
            y″_prev ← y″₀  
            x_curr  ← x₀ + h 
            y_curr  ← y₁ 
            y′_curr ← y′₁
            y″_curr ← y″₁
        which is a second-order accurate one-step method.
    else 
        The main PECE solver begins by advancing the independent variable
            x_next  = x_curr + h
        after which the following predictor integrates the dependent variables
            y_pred  = (1/3)*(4y_curr - y_prev)
                    + (h/6)*(3y′_curr + y′_prev) 
                    + (h²/36)*(31y″_curr - y″_prev)
            y′_pred = (1/3)*(4y′_curr - y′_prev) 
                    + (2h/3)*(2y″_curr - y″_prev) 
            y″_pred = ode(x_next, y_pred, y′_pred)
        while its corrector integrates as
            y_next  ← (1/3)*(4y_curr - y_prev)
                    + (h/24)*(y′_pred + 14y′_curr + y′_prev) 
                    + (h²/72)*(10y″_pred + 51y″_curr - y″_prev) 
            y′_next ← (1/3)*(4y′_curr - y′_prev) 
                    + (2h/3)*y″_pred
            y″_next ← ode(x_next, y_next, y′_next)
        whose local truncation error advances as:
            error = (h/24)*∥y′_pred + 2y′_curr - 3y′_prev
                    + (h/3)*(10y″_pred - 11y″_curr + y″_prev)∥
            ε_curr ← ε_next
            ε_next ← error / max(1, ∥y_next∥)
        This method is second-order accurate, with Gear's BDF2 solving y′.
    end
    
    A PI controller adjusts the local step size h according to the scheme:    
    if ε_next < tol and ε_curr < tol then  
        use the PI controller:
            C = (tol/ε_next)^(0.3/3) * (ε_curr/ε_next)^(0.4/3)
    else
        use the I controller:
            C = (tol/ε_next)^(1/2)
    end
    
    Manage the history by advancing the counters and variables:
    if C > 2 and s > 4 with s mod 2 == 1 then
        x_curr  ← x_next
        y_curr  ← y_next
        y′_curr ← y′_next
        y″_curr ← y″_next
        ε_curr  ← ε_next
        h ← 2h
        s ← (s - 1) ÷ 2
    else if C > 1 then
        x_prev  ← x_curr
        y_prev  ← y_curr
        y′_prev ← y′_curr
        y″_prev ← y″_curr
        x_curr  ← x_next
        y_curr  ← y_next
        y′_curr ← y′_next
        y″_curr ← y″_next
        ε_curr  ← ε_next
        s ← s - 1
    else if ε_next < tol then
        x_prev  ← (1/2)*(x_next + x_curr)
        y_prev  ← (1/2)*(y_next + y_curr) - (h/8)*(y′_next - y′_curr)
        y′_prev ← (1/2)*(y′_next + y′_curr) - (h/8)*(y″_next - y″_curr)
        y″_prev ← (1/8)*(3y″_next + 6y″_curr - y″_prev)
        x_curr  ← x_next 
        y_curr  ← y_next
        y′_curr ← y′_next
        y″_curr ← y″_next
        ε_curr  ← (1/2)*(ε_next + ε_curr)
        h ← h/2
        s ← 2(s - 1)
    else
        x_prev  ← (1/2)*(x_curr + x_prev)
        y_prev  ← (1/2)*(y_curr + y_prev) - (h/8)*(y′_curr - y′_prev)
        y′_prev ← (1/2)*(y′_curr + y′_prev) - (h/8)*(y″_curr - y″_prev)
        y″_prev ← (1/2)*(y″_curr + y″_prev)
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

struct SecondOrderPECE <: PECE
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
    y″_prev::PF.MVector     # Previous value for the derivative d²y/dx².
    y″_curr::PF.MVector     # Current value for the derivative d²y/dx².
    y″_next::PF.MVector     # Next value for the derivative d²y/dx².
    tol::Float64            # Error tolerance targetted by the PI controller.
    ε_curr::PF.MReal        # Truncation error at the current step.
    ε_next::PF.MReal        # Truncation error at the next step.
    steps::PF.MInteger      # Counter for successful steps taken.
    doubled::PF.MInteger    # Counter for times where step size was doubled.
    halved::PF.MInteger     # Counter for times where step size was halved.
    repeats::PF.MInteger    # Counter for times where a step was repeated.
    atNode::PF.MBoolean     # True if local step coincides with a global step.

    # internal constructors
    
    function SecondOrderPECE(my_ode::Function,     # differential equation
                             N::Integer,           # global steps to take
                             x₀::Real,             # solution begins at
                             x_N::Real,            # solution ends at
                             y₀::Vector{<:Real},   # initial condition for y
                             y′₀::Vector{<:Real},  # initial condition for y′
                             tol::Real)            # error tolerance
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
        if length(y′₀) ≠ length(y₀)
            msg = "The length of vectors y and y′ in y″ = f(x,y,y′) differ."
            throw(DimensionMismatch(msg))
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
        if eltype(y′₀) ≠ Float64
            len = length(y′₀)
            y′_0 = Vector{Float64}(undef, len)
            for i in 1:len
                y′_0[i] = convert(Float64, y′₀[i])
            end
            y′₀ = y′_0
        end
        if 3 ≠ (only(methods(my_ode)).nargs - 1)
            error("Function my_ode must have three arguments, viz., x, y and y′.")
        end
        y″₀ = my_ode(x₀, y₀, y′₀)
        if y″₀ isa Vector{<:Real}
            if length(y″₀) ≠ length(y₀)
                msg = "The length of vectors y, y′ and y″ in y″ = f(x,y,y′) differ."
                throw(DimensionMismatch(msg))
            end
            if eltype(y″₀) ≠ Float64
                len = length(y″₀)
                y″_0 = Vector{Float64}(undef, len)
                for i in 1:len
                    y″_0[i] = convert(Float64, y′₀[i])
                end
                y″₀ = y″_0
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
        y₁  = y₀ + h₀*y′₀ + (h₀*h₀/2)*y″₀
        y′₁ = y′₀ + h₀*y″₀
        y″₁ = my_ode(x₁, y₁, y′₁)
        y₁  = y₀ + (h₀/2)*(y′₁ + y′₀) - (h₀*h₀/12)*(y″₁ - y″₀)
        y′₁ = y′₀ + (h₀/2)*(y″₁ + y″₀)
        y″₁ = my_ode(x₁, y₁, y′₁)
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
        y₁  = y₀ + h*y′₀ + (h*h/2)*y″₀
        y′₁ = y′₀ + h*y″₀
        y″₁ = my_ode(x₁, y₁, y′₁)
        err = (h*h/6)*LA.norm(y″₁ - y″₀)
        
        # Finish integration with the corrector.
        y₁  = y₀ + (h/2)*(y′₁ + y′₀) - (h*h/12)*(y″₁ - y″₀)
        y′₁ = y′₀ + (h/2)*(y″₁ + y″₀)
        y″₁ = my_ode(x₁, y₁, y′₁)
        
        # Assign the history variables.
        x_prev  = PF.MReal(x₀)
        x_curr  = PF.MReal(x₀+h)
        x_next  = PF.MReal(0)
        y_prev  = PF.MVector(y₀)
        y_curr  = PF.MVector(y₁)
        y_next  = PF.MVector(length(y₀))
        y′_prev = PF.MVector(y′₀)
        y′_curr = PF.MVector(y′₁)
        y′_next = PF.MVector(length(y₀))
        y″_prev = PF.MVector(y″₀)
        y″_curr = PF.MVector(y″₁)
        y″_next = PF.MVector(length(y₀))
        
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
            y_prev, y_curr, y_next, y′_prev, y′_curr, y′_next, 
            y″_prev, y″_curr, y″_next, tol, ε_curr, ε_next, 
            steps, doubled, halved, repeats, atNode)
    end 
    
    # constructor called by JSON3
    
    function SecondOrderPECE(ode::Function,
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
                             y″_prev::PF.MVector,
                             y″_curr::PF.MVector,
                             y″_next::PF.MVector,
                             tol::Float64,
                             ε_curr::PF.MReal,
                             ε_next::PF.MReal,
                             steps::PF.MInteger,
                             doubled::PF.MInteger,
                             halved::PF.MInteger,
                             repeats::PF.MInteger,
                             atNode::PF.MBoolean)
        
        new(my_ode, N, dx, h, n, s, x₀, y₀, x_prev, x_curr, x_next, 
            y_prev, y_curr, y_next, y′_prev, y′_curr, y′_next, 
            y″_prev, y″_curr, y″_next, tol, ε_curr, ε_next, 
            steps, doubled, halved, repeats, atNode)
    end
end # SecondOrderPECE
 
# Functions that exist for instances of SecondOrderPECE.

StructTypes.StructType(::Type{SecondOrderPECE}) = StructTypes.Struct()

function toFile(y::SecondOrderPECE, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, y)
        write(json_stream, '\n')
    else
        error("The supplied JSON stream is not open.")
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{SecondOrderPECE}, json_stream::IOStream)::SecondOrderPECE
    if isopen(json_stream)
        y = JSON3.read(readline(json_stream), SecondOrderPECE)
    else
        error("The supplied JSON stream is not open.")
    end
    return y
end

function advance!(pece::SecondOrderPECE)
    if pece.n > pece.N
        print("\nThe ODE has been solved.\n")
        return nothing
    end
            
    PF.set!(pece.atNode, false)
    # get history variables
    h       = PF.get(pece.h)
    x_prev  = PF.get(pece.x_prev)
    y_prev  = PF.Vector(pece.y_prev)
    y′_prev = PF.Vector(pece.y′_prev)
    y″_prev = PF.Vector(pece.y″_prev)
    x_curr  = PF.get(pece.x_curr)
    y_curr  = PF.Vector(pece.y_curr)
    y′_curr = PF.Vector(pece.y′_curr)
    y″_curr = PF.Vector(pece.y″_curr)
    ε_curr  = PF.get(pece.ε_curr)
    # advance the independent variable
    x_next = x_curr + h
    # advance the solution with a prediction of
    y_pred  = ((1/3)*(4y_curr - y_prev)
               + (h/6)*(3y′_curr + y′_prev) 
               + (h*h/36)*(31y″_curr - y″_prev))
    y′_pred = ((1/3)*(4y′_curr - y′_prev) 
               + (2h/3)*(2y″_curr - y″_prev))
    # evaluate the ODE
    y″_pred = pece.ode(x_next, y_pred, y′_pred)
    # correct the solution with a corrector of
    y_next  = ((1/3)*(4y_curr - y_prev)
               + (h/24)*(y′_pred + 14y′_curr + y′_prev) 
               + (h*h/72)*(10y″_pred + 51y″_curr - y″_prev))
    y′_next = (1/3)*(4y′_curr - y′_prev) + (2h/3)*y″_pred
    # re-evaluated the ODE
    y″_next = pece.ode(x_next, y_next, y′_next)
    # whose local truncation error advances as
    err = (h/24)*LA.norm(y′_pred + 2y′_curr - 3y′_prev
                    + (h/3)*(10y″_pred - 11y″_curr + y″_prev))
    ε_next = err / max(1, LA.norm(y_next))
    
    # apply the step controller:
    if (ε_curr < pece.tol) && (ε_next < pece.tol) 
        # use a PI controller:
        C = (pece.tol/ε_next)^(0.1) * (ε_curr/ε_next)^(0.4/3)
    else
        # use an I controller:
        C = sqrt(pece.tol/ε_next)
    end
    
    # Manage the history by advancing its counters and variables:
    if (C > 2) && (pece.s > 4) && (pece.s%2 == 1)
        PF.set!(pece.x_curr, x_next)
        PF.set!(pece.ε_curr, ε_next)
        for i in 1:pece.y_curr.len
            pece.y_curr[i]  = y_next[i]
            pece.y′_curr[i] = y′_next[i]
            pece.y″_curr[i] = y″_next[i]
        end
        PF.set!(pece.steps, PF.get(pece.steps)+1)
        PF.set!(pece.doubled, PF.get(pece.doubled)+1)
        PF.set!(pece.h, 2PF.get(pece.h))
        PF.set!(pece.s, (PF.get(pece.s)-1)÷2)
    elseif C > 1
        PF.set!(pece.x_prev, x_curr)
        PF.set!(pece.x_curr, x_next)
        PF.set!(pece.ε_curr, ε_next)
        for i in 1:pece.y_curr.len
            pece.y_prev[i]  = y_curr[i]
            pece.y′_prev[i] = y′_curr[i]
            pece.y″_prev[i] = y″_curr[i]
            pece.y_curr[i]  = y_next[i]
            pece.y′_curr[i] = y′_next[i]
            pece.y″_curr[i] = y″_next[i]
        end
        PF.set!(pece.steps, PF.get(pece.steps)+1)
        PF.set!(pece.s, PF.get(pece.s)-1)
    elseif ε_next < pece.tol    # with C ≤ 1
        PF.set!(pece.x_prev, (1/2)*(x_next+x_curr))
        PF.set!(pece.x_curr, x_next)
        PF.set!(pece.ε_curr, (1/2)*(ε_next+ε_curr))
        for i in 1:pece.y_curr.len
            pece.y_prev[i]  = ((1/2)*(y_next[i] + y_curr[i])
                               - (h/8)*(y′_next[i] - y′_curr[i]))
            pece.y′_prev[i] = ((1/2)*(y′_next[i] + y′_curr[i])
                               - (h/8)*(y″_next[i] - y″_curr[i]))
            pece.y″_prev[i] = (1/8)*(3y″_next[i] + 6y″_curr[i] - y″_prev[i])
            pece.y_curr[i]  = y_next[i]
            pece.y′_curr[i] = y′_next[i]
            pece.y″_curr[i] = y″_next[i]
        end
        PF.set!(pece.steps, PF.get(pece.steps)+1)
        PF.set!(pece.halved, PF.get(pece.halved)+1)
        PF.set!(pece.h, PF.get(pece.h)/2)
        PF.set!(pece.s, 2(PF.get(pece.s)-1))
    else    # ε_next ≥ tol and C ≤ 1
        PF.set!(pece.x_prev, (1/2)*(x_curr+x_prev))
        PF.set!(pece.ε_curr, 1)    # forces the I controller to be used
        for i in 1:pece.y_curr.len
            pece.y_prev[i]  = ((1/2)*(y_curr[i] + y_prev[i])
                               - (h/8)*(y′_curr[i] - y′_prev[i]))
            pece.y′_prev[i] = ((1/2)*(y′_curr[i] + y′_prev[i])
                               - (h/8)*(y″_curr[i] - y″_prev[i]))
            pece.y″_prev[i] = (1/2)*(y″_curr[i] + y″_prev[i])
        end
        PF.set!(pece.repeats, PF.get(pece.repeats)+1)
        PF.set!(pece.h, PF.get(pece.h)/2)
        PF.set!(pece.s, 2PF.get(pece.s))
        advance!(pece)      # repeat this integration step
    end
    if pece.s == 0
        PF.set!(pece.n, PF.get(pece.n)+1)
        PF.set!(pece.s, Int(PF.round(pece.dx/pece.h)))
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
