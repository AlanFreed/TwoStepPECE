To view this file, use a markdown viewer that renders LaTeX commands, e.g., Remarkable.

# Ordinary Differential Equations (ODEs)

This package provides a numerical method for solving systems of first- and second-order ODEs. Specifically, the predict/evaluate/correct/evaluate (PECE) methods of Freed (2017) have been implemented. Because these methods are two-step integrators, they are not self starting, so one-step methods are required to take the first integration step.

For a first-order ODE, one solves
$$
y′  = \mathrm{f}(x, y)
$$
subject to an initial condition of  $y₀$  at  $x₀$  so that
$$
y′₀ = \mathrm{f}(x₀, y₀)
$$
where  $x$  is the independent variable, with  $y$  being the dependent variable, and where  $y$  and  $y′$  are both either scalar, vector or matrix valued.

For a second-order ODE, one solves
$$
y″  = \mathrm{f}(x, y, y′)
$$
subject to initial conditions of  $y₀$  and  $y′₀$  at  $x₀$  so that
$$
y″₀ = \mathrm{f}(x₀, y₀, y′₀)
$$
where  $x$  is the independent variable, with  $y$  being the dependent variable along with  $y′$  being its first-order derivative, or rate, and where  $y$,  $y′$  and  $y″$ are all either scalar, vector or matrix valued.

Notation  $y′$  denotes  $\mathrm{d}y/\mathrm{d}x$, while notation  $y″$  denotes  $\mathrm{d}²y/\mathrm{d}x²$.

# Numerical Methodology

A local solution advances along a sub-grid with a local step size  $h$  that is finer than the global step size  $\mathrm{d}x$  in which  $h$  embeds.  This solution spans an interval  $[x₀, x_N]$  with
$$
\mathrm{d}x = (x_N - x₀) / N
$$ 
and as such there will be  $N$  nodes of integration where solutions are sought. Solution arrays are of length $N+1$ whose first entry associates with the initial condition. The global step $\mathrm{d}x$  is taken to be uniform over the entire span of integration, while the local step  $h$  dynamically adjusts itself to ensure truncation error remains less than an user specified error tolerance denoted as  $tol$.  Global nodes increment as $n = 0, 1, 2, ⋯, N$. Local nodes decrement as $s = S, S-1, S-2, ⋯, 0$, where  $s=S$ and $s=0$  associate with two, neighboring, global nodes, say $n$ and $n+1$.

## Determining an Initial Step Size

To provide an estimate for the initial step size  $h$  to be used when taking the first integration step, we adopt the algorithm of Freed & Iskovitz (1996). From the initial condition, assign
$$
    h₀ = \| y₀ \| / \| y′₀ \|
$$
where  $\| y \|$  is a norm for  $y$.  To help avoid a wind-down or a wind-up instability, constrain this interval so that
$$
\mathrm{d}x/100 < h₀ < \mathrm{d}x/10 .
$$

For a first-order ODE, proceed by integrating
$$
\begin{aligned}
    x₁  & = x₀ + h₀ \\
    y₁  & = y₀ + h₀y′₀ \\
    y′₁ & = \mathrm{f}(x₁, y₁) \\
    y₁  & ← y₀ + (h₀/2)(y′₁ + y′₀) \\
    y′₁ & ← \mathrm{f}(x₁, y₁)
\end{aligned}
$$
which is the predict/evaluate/correct/evaluate (PECE) method of Heun. 

For a second-order ODE, proceed by integrating
$$
\begin{aligned}
    x₁  & = x₀ + h₀ \\
    y₁  & = y₀ + h₀y′₀ + (h₀²/2)y″₀ \\
    y′₁ & = y′₀ + h₀y″₀ \\
    y″₁ & = \mathrm{f}(x₁, y₁, y′₁) \\
    y₁  & ←  y₀ + (h₀/2)(y′₁ + y′₀) - (h₀²/12)(y″₁ - y″₀) \\
    y′₁ & ← y′₀ + (h₀/2)(y″₁ + y″₀) \\
    y″₁ & ← \mathrm{f}(x₁, y₁, y′₁)
\end{aligned}
$$
where the solution for  $y′$  is the same as a solution for $y$ given a first-order ODE, viz., Heun's method.

Afterwords, for both first- and second-order ODEs, refine one's estimate for the initial step size  $h$ according to the formula
$$
h = 2 \left| \frac{\| y₁ \| - \| y₀ \|}{\| y′₁ \| + \| y′₀ \|} \right|
$$
now constrained by 
$$
\mathrm{d}x/1000 < h
$$ 
to help avoid a potential wind-down instability.

The number of local steps required to advance toward the first global node comes from
$$
    S = \mathrm{max}(2, \mathrm{round}(\mathrm{d}x/h))
$$
from which the initial, local, step size is determined to be
$$
h ← \mathrm{d}x / S
$$
where  $\mathrm{d}x$  advances the independent variable $x$  from a current global step to its next global step.

## Truncation Error

A local truncation error comes from taking a difference between corrected and predicted values, viz.,
$$
error = \| y_{corr} - y_{pred} \|
$$
where, given a PECE method, the predicted solution is denoted as $y_{pred}$ while the corrected solution is denoted as $y_{corr}$. 

For first-order ODEs, the local truncation error for the first step of integration is
$$
error = (h/2) \| y′_{pred} - y′₀ \|
$$
while all remaining steps have a local truncation error of
$$
error = (2h/3) \| y′_{pred} - 2y′_{curr} + y′_{prev} \|
$$
with information being kept for the previous, current and next steps of integration, e.g., $y_{prev}$, $y_{curr}$ and $y_{next}$. Note that $\mathrm{d}²f_n = (1/h²)(f_{n+1} - 2f_n + f_{n-1})$ is a finite difference for approximating the second derivative of $f$, where $f = y′$ in this case.

For second-order ODEs, the first solution step has a local truncation error of
$$
error = (h²/6) \| y″_{pred} - y″₀ \|
$$
while all remaining integration steps have a truncation error of
$$ 
error = (h/24) \| y′_{pred} + 2y′_{curr} - 3y′_{prev} + (h/3)(10y″_{pred} - 11y″_{curr} + y″_{prev}) \|
$$
where the coefficients for  $y′$  and  $y″$, when present,  sum to zero for all estimates of truncation error.

The local truncation error  $ε$  is then determined by
$$
ε_{next} = \frac{error}{\mathrm{max}(1, \| y_{next} \|)}
$$
which is an absolute error whenever $\| y_{next} \| < 1$ or a relative error otherwise.

## Two Step PECE Methods

The PECE methods of Freed (2017) are implemented here. Because his integrators are two-step methods, integration must start with a one-step method, with all steps thereafter being integrated with its paired two-step method.

### First-Order ODEs

The one-step method that starts an integration is
$$
\begin{aligned}
    x₁ & = x₀ + h \\
    y₁ & = y₀ + hy′₀ & & \text{P} \\
    y′₁ & = \mathrm{f}(x₁, y₁) & & \text{E} \\
    y₁ & ← y₀ + (h/2)(y′₁ + y′₀) & & \text{C} \\
    y′₁ & ← \mathrm{f}(x₁, y₁) & & \text{E}
\end{aligned}
$$
where the predicted value for  $y′₁$ (the first of two evaluations for  $y′₁$  to appear in the above expressions) is used in the computation of truncation error. Upon completion of a first step, assign
$$
\begin{align}
    x_{prev} & = x₀ \\
    y_{prev} & = y₀ \\
    y′_{prev} & = y′₀ \\
    x_{curr} & = x₀ + h \\
    y_{curr} & = y₁ \\
    y′_{curr} & = y′₁
\end{align}
$$

Integration then continues with a two-step PECE method, spacifically
$$
\begin{aligned}
    x_{next} & = x_{curr} + h & & \\
    y_{next} & = (1/3)(4y_{curr} - y_{prev}) + (2h/3)(2y′_{curr} - y′_{prev}) & & \text{P} \\
    y′_{next} & = \mathrm{f}(x_{next}, y_{next}) & & \text{E} \\
    y_{next} & ← (1/3)(4y_{curr} - y_{prev}) + (2h/3)y′_{next} & & \text{C} \\
    y′_{next} & ← \mathrm{f}(x_{next}, y_{next}) & & \text{E}
\end{aligned}
$$
where the corrector in this PECE method is the well-known BDF2 method (backward difference formula of second order) made popular by Gear.  This method is second-order accurate and, most importantly, A stable. The predicted value for  $y′₁$  is used in the evaluation of truncation error.

### Second-Order ODEs

The one-step method that starts an integration is
$$
\begin{aligned}
    x₁ & = x₀ + h \\
    y₁ & =  y₀ + hy′₀ + (h²/2)y″₀ & & \text{P} \\
    y′₁ & = y′₀ + hy″₀ & & \\
    y″₁ & = \mathrm{f}(x₁, y₁, y′₁) & & \text{E} \\
    y₁ & ←  y₀ + (h/2)(y′₁ + y′₀) - (h²/12)(y″₁ - y″₀) & &  \text{C} \\
    y′₁ & ← y′₀ + (h/2)(y″₁ + y″₀) & & \\
    y″₁ & ← \mathrm{f}(x₁, y₁, y′₁) & & \text{E}
\end{aligned}
$$
where the predicted values for  $y′₁$  and  $y″₁$ (their first appearance in the above formulæ) are used in the evaluation of truncation error. Upon completion of a first step, assign
$$
\begin{align}
    x_{prev} & = x₀ \\
    y_{prev} & = y₀ \\
    y′_{prev} & = y′₀ \\
    y″_{prev} & = y″₀ \\
    x_{curr} & = x₀ + h \\
    y_{curr} & = y₁ \\
    y′_{curr} & = y′₁ \\
    y″_{curr} & = y″₁ 
\end{align}
$$
Integration then continues with a two-step PECE method, specifically
$$
\begin{aligned}
    x_{next} & = x_{curr} + h & & \\
    y_{next} & = (1/3)(4y_{curr} - y_{prev})
                    + (h/6)(3y′_{curr} + y′_{prev}) 
                    + (h²/36)(31y″_{curr} - y″_{prev}) & &  \text{P} \\
    y′_{next} & = (1/3)(4y′_{curr} - y′_{prev}) + (2h/3)(2y″_{curr} - y″_{prev}) & & \\
    y″_{next} & = \mathrm{f}(x_{next}, y_{next}, y′_{next}) & & \text{E} \\
    y_{next} & ← (1/3)*(4y_{curr} - y_{prev}) + (h/24)(y′_{next} + 14y′_{curr} + y′_{prev}) + (h²/72)(10y″_{next} + 51y″_{curr} - y″_{prev}) & & \text{C} \\
    y′_{next} & ← (1/3)(4y′_{curr} - y′_{prev}) + (2h/3)y″_{next} & & \\
    y″_{next} & ← \mathrm{f}(x_{next}, y_{next}, y′_{next}) & & \text{E}
\end{aligned}
$$
where the PECE method for integrating $y′_{next} $ is the same as the PECE method used for solving $y_{next}$ in a first-order ODE. 

Note that the coefficients for $y$ have a weight of 1, while the coefficients for $y′$ have a weight of 2/3, and the coefficients for $y″$, when present, have a weight of 5/6 for both predictor and corrector.

## PI Controller for Managing Step Size

The PI controller of Sőderlind (2002) is used to adjust the local step size  $h$  according to the following strategy.

If $ε_{next} < tol$ and $ε_{curr} < tol$, then use a PI controller, and compute
$$
    C = ( tol / ε_{next} )^{0.3/3} ( ε_{curr} /ε_{next} )^{0.4/3}
$$
otherwise use an I controller, and compute
$$
    C = ( tol / ε_{next})^{1/2}
$$
where $C$ scales the step size $h$ according to the scheme outlined below. The denominator in the exponents is $p+1$ for the PI controller and $p$ for the I controller, wherein $p$ designates the order of the method, which is 2 for the PECE methods presented here.

A conservative strategy is implemented to aid in avoiding a potential wind-down or wind-up instability occurring in the controller; specifically, the step size $h$ will double whenever the truncation error becomes too small, and it will halve whenever the truncation error becomes too large. When on target, the controller will maintain its step size going forward.

### First-Order ODEs

Whenever the scaling factor $C > 2$ and the local step counter $s > 4$ with $s \: \mathrm{mod} \: 2 = 1$, i.e., $s$ is odd, then the ensuing step size will double with the history updating as
$$
\begin{aligned}
    x_{curr}  & ← x_{next} \\
    y_{curr} & ← y_{next} \\
    y′_{curr} & ← y′_{next} \\
    ε_{curr} & ← ε_{next} \\
    h & ← 2h \\
    s & ← (s - 1) ÷ 2
\end{aligned}
$$
Otherwise, if $C > 1$, then the step size is maintained with the history updating as
$$
\begin{aligned}
    x_{prev} & ← x_{curr} \\
    y_{prev} & ← y_{curr} \\
    y′_{prev} & ← y′_{curr} \\
    x_{curr} & ← x_{next} \\
    y_{curr} & ← y_{next} \\
    y′_{curr} & ← y′_{next} \\
    ε_{curr} & ← ε_{next} \\
    s & ← s - 1
\end{aligned}
$$
Otherwise, if $C \le 1$ and $ε_{next} < tol$, then the step size is halved, with previous values being interpolated, and as such the history updates as
$$
\begin{aligned}
    x_{prev} & ← (1/2)(x_{next} + x_{curr}) \\
    y_{prev} & ← (1/2)(y_{next} + y_{curr}) - (h/8)*(y′_{next} - y′_{curr}) \\
    y′_{prev} & ← (1/8)(3y′_{next} + 6y′_{curr} - y′_{prev}) \\
    x_{curr} & ← x_{next} \\
    y_{curr} & ← y_{next} \\
    y′_{curr} & ← y′_{next} \\
    ε_{curr} & ← ε_{next} \\
    h & ← h/2 \\
    s & ← 2(s - 1)
\end{aligned}
$$
Otherwise $C \le 1$ and $ε_{next} \ge tol$; consequently, the step size is halved and this integration step must be repeated, with previous values being interpolated, and as such the history updates as
$$
\begin{aligned}
    x_{prev} & ← (1/2)(x_{curr} + x_{prev}) \\
    y_{prev} & ← (1/2)(y_{curr} + y_{prev}) - (h/8)(y′_{curr} - y′_{prev}) \\
    y′_{prev} & ← (1/2)(y′_{curr} + y′_{prev}) \\
    ε_{curr} & ← 1 \\
    h & ← h/2 \\
    s & ← 2s \\
    \text{Repeat} & \text{ the local integration step.}
\end{aligned}
$$
Whenever $s = 0$, the solution advances to the next global step $n$, specifically
$$
\begin{aligned}
    n & ← n + 1 \\
    s & ← \mathrm{round} ( \mathrm{d}x / h )
\end{aligned}
$$
with integration terminating when $n = N+1$.

### Second-Order ODEs

Whenever the scaling factor $C > 2$ and the local step counter $s > 4$ with $s \: \mathrm{mod} \: 2 = 1$, i.e., $s$ is odd, then the ensuing step size will double with the history updating as
$$
\begin{aligned}
    x_{curr}  & ← x_{next} \\
    y_{curr} & ← y_{next} \\
    y′_{curr} & ← y′_{next} \\
    y″_{curr} & ← y″_{next} \\
    h & ← 2h \\
    s & ← (s - 1) ÷ 2
\end{aligned}
$$
Otherwise, if $C > 1$, then the step size is maintained with the history updating as
$$
\begin{aligned}
    x_{prev} & ← x_{curr} \\
    y_{prev} & ← y_{curr} \\
    y′_{prev} & ← y′_{curr} \\
    y″_{prev} & ← y″_{curr} \\
    x_{curr} & ← x_{next} \\
    y_{curr} & ← y_{next} \\
    y′_{curr} & ← y′_{next} \\
    y″_{curr} & ← y″_{next} \\
    s & ← s - 1
\end{aligned}
$$
Otherwise, if $C \le 1$ and $ε_{next} < tol$, then the step size is halved, with previous values being interpolated, and as such the history updates as
$$
\begin{aligned}
    x_{prev} & ← (1/2)(x_{next} + x_{curr}) \\
    y_{prev} & ← (1/2)(y_{next} + y_{curr}) - (h/8)*(y′_{next} - y′_{curr}) \\
    y′_{prev} & ← (1/2)(y′_{next} + y′_{curr}) - (h/8)(y″_{next} - y″_{curr}) \\
    y″_{prev} & ← (1/8)(3y″_{next} + 6y″_{curr} - y″_{prev}) \\
    x_{curr} & ← x_{next} \\
    y_{curr} & ← y_{next} \\
    y′_{curr} & ← y′_{next} \\
    y″_{curr} & ← y″_{next} \\
    h & ← h/2 \\
    s & ← 2(s - 1)
\end{aligned}
$$
Otherwise $C \le 1$ and $ε_{next} \ge tol$; consequently, the step size is halved and this integration step must be repeated, with previous values being interpolated, and as such the history updates as
$$
\begin{aligned}
    x_{prev} & ← (1/2)(x_{curr} + x_{prev}) \\
    y_{prev} & ← (1/2)(y_{curr} + y_{prev}) - (h/8)(y′_{curr} - y′_{prev}) \\
    y′_{prev} & ← (1/2)(y′_{curr} + y′_{prev}) - (h/8)(y″_{curr} - y″_{prev}) \\
    y″_{prev} & ← (1/2)(y″_{curr} + y″_{prev}) \\
    h & ← h/2 \\
    s & ← 2s \\
    \text{Repeat} & \text{ the local integration step.}
\end{aligned}
$$
Whenever $s = 0$, the solution advances to the next global step $n$, specifically
$$
\begin{aligned}
    n & ← n + 1 \\
    s & ← \mathrm{round} ( \mathrm{d}x / h )
\end{aligned}
$$
with integration terminating when $n = N+1$.

# Software

Software has been written in the julia programming language. To use this module, you will need to add the following repositories to your project:

```julia
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/TwoStepPECEs.jl")
```

All integrator types are concrete implementations of the abstract type
```julia
abstract type PECE end 
```
whose implementations have three common functions, with the correct implementation being selected through multiple dispatch. Two of these functions make the object persistent, i.e.,
```julia
function toFile(pece::PECE, json_stream::IOStream)
```
writes object `pece` to the IOStream `json_stream,` while
```julia
function fromFile(::Type{PECE}, json_stream::IOStream)::PECE
```
reads a PECE object from `json_stream.` Stream `json_stream` can is managed by function calls to
```julia
function PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```
which are exported by the module `PhysicalFields.`

The third common function advances a solution by a single local step, e.g., from step $s$ to step $s-1$, via
```julia
function advance!(pece::PECE)
```
where `pece` is a concrete object that implements a PECE solver.

## First-Order ODEs

Three data structures are available to solve first-order ODEs. One for solving scalar valued ODEs. Another for solving vector valued ODEs. And a third for solving matrix valued ODEs. The PECE methods of Freed (2017) are based upon Gear's, implicit, BDF2 method (backward difference formula of second order). BDF2 is an A stable integrator.

### Scalar Valued ODEs

Given a function ode, i.e., dy/dx, that has an interface of
```julia
    ode = function(x::Real, y::Real)::Real
```
then its associated data structure will be
```julia
struct ScalarGearPECE <: PECE
    ode::Function           # The differential equation that is to be solved.
    N::UInt32               # Number of global steps to be integrated over.
    dx::Float64             # Step size for for the global integrator.
    h::PF.MReal             # Current step size for the local integrator.
    n::PF.MInteger          # Global step counter, increments to N.
    s::PF.MInteger          # Local step counter, decrements to 0.
    x₀::Float64             # Initial value for the independent variable.
    y₀::Float64             # Initial condition for the dependent variable.
    x_prev::PF.MReal        # Previous value for the independent variable.
    x_curr::PF.MReal        # Current value for the independent variable.
    x_next::PF.MReal        # Next value for the independent variable.
    y_prev::PF.MReal        # Previous value for the dependent variable.
    y_curr::PF.MReal        # Current value for the dependent variable.
    y_next::PF.MReal        # Next value for the dependent variable.
    y′_prev::PF.MReal       # Previous value for the derivative dy/dx.
    y′_curr::PF.MReal       # Current value for the derivative dy/dx.
    y′_next::PF.MReal       # Next value for the derivative dy/dx.
    tol::Float64            # Error tolerance targetted by the PI controller.
    ε_curr::PF.MReal        # Truncation error at the current step.
    ε_next::PF.MReal        # Truncation error at the next step.
    steps::PF.MInteger      # Counter for successful steps taken.
    doubled::PF.MInteger    # Counter for times where step size was doubled.
    halved::PF.MInteger     # Counter for times where step size was halved.
    repeats::PF.MInteger    # Counter for times where a step was repeated.
    atNode::PF.MBoolean     # True if the local step coincides with a global step.
end           
```
where `PF` is an alias for `PhysicalFields.` Types `MBoolean,` `MInteger,` `MReal,` `MVector` and `MMatrix` provide mutable types (their values can be changed within an otherwise immutable data structure) for objects that are: a boolean, an integer, a real, a real-valued vector, and a real-valued matrix, respectively.

Its constructor looks like
```julia
function ScalarGearPECE(my_ode::Function,    # The differential equation.
                        N::Integer,          # Global steps to take.
                        x₀::Real,            # Solution begins at.
                        X::Real,             # Solution ends at.
                        y₀::Real,            # Initial condition.
                        tol::Real)           # Error tolerance.
end
```

### Vector Valued ODEs

Given a function ode, i.e., dy/dx, that has an interface of
```julia
    ode = function(x::Real, y::Vector{<:Real})::Vector{<:Real}
```
then its associated data structure will be
```julia
struct VectorGearPECE <: PECE
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
    atNode::PF.MBoolean     # True if the local step coincides with a global step.
end
```
where its constructor looks like
```julia
function VectorGearPECE(my_ode::Function,    # The differential equation.
                        N::Integer,          # Global steps to take.
                        x₀::Real,            # Solution begins at.
                        X::Real,             # Solution ends at.
                        y₀::Vector{<:Real},  # Initial condition.
                        tol::Real)           # Error tolerance.
```

### Matrix Valued ODEs

Given a function ode, i.e., dy/dx, that has an interface of
```julia
    ode = function(x::Real, y::Matrix{<:Real})::Matrix{<:Real}
```
then its associated data structure will be
```julia
struct MatrixGearPECE <: PECE
    ode::Function           # The differential equation that is to be solved.
    N::UInt32               # Number of global steps to be integrated over.
    dx::Float64             # Step size for for the global integrator.
    h::PF.MReal             # Current step size for the local integrator.
    n::PF.MInteger          # Global step counter, increments to N.
    s::PF.MInteger          # Local step counter, decrements to 0.
    x₀::Float64             # Initial value for the independent variable.
    y₀::Matrix{<:Real}      # Initial condition for the dependent variable.
    x_prev::PF.MReal        # Previous value for the independent variable.
    x_curr::PF.MReal        # Current value for the independent variable.
    x_next::PF.MReal        # Next value for the independent variable.
    y_prev::PF.MMatrix      # Previous value for the dependent variable.
    y_curr::PF.MMatrix      # Current value for the dependent variable.
    y_next::PF.MMatrix      # Next value for the dependent variable.
    y′_prev::PF.MMatrix     # Previous value for the derivative dy/dx.
    y′_curr::PF.MMatrix     # Current value for the derivative dy/dx.
    y′_next::PF.MMatrix     # Next value for the derivative dy/dx.
    tol::Float64            # Error tolerance targetted by the PI controller.
    ε_curr::PF.MReal        # Truncation error at the current step.
    ε_next::PF.MReal        # Truncation error at the next step.
    steps::PF.MInteger      # Counter for successful steps taken.
    doubled::PF.MInteger    # Counter for times where step size was doubled.
    halved::PF.MInteger     # Counter for times where step size was halved.
    repeats::PF.MInteger    # Counter for times where a step was repeated.
    atNode::PF.MBoolean     # True if the local step coincides with a global step.
end
```
where its constructor looks like
```julia
function MatrixGearPECE(my_ode::Function,    # The differential equation.
                        N::Integer,          # Global steps to take.
                        x₀::Real,            # Solution begins at.
                        X::Real,             # Solution ends at.
                        y₀::Matrix{<:Real},  # Initial condition.
                        tol::Real)           # Error tolerance.
```

## Second-Order ODEs

# Examples


## The Brusselator

What is known as the Brusselator, which is a chemical kinetics problem, is described by a vector valued ODE whose components are
$$
\begin{aligned}
    \mathrm{d}y_1 / \mathrm{d}t & = A - (B + 1) y_1 + y_1^2 y_2 \\
    \mathrm{d}y_2 / \mathrm{d}t & = B y_1 - y_1^2 y_2
\end{aligned}
$$
Solutions will have a limit cycle whenever, e.g., $A = 1$ and $B = 3$. In contrast, solutions will be stiff whenever, e.g., $A = 1$ and $B = 100$.

For the limit cycle case, let $x₀ = 0$ with initial conditions of $(0.1, 0.1)$, $(1.5, 3)$, $(3, 1)$ and $(3.25, 2.5)$ being considered, running out to $X = 20$.
  
For the stiff case, one can use the same initial conditions, but it is useful to set $X = 0.1$.
    
## An FSAE Race Car

To illustrate a class of problems governed by a second-order ODE that is matrix valued, consider a simple vibration model for a car in three degrees of freedom: bounce, pitch and roll, all measured at the center of gravity of a car.  This example simulates an FSAE race car. (The author of this software derived a variety of vibration models for an FSAE race car in his undergraduate course on vibrations. This is one of those models.)
$$
\begin{aligned}
    x & = \{ b, p, r \}^T   \\
    v & = \{ \mathrm{d}b/\mathrm{d}t, \mathrm{d}p/\mathrm{d}t, \mathrm{d}r/\mathrm{d}t \}^T \\
    a & = \{ \mathrm{d}²b/\mathrm{d}t², \mathrm{d}²p/\mathrm{d}t², \mathrm{d}²r/\mathrm{d}t² \}^T  
\end{aligned}
$$
where $t$ is time (the independent variable) and where $b$ denotes bounce, $p$ denotes pitch, and  $r$ denotes roll (the dependent variables), with $x$ being a displacement vector, $v$ being a velocity vector, and $a$ being an acceleration vector. Bounce is in feet, while pitch and roll are in radians, per FSAE rules. Bounce is positive downward (towards the ground). Pitch is positive whenever the nose is up and the tail is down. Roll is positive whenever the driver side is up and the passenger side is down.

Newton's second law of motion then takes on the form of a matrix equation, viz.,
$$
f(t) = Ma + Cv + Kx
\qquad \therefore \qquad
a = M^{-1} \bigl( f(t) - Cv - Kx \bigr)
$$
where $f$ is a forcing function (a vector that depends upon time), $M$ is a mass matrix, $C$ is a damping matrix, and $K$ is a stiffness matrix.

The mass matrix  $M$ for this 3 degree-of-freedom (DOF) problem is
$$
M = \begin{bmatrix}
    m & 0 & 0 \\
    0 & J_y & 0 \\
    0 & 0 & J_x
\end{bmatrix}
\qquad \text{so that} \qquad
M^{-1} = \begin{bmatrix}
    1/m & 0 & 0 \\
    0 & 1/J_y & 0 \\
    0 & 0 & 1/J_x
\end{bmatrix}
$$
where  $m$  is the mass of the vehicle and driver in slugs, while  $J_x$  and  $J_y$  are the moments of inertia in units of  ft lbs/(rad/sec^2)  about the $x$ and $y$ axes, per FSAE rules. 
    
The symmetric damping matrix  $C$  for this 3 DOF car simulation is
$$
C = \begin{bmatrix}
    C_{11} & C_{12} & C_{13} \\
    C_{12} & C_{22} & C_{23} \\
    C_{13} & C_{23} & C_{33}
\end{bmatrix}
$$
wherein
$$
\begin{aligned}
    C_{11} & = c_1 + c_2 + c_3 + c_4 \\
    C_{12} & = −(c_1 + c_2) l_f + (c_3 + c_4) l_r \\
    C_{13} & = −(c_1 − c_2) r_f + (c_3 − c_4) r_r \\
    C_{22} & = (c_1 + c_2) l_f^2 + (c_3 + c_4) l_r^2 \\
    C_{23} & = -(c_1 − c_2) l_f r_f + (c_3 − c_4) l_r r_r \\
    C_{33} & = (c_1 + c_2) r_f^2 + (c_3 + c_4) r_r^2
\end{aligned}
$$
where $c_1$ is the damping of the driver-front shock absorber, $c_2$ is the damping of the passenger-front shock absorber, $c_3$ is the damping of the passenger-rear shock absorber, $c_4$ is the damping of the driver-rear shock absorber, all of which have units of lb/(ft/sec).  Parameter $l_f$ is the distance from the center of gravity (CG) to the front axle, $l_r$ is the distance from the CG to the rear axle, $r_f$ is the radial distance from the axial centerline (CL) to the center of the tire patch at the front axle, and $r_r$ is the radial distance from the CL to the center of the tire patch at the rear axle, with all distances being measured in feet, per FSAE race rules.
    
The symmetric stiffness matrix  $K$  for this 3 DOF car simulation is
$$
K = \begin{bmatrix}
    K_{11} & K_{12} & K_{13} \\
    K_{12} & K_{22} & K_{23} \\
    K_{13} & K_{23} & K_{33}
\end{bmatrix}
$$
wherein
$$
\begin{aligned}
    K_{11} & = k_1 + k_2 + k_3 + k_4 \\
    K_{12} & = −(k_1 + k_2) l_f + (k_3 + k_4) l_r \\
    K_{13} & = −(k_1 − k_2) r_f + (k_3 − k_4) r_r \\
    K_{22} & = (k_1 + k_2) l_f^2 + (k_3 + k_4) l_r^2 \\
    K_{23} & = -(k_1 − k_2) l_f r_f + (k_3 − k_4) l_r r_r \\
    K_{33} & = (k_1 + k_2) r_f^2 + (k_3 + k_4) r_r^2
\end{aligned}
$$
where $k_1$ is the stiffness of the driver-front spring, $k_2$ is the stiffness of the passenger-front spring, $k_3$ is the stiffness of the passenger-rear spring, $k_4$ is the stiffness of the driver-rear spring, all of which have units of lb/ft, per FSAE rules. The other parameters are as defined for the damping matrix.
        
The forcing function $f(t)$ for thie 3 DOF car simulation is
$$
f = \left\{ \begin{matrix}
    f_1 \\ f_2 \\ f_3
\end{matrix} \right\}
$$
wherein
$$
\begin{aligned}
   f_1 & = w − c_1 Ṙ_1 − c_2 Ṙ_2 − c_3 Ṙ_3 − c_4 Ṙ_4 
          −  k_1 R_1 − k_2 R_2 − k_3 R_3 − k_4 R_4 \\
   f_2 & = ( c_1 Ṙ_1 + c_2 Ṙ_2 + k_1 R_1 + k_2 R_2 ) l_f 
          − (c_3 Ṙ_3 + c_4 Ṙ_4 + k_3 R_3 + k_4 R_4) l_r \\
   f_3 & = (c_1 Ṙ_1 − c_2 Ṙ_2 + k_1 R_1 − k_2 R_2) r_f 
          − (c_3 Ṙ_3 − c_4 Ṙ_4 + k_3 R_3 − k_4 R_4) r_r
\end{aligned}
$$
where  $w$  is the weight of the car and driver in pounds,  $R_1$, $R_2$, $R_3$, $R_4$ are the upward displacements of the roadway, which are functions of time, and $Ṙ_1$, $Ṙ_2$, $Ṙ_3$, $Ṙ_4$  are their time rates of change. Units are in ft and ft/sec, respectively.  The other parameters are as defined for the damping and stiffness matrices.
    
Representative parameters for a typical FSAE race car with driver are:
$$
\begin{aligned}
    m & = 14       & & \text{slug} \\
    w & = 450      & & \text{lbf} \\
    J_x & = 20     & & \text{ft}\cdot\text{lbf/(rad/sec^2)} \\
    J_y & = 45     & & \text{ft}\cdot\text{lbf/(rad/sec^2)} \\
    l_f & = 3.2    & & \text{ft} \\
    l_r & = 1.8    & & \text{ft} \\
    r_f & = 2.1    & & \text{ft} \\
    r_r & = 2      & & \text{ft} \\
    c_1 & = 120    & & \text{lbf/(ft/sec)} \\
    c_2 & = 120    & & \text{lbf/(ft/sec)} \\
    c_3 & = 180    & & \text{lbf/(ft/sec)} \\
    c_4 & = 180    & & \text{lbf/(ft/sec)} \\
    k_1 & = 1800   & & \text{lbf/ft} \\
    k_2 & = 1800   & & \text{lbf/ft} \\
    k_3 & = 3600   & & \text{lbf/ft} \\
    k_4 & = 3600   & & \text{lbf/ft}
\end{aligned}
$$
where units are per FSAE race rules, with $w = mg$ where $g = 32.2 \; \text{ft/sec}^2$ is the acceleration due to gravity, and where a $\text{slug} = \text{1 lbf}\cdot\text{sec}^2\text{/ft}$ and a moment of inertia $J = \text{1 slug}\cdot\text{ft}^2$ or $J = \text{1 slug}\cdot\text{ft}^2\text{/rad} = \text{1 ft}\cdot\text{lbf/(rad/sec}^2\text{)}$ in polar coordinates.

### Roadway Forcing Function

The forcing function $f$ that excites a vehicle's motion comes from the roadway being driven; specifically, four functions are considered
$$
\begin{aligned}
    (R_1, Ṙ_1) & = \mathrm{roadwayDF}(t) \\
    (R_2, Ṙ_2) & = \mathrm{roadwayPF}(t) \\
    (R_3, Ṙ_3) & = \mathrm{roadwayPR}(t) \\
    (R_4, Ṙ_4) & = \mathrm{roadwayDR}(t)
\end{aligned}
$$
each of which returns a tuple comprised of the roadway's vertical position and vertical velocity at its respective tire, with `DR` implying drive front, `PF` implying passenger front, `PR` implying passenger rear, and `DR` implying driver rear.

A familiar roadway roadway obstacle is the speed bump, which is modeled as

```julia
struct bump
    height::Real    # height of a bump
    width::Real     # width of a bump
    top::Real       # length of flat region atop a bump
end
```

def bump(height, length, top, x, v):
    # geometric properties of the bump: dimensions are in ft
    height = 1. / 6.   # height of the bump
    length = 2.        # length of the bump
    top = 0.5          # length of flat region on top of bump
    
    # a haversine bump
    if (x <= 0.) or (x >= length):
        # located either ahead of or behind the bump
        R = 0.
        dR = 0.
    else:
        # locate where your position is on the bump
        if x < (length - top) / 2.:
            phi = 2. * math.pi * x / (length - top)
        elif x < (length + top) / 2.:
            phi = math.pi
        else:
            phi = 2 * math.pi * (x - top) / (length - top)
        R = (1. - math.cos(phi)) * height / 2.
        dR = (math.pi * height / (length - top)) * math.sin(phi) * v
        
    return R, dR


def trajectory(t):
    # consider constant speed
    mph2fps = 1.467
    speed = 10 * mph2fps
    
    x = speed * t
    v = speed
    return x, v
    
def mogul(position, speed):
    # properties of the mogul field
    height = 0.25     # height of each bump in the moguls
    length = 5.       # wavelength of each bump
    top = 2.          # length of flat reagion on top of each bump

    if position >= 0. and position < length:
        location = position
    elif position >= length and position < 2. * length:
        location = position - length
    elif position >= 2. * length and position < 3. * length:
        location = position - 2. * length
    elif position >= 3. * length and position < 4. * length:
        location = position - 3. * length
    elif position >= 4. * length and position < 5. * length:
        location = position - 4. * length
    else:
       location = -1.
    R, dR = bump(height, length, top, location, speed)
    return R, dR
    

def roadwayDF(t):
    position, speed = trajectory(t)
    return mogul(position, speed)
    
    
def roadwayPF(t):
    offset = 0.5         # distance passenger side trails the driver's side
    
    position, speed = trajectory(t)
    position = position - offset
    return mogul(position, speed)
    
    
def roadwayPR(t):
    offset = 0.5         # distance passenger side trails the driver's side
    wheelbase = 6.       # distance rear axle is behind the front axle
    
    position, speed = trajectory(t)
    position = position - wheelbase - offset
    return mogul(position, speed)
    
    
def roadwayDR(t):
    wheelbase = 6.       # distance rear axle is behind the front axle
    
    position, speed = trajectory(t)
    position = position - wheelbase
    return mogul(position, speed)
    

# References

1. G. Sőderlind, "Automatic control and adaptive time-stepping", *Numerical Algorithms*, **31**, 2002, 281-310.

2. A. D. Freed, "A Technical Note: Two-Step PECE Methods for Approximating Solutions to First- and Second-Order ODEs", arXiv:1707.02125 [cs.{NA}], 2017.

3. A. D. Freed and I. Iskovitz, "Development and Applications of a Rosenbrock Integrator," NASA TM 4709, 1996.
