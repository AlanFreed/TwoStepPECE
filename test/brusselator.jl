"""
-------------------------------------------------------------------------------

    author:  Alan Freed
    date:    March 20, 2026
    updated: Jun 16, 2026

    This problem is known as the Brusselator:
    
        dy1/dt = A - (B + 1) y1 + y1^2 y2
        dy2/dt = B y1 - y1^2 y2
        
    Solutions will have a limit cycle when A = 1 and B = 3.  For x₀ = 0, use
    initial conditions of (0.1, 0.1), (1.5, 3), (3, 1), (3.25, 2.5) running 
    out to X = 20.

    Solutions will be stiff when A = 1 and B = 100.  Here one can use the same
    initial conditions, but it is useful to set X = 0.1.
    
    Reference:
    A. D. Freed and I. Iskovitz, "Development and Applications of a Rosenbrock
    Integrator," NASA TM 4709, 1996.
"""
module brusselator

using
    Measures,        # needed to pad white space (margins) around a plot
    PhysicalFields,
    Plots,
    TwoStepPECE
    
import
    PhysicalFields as PF,
    TwoStepPECE    as PECE
    
export
    run
    
struct Brusselator <: PECE.Model
    # the model
    #   y₁′ = A - (B+1)y₁ + y₁²y₂
    #   y₂′ = By₁ - y₁²y₂
    # its parameters
    A::Real   # model parameter A
    B::Real   # model parameter B
    
    function Brusselator(A::Real, B::Real)
        new(A, B)
    end
end
       
# Functions that exist for instances of Brusselator.

# The ODE governing control.
function fₓ(model::Brusselator, t::Real, x::Vector{<:Real}, y::Vector{<:Real})::Vector{<:Real}
    # arguments
    #    model  the model's constants wrt control.
    #    t      the independent variable, typically time.
    #    x      the control variables.
    #    y      the response variables.
    # returned
    #    x′     ODE that describes the rate-of-change in control wrt t.
    
    # There is no control function x for the brusselator.
    x′ = Vector{Float64}(undef, 0)
    
    return x′
end # fₓ

# The ODE governing response.
function fᵧ(model::Brusselator, t::Real, x::Vector{<:Real}, y::Vector{<:Real})::Vector{<:Real}
    # arguments
    #    model  the model's constants wrt response.
    #    t      the independent variable, usually time.
    #    x      the control variables.
    #    y      the response variables.
    # returned
    #    y′     ODE that describes the rate-of-change in response wrt t.
    
    # The response function for the Brusselator.
    A = model.A
    B = model.B
    y′ = Vector{Float64}(undef, 2)
    y′[1] = A - (B + 1)*y[1] + y[1]*y[1]*y[2]
    y′[2] = B*y[1] - y[1]*y[1]*y[2]
    
    return y′
end # fᵧ

function run()
    N   = 250
    t₀  = 0.0
    T   = 20.0
    tol = 0.0001
    
    # Model parameters for the response: [A, B]
    limitCycle = Brusselator(1.0, 3.0)
    stiffModel = Brusselator(100.0, 3.0)
    
    # Initial conditions:
    # control
    x₀ = Vector{Float64}(undef, 0)
    # repsonse (Steady state is at [1, 3], so y₀ cannot be [1, 3]!)
    y₀_1 = [0.0, 1.0]
    y₀_2 = [1.0, 4.0]
    y₀_3 = [3.0, 0.0]
    y₀_4 = [4.0, 2.0]
    
    # The first initial condition for the cyclic solution.
    solver1 = FirstOrderPECE(fₓ, fᵧ, limitCycle, N, t₀, T, x₀, y₀_1, tol)
    t1 = zeros(Float64, N+1)
    t1[1] = t₀
    for i = 2:N+1
        t1[i] = t1[i-1] + solver1.dt
    end
    err1 = zeros(Float64, N+1)
    y1_1 = zeros(Float64, N+1)
    y1_2 = zeros(Float64, N+1)
    y1_1[1] = y₀_1[1]
    y1_2[1] = y₀_1[2]
    i = 1
    while solver1.n < solver1.N
        advance!(solver1)
        if solver1.atNode == true
            i = i + 1
            err1[i] = PF.get(solver1.ε_curr)
            y1_1[i] = solver1.y_curr[1]
            y1_2[i] = solver1.y_curr[2]
        end
    end
    err1[1] = err1[2]
    print("\nThe Brusselator with ICs (0, 1) ran with statistics:\n")
    print("   ", PF.toString(solver1.steps), " steps taken with ",
          PF.toString(solver1.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver1.doubled), " were doubled and ",
          PF.toString(solver1.halved), " were halved.\n")
    
    # The second initial condition for the cyclic solution.
    solver2 = FirstOrderPECE(fₓ, fᵧ, limitCycle, N, t₀, T, x₀, y₀_2, tol)
    t2 = zeros(Float64, N+1)
    t2[1] = t₀
    for i = 2:N+1
        t2[i] = t2[i-1] + solver2.dt
    end
    err2 = zeros(Float64, N+1)
    y2_1 = zeros(Float64, N+1)
    y2_2 = zeros(Float64, N+1)
    y2_1[1] = y₀_2[1]
    y2_2[1] = y₀_2[2]
    i = 1
    while solver2.n < solver2.N
        advance!(solver2)
        if solver2.atNode == true
            i = i + 1
            err2[i] = PF.get(solver2.ε_curr)
            y2_1[i] = solver2.y_curr[1]
            y2_2[i] = solver2.y_curr[2]
        end
    end
    err2[1] = err2[2]
    print("\nThe Brusselator with ICs (1, 4) ran with statistics:\n")
    print("   ", PF.toString(solver2.steps), " steps taken with ",
          PF.toString(solver2.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver2.doubled), " were doubled and ",
          PF.toString(solver2.halved), " were halved.\n")
    
    # The third initial condition for the cyclic solution.
    solver3 = FirstOrderPECE(fₓ, fᵧ, limitCycle, N, t₀, T, x₀, y₀_3, tol)
    t3 = zeros(Float64, N+1)
    t3[1] = t₀
    for i = 2:N+1
        t3[i] = t3[i-1] + solver3.dt
    end
    err3 = zeros(Float64, N+1)
    y3_1 = zeros(Float64, N+1)
    y3_2 = zeros(Float64, N+1)
    y3_1[1] = y₀_3[1]
    y3_2[1] = y₀_3[2]
    i = 1
    while solver3.n < solver3.N
        advance!(solver3)
        if solver3.atNode == true
            i = i + 1
            err3[i] = PF.get(solver3.ε_curr)
            y3_1[i] = solver3.y_curr[1]
            y3_2[i] = solver3.y_curr[2]
        end
    end
    err3[1] = err3[2]
    print("\nThe Brusselator with ICs (3, 0) ran with statistics:\n")
    print("   ", PF.toString(solver3.steps), " steps taken with ",
          PF.toString(solver3.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver3.doubled), " were doubled and ",
          PF.toString(solver3.halved), " were halved.\n")
    
    # The fourth initial condition for the cyclic solution.
    solver4 = FirstOrderPECE(fₓ, fᵧ, limitCycle, N, t₀, T, x₀, y₀_4, tol)
    t4 = zeros(Float64, N+1)
    t4[1] = t₀
    for i = 2:N+1
        t4[i] = t4[i-1] + solver4.dt
    end
    err4 = zeros(Float64, N+1)
    y4_1 = zeros(Float64, N+1)
    y4_2 = zeros(Float64, N+1)
    y4_1[1] = y₀_4[1]
    y4_2[1] = y₀_4[2]
    i = 1
    while solver4.n < solver4.N
        advance!(solver4)
        if solver4.atNode == true
            i = i + 1
            err4[i] = PF.get(solver4.ε_curr)
            y4_1[i] = solver4.y_curr[1]
            y4_2[i] = solver4.y_curr[2]
        end
    end
    err4[1] = err4[2]
    print("\nThe Brusselator with ICs (4, 2) ran with statistics:\n")
    print("   ", PF.toString(solver4.steps), " steps taken with ",
          PF.toString(solver4.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver4.doubled), " were doubled and ",
          PF.toString(solver4.halved), " were halved.\n")
    
    # for a stiff system of ODEs
    T  = 0.5
    
    # initial conditions:
    # control
    x₀ = Vector{Float64}(undef, 0)
    # repsonse (Steady state is at [1, 3], so y₀ can't be [1, 3]!)
    y₀_1 = [0.0, 1.0]
    y₀_2 = [1.0, 4.0]
    y₀_3 = [3.0, 0.0]
    y₀_4 = [4.0, 2.0]
    
    # The first initial condition for the stiff solution.
    solver5 = FirstOrderPECE(fₓ, fᵧ, stiffModel, N, t₀, T, x₀, y₀_1, tol)
    t5 = zeros(Float64, N+1)
    t5[1] = t₀
    for i = 2:N+1
        t5[i] = t5[i-1] + solver5.dt
    end
    err5 = zeros(Float64, N+1)
    y5_1 = zeros(Float64, N+1)
    y5_2 = zeros(Float64, N+1)
    y5_1[1] = y₀_1[1]
    y5_2[1] = y₀_1[2]
    i = 1
    while solver5.n < solver5.N
        advance!(solver5)
        if solver5.atNode == true
            i = i + 1
            err5[i] = PF.get(solver5.ε_curr)
            y5_1[i] = solver5.y_curr[1]
            y5_2[i] = solver5.y_curr[2]
        end
    end
    err5[1] = err5[2]
    print("\nThe stiff Brusselator with ICs (0, 1) ran with statistics:\n")
    print("   ", PF.toString(solver5.steps), " steps taken with ",
          PF.toString(solver5.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver5.doubled), " were doubled and ",
          PF.toString(solver5.halved), " were halved.\n")
    
    # The second initial condition for the stiff solution.
    solver6 = FirstOrderPECE(fₓ, fᵧ, stiffModel, N, t₀, T, x₀, y₀_2, tol)
    t6 = zeros(Float64, N+1)
    t6[1] = t₀
    for i = 2:N+1
        t6[i] = t6[i-1] + solver6.dt
    end
    err6 = zeros(Float64, N+1)
    y6_1 = zeros(Float64, N+1)
    y6_2 = zeros(Float64, N+1)
    y6_1[1] = y₀_2[1]
    y6_2[1] = y₀_2[2]
    i = 1
    while solver6.n < solver6.N
        advance!(solver6)
        if solver6.atNode == true
            i = i + 1
            err6[i] = PF.get(solver6.ε_curr)
            y6_1[i] = solver6.y_curr[1]
            y6_2[i] = solver6.y_curr[2]
        end
    end
    err6[1] = err6[2]
    print("\nThe stiff Brusselator with ICs (1, 4) ran with statistics:\n")
    print("   ", PF.toString(solver6.steps), " steps taken with ",
          PF.toString(solver6.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver6.doubled), " were doubled and ",
          PF.toString(solver6.halved), " were halved.\n")
    
    # The third initial condition for the stiff solution.
    solver7 = FirstOrderPECE(fₓ, fᵧ, stiffModel, N, t₀, T, x₀, y₀_3, tol)
    t7 = zeros(Float64, N+1)
    t7[1] = t₀
    for i = 2:N+1
        t7[i] = t7[i-1] + solver7.dt
    end
    err7 = zeros(Float64, N+1)
    y7_1 = zeros(Float64, N+1)
    y7_2 = zeros(Float64, N+1)
    y7_1[1] = y₀_3[1]
    y7_2[1] = y₀_3[2]
    i = 1
    while solver7.n < solver7.N
        advance!(solver7)
        if solver7.atNode == true
            i = i + 1
            err7[i] = PF.get(solver7.ε_curr)
            y7_1[i] = solver7.y_curr[1]
            y7_2[i] = solver7.y_curr[2]
        end
    end
    err7[1] = err7[2]
    print("\nThe stiff Brusselator with ICs (3, 0) ran with statistics:\n")
    print("   ", PF.toString(solver7.steps), " steps taken with ",
          PF.toString(solver7.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver7.doubled), " were doubled and ",
          PF.toString(solver7.halved), " were halved.\n")
    
    # The fourth initial condition for the stiff solution.
    solver8 = FirstOrderPECE(fₓ, fᵧ, stiffModel, N, t₀, T, x₀, y₀_4, tol)
    t8 = zeros(Float64, N+1)
    t8[1] = t₀
    for i = 2:N+1
        t8[i] = t8[i-1] + solver8.dt
    end
    err8 = zeros(Float64, N+1)
    y8_1 = zeros(Float64, N+1)
    y8_2 = zeros(Float64, N+1)
    y8_1[1] = y₀_4[1]
    y8_2[1] = y₀_4[2]
    i = 1
    while solver8.n < solver8.N
        advance!(solver8)
        if solver8.atNode == true
            i = i + 1
            err8[i] = PF.get(solver8.ε_curr)
            y8_1[i] = solver8.y_curr[1]
            y8_2[i] = solver8.y_curr[2]
        end
    end
    err8[1] = err8[2]
    print("\nThe stiff Brusselator with ICs (4, 2) ran with statistics:\n")
    print("   ", PF.toString(solver8.steps), " steps taken with ",
          PF.toString(solver8.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver8.doubled), " were doubled and ",
          PF.toString(solver8.halved), " were halved.\n")
          
    # set the graphics backend to GR
    ENV["QT_QPA_PLATFORM"] = "wayland"
    gr()
    
    plot(y1_1,  y1_2, label="y₀ = [0, 1]", linecolor=:black, linewidth=3)
    plot!(y2_1, y2_2, label="y₀ = [1, 4]", linecolor=:blue,  linewidth=3)
    plot!(y3_1, y3_2, label="y₀ = [3, 0]", linecolor=:red,   linewidth=3)
    plot!(y4_1, y4_2, label="y₀ = [4, 2]", linecolor=:cyan,  linewidth=3)
    plot!(size=(500, 500))
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!(string("Brusselator with Limit Cycle: Values"))
    xlabel!("y₁")
    ylabel!("y₂")
    dirPath = string(pwd(), "/figures/")
    if !isdir(dirPath)
        mkdir(dirPath)
    end
    figName = string("limitCycleBrusselator.png")
    figPath = string(dirPath, figName)
    savefig(figPath)

    plot(t1, err1,  label="y₀ = [0, 1]", linecolor=:black, linewidth=3)
    plot!(t2, err2, label="y₀ = [1, 4]", linecolor=:blue,  linewidth=3)
    plot!(t3, err3, label="y₀ = [3, 0]", linecolor=:red,   linewidth=3)
    plot!(t4, err4, label="y₀ = [4, 2]", linecolor=:cyan,  linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(yscale=:log10, minorgrid=true)
    ylims!(1e-5*tol, 10*tol)
    plot!(legend=:bottom)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!("Brusselator with Limit Cycle: Errors")
    xlabel!("t")
    ylabel!("Truncation Error, ε")
    figName = string("limitCycleBrusselatorError.png")
    figPath = string(dirPath, figName)
    savefig(figPath)

    plot(y5_1,  y5_2, label="y₀ = [0, 1]", linecolor=:black, linewidth=3)
    plot!(y6_1, y6_2, label="y₀ = [1, 4]", linecolor=:blue,  linewidth=3)
    plot!(y7_1, y7_2, label="y₀ = [3, 0]", linecolor=:red,   linewidth=3)
    plot!(y8_1, y8_2, label="y₀ = [4, 2]", linecolor=:cyan,  linewidth=3)
    plot!(size=(500, 500))
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!(string("Stiff Brusselator"))
    xlabel!("y₁")
    ylabel!("y₂")
    figName = string("stiffBrusselatorY1vsY2.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(t5, err5,  label="y₀ = [0, 1]", linecolor=:black, linewidth=3)
    plot!(t6, err6, label="y₀ = [1, 4]", linecolor=:blue,  linewidth=3)
    plot!(t7, err7, label="y₀ = [3, 0]", linecolor=:red,   linewidth=3)
    plot!(t8, err8, label="y₀ = [4, 2]", linecolor=:cyan,  linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(yscale=:log10, minorgrid=true)
    ylims!(1e-5*tol, 10*tol)
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!("Stiff Brusselator: Errors")
    xlabel!("t")
    ylabel!("Truncation Error, ε")
    figName = string("stiffBrusselatorError.png")
    figPath = string(dirPath, figName)
    savefig(figPath)

    return nothing
end # run

end # brusselator