"""
    author:  Alan Freed
    date:    March 20, 2026

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
    PhysicalFields as PF
    
export
    run
    
function limitCycle(x::Real, y::Vector{Float64})::Vector{Float64}
    A = 1
    B = 3
    ode = Vector{Float64}(undef, 2)
    ode[1] = A - (B + 1)*y[1] + y[1]*y[1]*y[2]
    ode[2] = B*y[1] - y[1]*y[1]*y[2]
    return ode
end # limitCycle

function stiff(x::Real, y::Vector{Float64})::Vector{Float64}
    A = 1
    B = 100
    ode = Vector{Float64}(undef, 2)
    ode[1] = A - (B + 1)*y[1] + y[1]*y[1]*y[2]
    ode[2] = B*y[1] - y[1]*y[1]*y[2]
    return ode
end # stiff

function run()
    N   = 250
    x₀  = 0.0
    X   = 20.0
    tol = 0.0001
    
    # the first initial condition
    y₀ = [0.1, 1.0]
    solver1 = FirstOrderPECE(limitCycle, N, x₀, X, y₀, tol)
    x1 = zeros(Float64, N+1)
    x1[1] = x₀
    for i = 2:N+1
        x1[i] = x1[i-1] + solver1.dx
    end
    err1 = zeros(Float64, N+1)
    y1_1 = zeros(Float64, N+1)
    y1_2 = zeros(Float64, N+1)
    err1[1] = tol
    y1_1[1] = y₀[1]
    y1_2[1] = y₀[2]
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
    print("\nThe Brusselator with ICs (0.1, 1.0) ran with statistics:\n")
    print("   ", PF.toString(solver1.steps), " steps taken with ",
          PF.toString(solver1.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver1.doubled), " were doubled and ",
          PF.toString(solver1.halved), " were halved.\n")
    
    # the second initial condition
    y₀ = [1.5, 3.0]
    solver2 = FirstOrderPECE(limitCycle, N, x₀, X, y₀, tol)
    x2 = zeros(Float64, N+1)
    x2[1] = x₀
    for i = 2:N+1
        x2[i] = x2[i-1] + solver2.dx
    end
    err2 = zeros(Float64, N+1)
    y2_1 = zeros(Float64, N+1)
    y2_2 = zeros(Float64, N+1)
    err2[1] = tol
    y2_1[1] = y₀[1]
    y2_2[1] = y₀[2]
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
    print("\nThe Brusselator with ICs (1.5, 3.0) ran with statistics:\n")
    print("   ", PF.toString(solver2.steps), " steps taken with ",
          PF.toString(solver2.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver2.doubled), " were doubled and ",
          PF.toString(solver2.halved), " were halved.\n")
    
    # the third initial condition
    y₀ = [3.0, 1.0]
    solver3 = FirstOrderPECE(limitCycle, N, x₀, X, y₀, tol)
    x3 = zeros(Float64, N+1)
    x3[1] = x₀
    for i = 2:N+1
        x3[i] = x3[i-1] + solver3.dx
    end
    err3 = zeros(Float64, N+1)
    y3_1 = zeros(Float64, N+1)
    y3_2 = zeros(Float64, N+1)
    err3[1] = tol
    y3_1[1] = y₀[1]
    y3_2[1] = y₀[2]
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
    print("\nThe Brusselator with ICs (3.0, 1.0) ran with statistics:\n")
    print("   ", PF.toString(solver3.steps), " steps taken with ",
          PF.toString(solver3.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver3.doubled), " were doubled and ",
          PF.toString(solver3.halved), " were halved.\n")
    
    # the fourth initial condition
    y₀ = [3.5, 2.5]
    solver4 = FirstOrderPECE(limitCycle, N, x₀, X, y₀, tol)
    x4 = zeros(Float64, N+1)
    x4[1] = x₀
    for i = 2:N+1
        x4[i] = x4[i-1] + solver4.dx
    end
    err4 = zeros(Float64, N+1)
    y4_1 = zeros(Float64, N+1)
    y4_2 = zeros(Float64, N+1)
    err4[1] = tol
    y4_1[1] = y₀[1]
    y4_2[1] = y₀[2]
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
    print("\nThe Brusselator with ICs (3.5, 2.5) ran with statistics:\n")
    print("   ", PF.toString(solver4.steps), " steps taken with ",
          PF.toString(solver4.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver4.doubled), " were doubled and ",
          PF.toString(solver4.halved), " were halved.\n")
    
    # for a stiff system of ODEs
    # the first initial condition
    X  = 0.1
    y₀ = [0.1, 1.0]
    solver5 = FirstOrderPECE(stiff, N, x₀, X, y₀, tol)
    x5 = zeros(Float64, N+1)
    x5[1] = x₀
    for i = 2:N+1
        x5[i] = x5[i-1] + solver5.dx
    end
    err5 = zeros(Float64, N+1)
    y5_1 = zeros(Float64, N+1)
    y5_2 = zeros(Float64, N+1)
    x5[1] = x₀
    for i = 2:N+1
        x5[i] = x5[i-1] + solver5.dx
    end
    err5[1] = tol
    y5_1[1] = y₀[1]
    y5_2[1] = y₀[2]
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
    print("\nThe stiff Brusselator with ICs (0.1, 1.0) ran with statistics:\n")
    print("   ", PF.toString(solver5.steps), " steps taken with ",
          PF.toString(solver5.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver5.doubled), " were doubled and ",
          PF.toString(solver5.halved), " were halved.\n")
    
    # the second initial condition
    y₀ = [1.5, 3.0]
    solver6 = FirstOrderPECE(stiff, N, x₀, X, y₀, tol)
    x6 = zeros(Float64, N+1)
    x6[1] = x₀
    for i = 2:N+1
        x6[i] = x6[i-1] + solver6.dx
    end
    err6 = zeros(Float64, N+1)
    y6_1 = zeros(Float64, N+1)
    y6_2 = zeros(Float64, N+1)
    err6[1] = tol
    y6_1[1] = y₀[1]
    y6_2[1] = y₀[2]
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
    print("\nThe stiff Brusselator with ICs (1.5, 3.0) ran with statistics:\n")
    print("   ", PF.toString(solver6.steps), " steps taken with ",
          PF.toString(solver6.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver6.doubled), " were doubled and ",
          PF.toString(solver6.halved), " were halved.\n")
    
    # the third initial condition
    y₀ = [3.0, 1.0]
    solver7 = FirstOrderPECE(stiff, N, x₀, X, y₀, tol)
    x7 = zeros(Float64, N+1)
    x7[1] = x₀
    for i = 2:N+1
        x7[i] = x7[i-1] + solver7.dx
    end
    err7 = zeros(Float64, N+1)
    y7_1 = zeros(Float64, N+1)
    y7_2 = zeros(Float64, N+1)
    err7[1] = tol
    y7_1[1] = y₀[1]
    y7_2[1] = y₀[2]
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
    print("\nThe stiff Brusselator with ICs (3.0, 1.0) ran with statistics:\n")
    print("   ", PF.toString(solver7.steps), " steps taken with ",
          PF.toString(solver7.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver7.doubled), " were doubled and ",
          PF.toString(solver7.halved), " were halved.\n")
    
    # the fourth initial condition
    y₀ = [3.5, 2.5]
    solver8 = FirstOrderPECE(stiff, N, x₀, X, y₀, tol)
    x8 = zeros(Float64, N+1)
    x8[1] = x₀
    for i = 2:N+1
        x8[i] = x8[i-1] + solver8.dx
    end
    err8 = zeros(Float64, N+1)
    y8_1 = zeros(Float64, N+1)
    y8_2 = zeros(Float64, N+1)
    err8[1] = tol
    y8_1[1] = y₀[1]
    y8_2[1] = y₀[2]
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
    print("\nThe stiff Brusselator with ICs (3.5, 2.5) ran with statistics:\n")
    print("   ", PF.toString(solver8.steps), " steps taken with ",
          PF.toString(solver8.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver8.doubled), " were doubled and ",
          PF.toString(solver8.halved), " were halved.\n")
          
    # set the graphics backend to GR
    ENV["QT_QPA_PLATFORM"] = "wayland"
    gr()
    
    plot(y1_1,  y1_2, label="y₀ = [0.1, 0.1]", linecolor=:black, linewidth=3)
    plot!(y2_1, y2_2, label="y₀ = [1.5, 3.0]", linecolor=:blue,  linewidth=3)
    plot!(y3_1, y3_2, label="y₀ = [3.0, 1.0]", linecolor=:red,   linewidth=3)
    plot!(y4_1, y4_2, label="y₀ = [3.5, 2.5]", linecolor=:cyan,  linewidth=3)
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

    plot(x1, err1,  label="y₀ = [0.1, 0.1]", linecolor=:black, linewidth=3)
    plot!(x2, err2, label="y₀ = [1.5, 3.0]", linecolor=:blue,  linewidth=3)
    plot!(x3, err3, label="y₀ = [3.0, 1.0]", linecolor=:red,   linewidth=3)
    plot!(x4, err4, label="y₀ = [3.5, 2.5]", linecolor=:cyan,  linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(yscale=:log10, minorgrid=true)
    ylims!(1e-5*tol, 10*tol)
    plot!(legend=:bottom)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!("Brusselator with Limit Cycle: Errors")
    xlabel!("x")
    ylabel!("Truncation Error, ε")
    figName = string("limitCycleBrusselatorError.png")
    figPath = string(dirPath, figName)
    savefig(figPath)

    plot(x5, y5_1,  label="y₀ = [0.1, 0.1]", linecolor=:black, linewidth=3)
    plot!(x6, y6_1, label="y₀ = [1.5, 3.0]", linecolor=:blue,  linewidth=3)
    plot!(x7, y7_1, label="y₀ = [3.0, 1.0]", linecolor=:red,   linewidth=3)
    plot!(x8, y8_1, label="y₀ = [3.5, 2.5]", linecolor=:cyan,  linewidth=3)
    plot!(size=(500, 500))
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!(string("Stiff Brusselator"))
    xlabel!("x")
    ylabel!("y₁")
    figName = string("stiffBrusselatorY1.png")
    figPath = string(dirPath, figName)
    savefig(figPath)

    plot(x5, y5_2,  label="y₀ = [0.1, 0.1]", linecolor=:black, linewidth=3)
    plot!(x6, y6_2, label="y₀ = [1.5, 3.0]", linecolor=:blue,  linewidth=3)
    plot!(x7, y7_2, label="y₀ = [3.0, 1.0]", linecolor=:red,   linewidth=3)
    plot!(x8, y8_2, label="y₀ = [3.5, 2.5]", linecolor=:cyan,  linewidth=3)
    plot!(size=(500, 500))
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!(string("Stiff Brusselator"))
    xlabel!("x")
    ylabel!("y₂")
    figName = string("stiffBrusselatorY2.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(x5, err5,  label="y₀ = [0.1, 0.1]", linecolor=:black, linewidth=3)
    plot!(x6, err6, label="y₀ = [1.5, 3.0]", linecolor=:blue,  linewidth=3)
    plot!(x7, err7, label="y₀ = [3.0, 1.0]", linecolor=:red,   linewidth=3)
    plot!(x8, err8, label="y₀ = [3.5, 2.5]", linecolor=:cyan,  linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(yscale=:log10, minorgrid=true)
    ylims!(1e-5*tol, 10*tol)
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!("Stiff Brusselator: Errors")
    xlabel!("x")
    ylabel!("Truncation Error, ε")
    figName = string("stiffBrusselatorError.png")
    figPath = string(dirPath, figName)
    savefig(figPath)

    return nothing
end # run

end # brusselator