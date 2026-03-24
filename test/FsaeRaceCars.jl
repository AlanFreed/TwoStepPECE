"""
-------------------------------------------------------------------------------

    author:  Alan Freed
    date:    March 18, 2026
    
    To illustrate this class of problems, consider a vibration model for a car 
    in three degrees of freedom: heave, pitch and roll, all measured at the 
    center of gravity of a car and its driver.  This example simulates an FSAE 
    race car.
    
    x = {h, p, r}^T   where  h = heave,  p = pitch,  r = roll
    v = {dh/dt, dp/dt, dr/dt}^T
    a = {d²h/dt², d²p/dt², d²r/dt²}^T  and this is given by the equation:
    
    a = M^{-1} [fFn(t) - C*v - K*x]

    Heave is in feet, while pitch and roll are in radians, per FSAE rules. 
    Heave is positive downward (towards the ground). Pitch is positive when 
    the nose is up and the tail is down. Roll is positive when the driver is
    up and the passenger is down.
    
    The mass matrix  M  for this 3 degree-of-freedom (DOF) problem is
            ⌈ m   0  0  ⌉                        ⌈ 1/m  0    0   ⌉
        M = | 0  Jy  0  |    so that    M^{-1} = |  0  1/Jy      |
            ⌊ 0   0  Jx ⌋                        ⌊  0   0   1/Jx ⌋
    where m is the mass of the vehicle in slugs, while Jx and Jy are the
    moments of inertia in units of slugs.ft² about the x and y axes, per
    FSAE rules.  
    
    The symmetric damping matrix C for this 3 DOF car simulation is
            ⌈ c11 c12 c13 ⌉
        C = | c12 c22 c23 |
            ⌊ c13 c23 c33 ⌋
    wherein
        c11 = c1 + c2 + c3 + c4
        c12 = −(c1 + c2) lf + (c3 + c4) lr
        c13 = −(c1 − c2) rf + (c3 − c4) rr 
        c22 = (c1 + c2) lf^2 + (c3 + c4) lr^2
        c23 = -(c1 − c2) lf rf + (c3 − c4) lr rr
        c33 = (c1 + c2) rf^2 + (c3 + c4) rr^2
    where c1 is the damping of the driver front shock absorber, c2 is the
    damping of the passenger front shock absorber, c3 is the damping of the
    passenger rear shock absorber, c4 is the damping of the driver rear 
    shock absorber, all of which have units of lbf/(ft/sec).  Parameter lf  
    is the distance from the center of gravity (CG) to the front axle, lr 
    is the distance from the CG to the rear axle, rf  is the radial distance 
    from the axial centerline (CL) to the center of the tire patch at the front 
    axle, and rr is the radial distance from the CL to the center of the tire 
    patch at the rear axle, with distances being in feet, per FSAE rules.
    
    The symmetric stiffness matrix K for this 3 DOF car simulation is
            ⌈ k11 k12 k13 ⌉
        K = | k12 k22 k23 |
            ⌊ k13 k23 k33 ⌋
    wherein
        k11 = k1 + k2 + k3 + k4
        k12 = −(k1 + k2) lf + (k3 + k4) lr
        k13 = −(k1 − k2) rf + (k3 − k4) rr 
        k22 = (k1 + k2) lf^2 + (k3 + k4) lr^2
        k23 = -(k1 − k2) lf rf + (k3 − k4) lr rr
        k33 = (k1 + k2) rf^2 + (k3 + k4) rr^2
    where k1 is the stiffness of the driver front spring, k2 is the
    stiffness of the passenger front spring, k3 is the stiffness of the
    passenger rear spring, k4 is the stiffness of the driver rear spring,
    all of which have units of lb/ft, per FSAE rules. The other parameters 
    are as defined for the damping matrix.
        
    The forcing function  fFn  for thie 3 DOF car simulation is
              ⌈ f1 ⌉
        fFn = | f2 |
              ⌊ f3 ⌋
    wherein
       f1 = w − c1 Ṙ1 − c2 Ṙ2 − c3 Ṙ3 − c4 Ṙ4 
          −  k1 R1 − k2 R2 − k3 R3 − k4 R4
       f2 = (c1 Ṙ1 + c2 Ṙ2 + k1 R1 + k2 R2) lf 
          − (c3 Ṙ3 + c4 Ṙ4 + k3 R3 + k4 R4) lr
       f3 = (c1 Ṙ1 − c2 Ṙ2 + k1 R1 − k2 R2) rf 
          − (c3 Ṙ3 − c4 Ṙ4 + k3 R3 − k4 R4) rr
    where w is the weight of the car in pounds, R1, R2, R3, R4 are the 
    upward displacements of the roadway, which are functions of time, and
    Ṙ1, Ṙ2, Ṙ3, Ṙ4 are their rates of change, which are also functions of 
    time. Units are in ft and ft/sec, respectively. The other parameters 
    are as defined for the damping and stiffness matrices.
    
    Representative parameters for a typical FSAE race car with driver are:
        m = 14       in slug = lbf.sec²/ft
        w = 450      in lbf
        Jx = 20      in slug.ft²
        Jy = 45      in slub.ft²
        lf = 3.2     in ft
        lr = 1.8     in ft
        rf = 2.1     in ft
        rr = 2       in ft
        c1 = 120     in lbf.sec/ft
        c2 = 120     in lbf.sec/ft
        c3 = 180     in lbf.sec/ft
        c4 = 180     in lbf.sec/ft
        k1 = 1800    in lbf/ft
        k2 = 1800    in lbf/ft
        k3 = 3600    in lbf/ft
        k4 = 3600    in lbf/ft
"""
module fsaeRaceCars

using
    Measures,        # needed to pad white space (margins) around a plot
    PhysicalFields,
    Plots,
    TwoStepPECE
    
import
    PhysicalFields as PF
    
export
    run
    
struct FsaeRaceCar
    m::Float64     # mass of car and driver
    w::Float64     # weight of car and driver
    Jx::Float64    # moment of inertia resisting roll
    Jy::Float64    # moment of inertia resisting pitch
    lf::Float64    # distance from CG to front axle
    lr::Float64    # distance from CG to rear axle
    rf::Float64    # distance from CL to center of tire patch at front axle
    rr::Float64    # distance from CL to center of tire patch at rear axle
    c1::Float64    # damping from shock absorber at driver front
    c2::Float64    # damping from shock absorber at passenger front
    c3::Float64    # damping from shock absorber at passenger rear
    c4::Float64    # damping from shock absorber at driver rear
    k1::Float64    # stiffness from spring at driver front
    k2::Float64    # stiffness from spring at passenger front
    k3::Float64    # stiffness from spring at passenger rear
    k4::Float64    # stiffness from spring at driver rear
    
    function FsaeRaceCar(m::Float64, w::Float64, Jx::Float64, Jy::Float64,
                         lf::Float64, lr::Float64, rf::Float64, rr::Float64,
                         c1::Float64, c2::Float64, c3::Float64, c4::Float64,
                         k1::Float64, k2::Float64, k3::Float64, k4::Float64)
        new(m, w, Jx, Jy, lf, lr, rf, rr, c1, c2, c3, c4, k1, k2, k3, k4)
    end
end # FsaeRaceCar

function massMtx(car::FsaeRaceCar)::Matrix{Float64}
    M = zeros(Float64, 3, 3)
    M[1,1] = car.m
    M[2,2] = car.Jy
    M[3,3] = car.Jx 
    return M
end # massMtx

function dampingMtx(car::FsaeRaceCar)::Matrix{Float64}
    C = zeros(Float64, 3, 3)
    C[1,1] = car.c1 + car.c2 + car.c3 + car.c4
    C[1,2] = -(car.c1 + car.c2) * car.lf + (car.c3 + car.c4) * car.lr
    C[1,3] = -(car.c1 - car.c2) * car.rf + (car.c3 - car.c4) * car.rr
    C[2,1] = C[1,2]
    C[2,2] = (car.c1 + car.c2) * car.lf^2 + (car.c3 + car.c4) * car.lr^2
    C[2,3] = (-(car.c1 - car.c2) * car.lf * car.rf 
              + (car.c3 - car.c4) * car.lr * car.rr)
    C[3,1] = C[1,3]
    C[3,2] = C[2,3]
    C[3,3] = (car.c1 + car.c2) * car.rf^2 + (car.c3 + car.c4) * car.rr^2
    return C
end # dampingMtx

function stiffnessMtx(car::FsaeRaceCar)::Matrix{Float64}
    K = zeros(Float64, 3, 3)
    K[1,1] = car.k1 + car.k2 + car.k3 + car.k4
    K[1,2] = -(car.k1 + car.k2) * car.lf + (car.k3 + car.k4) * car.lr
    K[1,3] = -(car.k1 - car.k2) * car.rf + (car.k3 - car.k4) * car.rr 
    K[2,1] = K[1,2]
    K[2,2] = (car.k1 + car.k2) * car.lf^2 + (car.k3 + car.k4) * car.lr^2
    K[2,3] = (-(car.k1 - car.k2) * car.lf * car.rf 
              + (car.k3 - car.k4) * car.lr * car.rr)
    K[3,1] = K[1,3]
    K[3,2] = K[2,3]
    K[3,3] = (car.k1 + car.k2) * car.rf^2 + (car.k3 + car.k4) * car.rr^2
    return K
end # stiffnessMtx

struct Bump
    height::Float64    # verticle height of a bump in ft
    width::Float64     # horizontal width of a bump in ft
    top::Float64       # width of flat region at top of a bump in ft
    
    function Bump(height::Float64, width::Float64, top::Float64)
        new(height, abs(width), abs(top))
    end
end # Bump

function speedBump(bump::Bump, position::Float64, speed::Float64)::Tuple
    if position < 0.0 || position > bump.width
        R = 0.0
        Ṙ = 0.0
    else
        if position < (bump.width - bump.top) / 2.0
            φ = 2π * position / (bump.width - bump.top)
        elseif position < (bump.width + bump.top) / 2.0
            φ = π
        else
            φ = 2π * (position - bump.top) / (bump.width - bump.top)
        end
        R = (1.0 - cos(φ)) * bump.height / 2.0
        Ṙ = (π * bump.height / (bump.width - bump.top)) * sin(φ) * speed
    end
    return (R, Ṙ)
end # speedBump
    
function mogul(bump::Bump, position::Float64, speed::Float64)::Tuple
    if position ≥ 0.0 && position < bump.width
        location = position
    elseif position ≥ bump.width && position < 2bump.width
        location = position - bump.width
    elseif position ≥ 2bump.width && position < 3bump.width
        location = position - 2bump.width
    elseif position ≥ 3bump.width && position < 4bump.width
        location = position - 3bump.width
    elseif position ≥ 4bump.width && position < 5bump.width
        location = position - 4bump.width
    else
        location = -1.0
    end
    return speedBump(bump, location, speed)
end # mogul

function trajectory(time::Float64, speed::Float64)::Tuple
    mph2fps = 1.467
    speed = speed * mph2fps
    x = speed * time
    v = speed
    return (x, v)
end # trajectory

function roadwayDF(car::FsaeRaceCar, bump::Bump, time::Float64, speed::Float64)::Tuple
    (position, speed) = trajectory(time, speed)
    (R, Ṙ) = mogul(bump, position, speed)
    return (R, Ṙ)
end # roadwayDF

function roadwayPF(car::FsaeRaceCar, bump::Bump, time::Float64, speed::Float64)::Tuple
    offset = 0.5         # distance passenger side trails the driver's side
    
    (position, speed) = trajectory(time, speed)
    position = position - offset
    (R, Ṙ) = mogul(bump, position, speed)
    return (R, Ṙ)
end # roadwayDF

function roadwayPR(car::FsaeRaceCar, bump::Bump, time::Float64, speed::Float64)::Tuple
    offset = 0.5         # distance passenger side trails the driver's side
    wheelbase = car.lf + car.lr

    (position, speed) = trajectory(time, speed)
    position = position - wheelbase - offset
    (R, Ṙ) = mogul(bump, position, speed)
    return (R, Ṙ)
end # roadwayDF

function roadwayDR(car::FsaeRaceCar, bump::Bump, time::Float64, speed::Float64)::Tuple
    wheelbase = car.lf + car.lr

    (position, speed) = trajectory(time, speed)
    position = position - wheelbase
    (R, Ṙ) = mogul(bump, position, speed)
    return (R, Ṙ)
end # roadwayDF

function forcingFn(car::FsaeRaceCar, time::Float64, speed::Float64)::Vector{Float64}
    height = 1.0 / 6.0
    width = 2.0
    top = 0.5
    bump = Bump(height, width, top)
    
    (R1, Ṙ1) = roadwayDF(car, bump, time, speed)
    (R2, Ṙ2) = roadwayPF(car, bump, time, speed)
    (R3, Ṙ3) = roadwayPR(car, bump, time, speed)
    (R4, Ṙ4) = roadwayDR(car, bump, time, speed)
    
    F = Vector{Float64}(undef, 3)
    F[1] = (car.w - car.c1 * Ṙ1 - car.c2 * Ṙ2 - car.c3 * Ṙ3 - car.c4 * Ṙ4 
                  - car.k1 * R1 - car.k2 * R2 - car.k3 * R3 - car.k4 * R4)
    F[2] = ((car.c1 * Ṙ1 + car.c2 * Ṙ2 + car.k1 * R1 + car.k2 * R2) * car.lf 
           -(car.c3 * Ṙ3 + car.c4 * Ṙ4 + car.k3 * R3 + car.k4 * R4) * car.lr)
    F[3] = ((car.c1 * Ṙ1 - car.c2 * Ṙ2 + car.k1 * R1 - car.k2 * R2) * car.rf 
           -(car.c3 * Ṙ3 - car.c4 * Ṙ4 + car.k3 * R3 - car.k4 * R4) * car.rr)
    return F
end # forcingFn

function ICs(car::FsaeRaceCar)::Tuple
    f0 = zeros(Float64, 3)
    x0 = zeros(Float64, 3)
    v0 = zeros(Float64, 3)
    
    f0[1] = car.w
    
    K = stiffnessMtx(car)
    x0 = K \ f0
    
    v0[1] = 0.0
    v0[2] = 0.0
    v0[3] = 0.0
    return (x0, v0)
end # ICs

function acceleration(car::FsaeRaceCar, time::Float64, speed::Float64, x::Vector{Float64}, v::Vector{Float64})::Vector{Float64}
    a = Vector{Float64}(undef, 3)
    f = forcingFn(car, time, speed)
    M = massMtx(car)
    C = dampingMtx(car)
    Cv = C * v
    K = stiffnessMtx(car)
    Kx = K * x
    a = M \ (f - Kx - Cv)
    return a
end # acceleration

function run()
    # asign the parameters that define a car
    m  = 14.0     # mass in slugs
    w  = 450.0    # weight in lbs
    Jx = 20.0     # moment of inertia resisting roll  in ft.lbs/(rad/sec^2)
    Jy = 45.0     # moment of inertia resisting pitch in ft.lbs/(rad/sec^2)
    lf = 3.2      # distance from CG to front axle in ft
    lr = 1.8      # distance from CG to rear  axle in ft
    rf = 2.1      # distance from CL to center tire patch at front axle in ft
    rr = 2.0      # distance from CL to center tire patch at rear  axle in ft
    c1 = 120.0    # driver front damping from shock absorber in lbs/(ft/sec)
    c2 = 120.0    # passenger front damping from shock absorber in lbs/(ft/sec)
    c3 = 180.0    # passenger rear damping from shock absorber in lbs/(ft/sec)
    c4 = 180.0    # driver rear damping from shock absorber in lbs/(ft/sec)
    k1 = 1800.0   # driver front spring stiffness in lbs/ft
    k2 = 1800.0   # passenger front spring stiffness in lbs/ft
    k3 = 3600.0   # passenger rear spring stiffness in lbs/ft
    k4 = 3600.0   # driver rear spring stiffness in lbs/ft
    car = FsaeRaceCar(m, w, Jx, Jy, lf, lr, rf, rr, 
                      c1, c2, c3, c4, k1, k2, k3, k4)
                      
    # properties for the integrator
    tol = 0.0001   # upper bound on the local truncation error
    N = 500        # number of global steps
    T = 1.5        # time at the end of the run/analysis
    
    # establish the initial state
    t0 = 0.0
    speed = 10.0   # speed of the car in mph
    (x0, v0) = ICs(car)
    a0 = acceleration(car, t0, speed, x0, v0)
    
    print("\nThe static deflection is:\n")
    print("  z = ", PF.toString(12x0[1]), " inches\n")
    print("  θ = ", PF.toString(180x0[2]/π), " degrees\n") 
    print("  φ = ", PF.toString(180x0[3]/π), " degrees\n")
    
    function ode(x::Float64, y::Vector{Float64}, y′::Vector{Float64})
        y″ = acceleration(car, x, speed, y, y′)
        return y″
    end # ode
    solver = SecondOrderPECE(ode, N, t0, T, x0, v0, tol)
    
    t  = zeros(Float64, N+1)    # time
    ε  = zeros(Float64, N+1)    # local truncation error
    z  = zeros(Float64, N+1)    # heave in inches
    θ  = zeros(Float64, N+1)    # pitch in degrees
    φ  = zeros(Float64, N+1)    # roll in degrees
    z′ = zeros(Float64, N+1)    # rate of heave in in/sec
    θ′ = zeros(Float64, N+1)    # rate of pitch in deg/sec
    φ′ = zeros(Float64, N+1)    # rate of roll in deg/sec
    z″ = zeros(Float64, N+1)    # acceleration of heave in in/sec²
    θ″ = zeros(Float64, N+1)    # acceleration of pitch in deg/sec²
    φ″ = zeros(Float64, N+1)    # acceleration of roll in deg/sec²
    
    t[1]  = t0
    ε[1]  = tol
    z[1]  = 12x0[1]
    θ[1]  = 180x0[2] / π
    φ[1]  = 180x0[3] / π
    z′[1] = 12v0[1]
    θ′[1] = 180v0[2] / π
    φ′[1] = 180v0[3] / π
    a0    = ode(t0, x0, v0)
    z″[1] = 12a0[1]
    θ″[1] = 180a0[2] / π
    φ″[1] = 180a0[3] / π
    
    i = 1
    while solver.n < solver.N
        advance!(solver)
        if solver.atNode == true
            i = i + 1
            t[i]  = PF.get(solver.x_curr)
            ε[i]  = PF.get(solver.ε_curr)
            z[i]  = 12solver.y_curr[1]
            θ[i]  = 180solver.y_curr[2] / π
            φ[i]  = 180solver.y_curr[3] / π
            z′[i] = 12solver.y′_curr[1]
            θ′[i] = 180solver.y′_curr[2] / π
            φ′[i] = 180solver.y′_curr[3] / π
            z″[i] = 12solver.y″_curr[1]
            θ″[i] = 180solver.y″_curr[2] / π
            φ″[i] = 180solver.y″_curr[3] / π
        end
    end # while
    print("\nThe FSAE race car analysis ran with statistics:\n")
    print("   ", PF.toString(solver.steps), " steps taken with ",
          PF.toString(solver.repeats), " steps repeated\n")
    print("   of which ", PF.toString(solver.doubled), " were doubled and ",
          PF.toString(solver.halved), " were halved.\n")
    
    # set the graphics backend to GR
    ENV["QT_QPA_PLATFORM"] = "wayland"
    gr()
    
    plot(t, z, label="z", linecolor=:black, linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!("FSAE race car: heave")
    xlabel!("time, t  (sec)")
    ylabel!("heave, z  (in)")
    dirPath = string(pwd(), "/figures/")
    if !isdir(dirPath)
        mkdir(dirPath)
    end
    figName = string("FSAE_z.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(t, z′, label="dz/dt", linecolor=:black, linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    xlabel!("time, t  (sec)")
    ylabel!("rate of heave, z′  (in/sec)")
    figName = string("FSAE_z′.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(t, z″, label="d²z/dt²", linecolor=:black, linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    xlabel!("time, t  (sec)")
    ylabel!("rate of heave, z″  (in/sec²)")
    figName = string("FSAE_z″.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(t, θ,  label="θ", linecolor=:blue, linewidth=3)
    plot!(t, φ, label="φ", linecolor=:red,  linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!("FSAE race car: pitch θ and roll φ")
    xlabel!("time, t  (sec)")
    ylabel!("angles, θ and φ  (deg)")
    dirPath = string(pwd(), "/figures/")
    figName = string("FSAE_θ.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(t, θ′,  label="dθ/dt", linecolor=:blue, linewidth=3)
    plot!(t, φ′, label="dφ/dt", linecolor=:red,  linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    xlabel!("time, t  (sec)")
    ylabel!("angle rates, dθ/dt and dφ/dt  (deg/sec)")
    dirPath = string(pwd(), "/figures/")
    figName = string("FSAE_θ′.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(t, θ″,  label="d²θ/dt²", linecolor=:blue, linewidth=3)
    plot!(t, φ″, label="d²φ/dt²", linecolor=:red,  linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    xlabel!("time, t  (sec)")
    ylabel!("angle rates, d²θ/dt² and d²φ/dt²  (deg/sec²)")
    dirPath = string(pwd(), "/figures/")
    figName = string("FSAE_θ″.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(t, ε, label="ε", linecolor=:black, linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(yscale=:log10, minorgrid=true)
    ylims!(1e-5*tol, 10*tol)
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    title!("FSAE race car: truncation error")
    xlabel!("time, t  (sec)")
    ylabel!("truncation error, ε")
    figName = string("FSAE_Error.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
end # run

end # FsaeRaceCars
