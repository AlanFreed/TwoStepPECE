"""
-------------------------------------------------------------------------------

    author:  Alan Freed
    date:    Mar 18 2026
    updated: Jun 18 2026
    
    To illustrate this class of problems, consider a vibration model for a car 
    in three degrees of freedom: heave, pitch and roll, all measured at the 
    center of gravity of a car and its driver.  This example simulates an FSAE 
    race car.
    
    x = {h, p, r}^T   where  h = heave,  p = pitch,  r = roll
    v = {dh/dt, dp/dt, dr/dt}^T
    a = {d²h/dt², d²p/dt², d²r/dt²}^T  and this is given by the equation:
    
    a = M^{-1} [ff(t) - C*v - K*x]

    Heave is in feet, while pitch and roll are in radians, per FSAE rules. 
    Heave is positive downward (towards the ground). Pitch is positive when 
    the nose is up and the tail is down. Roll is positive when the driver is
    up and the passenger is down.
    
    The mass matrix  M  for this 3 degree-of-freedom (DOF) problem is
            ⌈ m   0  0  ⌉                        ⌈ 1/m  0    0   ⌉
        M = | 0  Jᵧ  0  |    so that    M^{-1} = |  0  1/Jᵧ      |
            ⌊ 0   0  Jₓ ⌋                        ⌊  0   0   1/Jₓ ⌋
    where m is the mass of the vehicle in slugs, while Jₓ and Jᵧ are the
    moments of inertia in units of slugs.ft² about the X and Y axes, 
    per FSAE rules.  The X axis is longitudinal, the Y axis is radial, 
    and the Z axis is verticle.  They form a right-handed triad.
    
    The symmetric damping matrix C for this 3 DOF car simulation is
            ⌈ C₁₁ C₁₂ C₁₃ ⌉
        C = | C₁₂ C₂₂ C₂₃ |
            ⌊ C₁₃ C₂₃ C₃₃ ⌋
    wherein
        C₁₁ = c₁ + c₂ + c₃ + c₄
        C₁₂ = −(c₁ + c₂) lf + (c₃ + c₄) lr
        C₁₃ = −(c₁ − c₂) rf + (c₃ − c₄) rr 
        C₂₂ = (c₁ + c₂) lf² + (c₃ + c₄) lr²
        C₂₃ = -(c₁ − c₂) lf rf + (c₃ − c₄) lr rr
        C₃₃ = (c₁ + c₂) rf² + (c₃ + c₄) rr²
    where c₁ is the damping of the driver-front shock absorber, c₂ is the
    damping of the passenger-front shock absorber, c₃ is the damping of the
    passenger-rear shock absorber, c₄ is the damping of the driver-rear 
    shock absorber, all of which have units of lbf/(ft/sec).  Parameter lf  
    is the distance from the center of gravity (CG) to the front axle, lr 
    is the distance from the CG to the rear axle, rf  is the radial distance 
    from the axial centerline (CL) to the center of the tire patch at the front 
    axle, and rr is the radial distance from the CL to the center of the tire 
    patch at the rear axle, with distances being in feet, per FSAE rules.
    
    The symmetric stiffness matrix K for this 3 DOF car simulation is
            ⌈ K₁₁ K₁₂ K₁₃ ⌉
        K = | K₁₂ K₂₂ K₂₃ |
            ⌊ K₁₃ K₂₃ K₃₃ ⌋
    wherein
        K₁₁ = k₁ + k₂ + k₃ + k₄
        K₁₂ = −(k₁ + k₂) lf + (k₃ + k₄) lr
        K₁₃ = −(k₁ − k₂) rf + (k₃ − k₄) rr 
        K₂₂ = (k₁ + k₂) lf² + (k₃ + k₄) lr²
        K₂₃ = -(k₁ − k₂) lf rf + (k₃ − k₄) lr rr
        K₃₃ = (k₁ + k₂) rf² + (k₃ + k₄) rr²
    where k₁ is the stiffness of the driver-front spring, k₂ is the
    stiffness of the passenger-front spring, k₃ is the stiffness of the
    passenger-rear spring, k₄ is the stiffness of the driver-rear spring,
    all of which have units of lb/ft, per FSAE rules. The other parameters 
    are as defined for the damping matrix.
    
    The forcing function ff for this 3 DOF car simulation is
             ⌈ ff₁ ⌉
        ff = | ff₂ |
             ⌊ ff₃ ⌋
    wherein
             ⌈ w ⌉   ⌈  c₁     c₂     c₃     c₄    ⌉ ⌈ Ṙ₁ ⌉
        ff = | 0 | - | -c₁ lf -c₂ lf  c₃ lr  c₄ lr | | Ṙ₂ |
             ⌊ 0 ⌋   ⌊ -c₁ rf  c₂ rf  c₃ rr -c₄ rr ⌋ | Ṙ₃ |
                                                     ⌊ Ṙ₄ ⌋
                     ⌈  k₁     k₂     k₃     k₄    ⌉ ⌈ R₁ ⌉
                   - | -k₁ lf -k₂ lf  k₃ lr  k₄ lr | | R₂ |
                     ⌊ -k₁ rf  k₂ rf  k₃ rr -k₄ rr ⌋ | R₃ |
                                                     ⌊ R₄ ⌋
    where w is the weight of the car in pounds, R₁, R₂, R₃, R₄ are the 
    upward displacements of the roadway, which are functions of time, and
    Ṙ₁, Ṙ₂, Ṙ₃, Ṙ₄ are their rates of change, which are also functions of 
    time. Units are in ft and ft/sec, respectively. These are the control
    variables and their time rates-of-change.  The other parameters are
    as defined for the damping and stiffness matrices.
       
    Given the general form for a second-order ODE of:
        x′ = fₓ(model, t, x, y)         controls: a first-order  ODE
        y″ = fᵧ(model, t, x, y, y′)    responses: a second-order ODE
    an instance of `FsaeRaceCar` is the `model,` `fₓ` is described by the
    function `roadway,` while `fᵧ` is described by the function `vehicle.`
    The control variables `x` define a vector {R₁, R₂, R₃, R₄}ᵀ, while the
    response variables `y` define a vector {h, p, r}ᵀ.
    
    Representative parameters for a typical FSAE race car with driver are:
        m = 14       in slug = lbf.sec²/ft
        w = 450      in lbf
        Jₓ = 20      in slug.ft²
        Jᵧ = 45      in slub.ft²
        lf = 3.2     in ft
        lr = 1.8     in ft
        rf = 2.1     in ft
        rr = 2       in ft
        c₁ = 120     in lbf.sec/ft
        c₂ = 120     in lbf.sec/ft
        c₃ = 180     in lbf.sec/ft
        c₄ = 180     in lbf.sec/ft
        k₁ = 1800    in lbf/ft
        k₂ = 1800    in lbf/ft
        k₃ = 3600    in lbf/ft
        k₄ = 3600    in lbf/ft
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
    
const
    mph2fps = 1.467     # 1 mph = 1.467 ft./sec.
    speed   = 10mph2fps # speed of car, ft./sec.

# Describe the car.

struct FsaeRaceCar <: TwoStepPECE.Model
    M::Matrix{Float64}      # the mass matrix
    C::Matrix{Float64}      # the damping matrix
    K::Matrix{Float64}      # the stiffness matrix
    ffF::Vector{Float64}    # the static contribution to the forcing function
    ffC::Matrix{Float64}    # the damping matrix for the forcing function
    ffK::Matrix{Float64}    # the stiffness matrix for the forcing function
        
    function FsaeRaceCar(m::Float64,  w::Float64,  Jₓ::Float64, Jᵧ::Float64,
                         lf::Float64, lr::Float64, rf::Float64, rr::Float64,
                         c₁::Float64, c₂::Float64, c₃::Float64, c₄::Float64,
                         k₁::Float64, k₂::Float64, k₃::Float64, k₄::Float64)
        # `m`  is the mass of car and driver
        # `w`  is the weight of car and driver
        # `Jₓ` is the moment of inertia resisting roll
        # `Jᵧ` is the moment of inertia resisting pitch
        # `lf` is the distance from CG to front axle
        # `lr` is the distance from CG to rear axle
        # `rf` is the distance from CL to center of tire patch at front axle
        # `rr` is the distance from CL to center of tire patch at rear axle
        # `c₁` is the damping from shock absorber at driver front
        # `c₂` is the damping from shock absorber at passenger front
        # `c₃` is the damping from shock absorber at passenger rear
        # `c₄` is the damping from shock absorber at driver rear
        # `k₁` is the stiffness from spring at driver front
        # `k₂` is the stiffness from spring at passenger front
        # `k₃` is the stiffness from spring at passenger rear
        # `k₄` is the stiffness from spring at driver rear
        
        # Create the mass matrix.
        M = zeros(Float64, 3, 3)
        M[1,1] = m
        M[2,2] = Jᵧ
        M[3,3] = Jₓ 
        # Create the damping matrix.
        C = zeros(Float64, 3, 3)
        C[1,1] = c₁ + c₂ + c₃ + c₄
        C[1,2] = -(c₁ + c₂) * lf + (c₃ + c₄) * lr
        C[1,3] = -(c₁ - c₂) * rf + (c₃ - c₄) * rr
        C[2,1] = C[1,2]
        C[2,2] = (c₁ + c₂) * lf^2 + (c₃ + c₄) * lr^2
        C[2,3] = -(c₁ - c₂) * lf * rf  + (c₃ - c₄) * lr * rr
        C[3,1] = C[1,3]
        C[3,2] = C[2,3]
        C[3,3] = (c₁ + c₂) * rf^2 + (c₃ + c₄) * rr^2
        # Create the stiffness matrix.
        K = zeros(Float64, 3, 3)
        K[1,1] = k₁ + k₂ + k₃ + k₄
        K[1,2] = -(k₁ + k₂) * lf + (k₃ + k₄) * lr
        K[1,3] = -(k₁ - k₂) * rf + (k₃ - k₄) * rr 
        K[2,1] = K[1,2]
        K[2,2] = (k₁ + k₂) * lf^2 + (k₃ + k₄) * lr^2
        K[2,3] = -(k₁ - k₂) * lf * rf + (k₃ - k₄) * lr * rr
        K[3,1] = K[1,3]
        K[3,2] = K[2,3]
        K[3,3] = (k₁ + k₂) * rf^2 + (k₃ + k₄) * rr^2
        # Create the static force vector for the forcing function.
        ffF = zeros(Float64, 3)
        ffF[1] = w
        # Create the damping matrix for the forcing function.
        ffC = zeros(Float64, 3, 4)
        ffC[1,1] = c₁
        ffC[1,2] = c₂
        ffC[1,3] = c₃
        ffC[1,4] = c₄
        ffC[2,1] = -c₁ * lf
        ffC[2,2] = -c₂ * lf
        ffC[2,3] = c₃ * lr
        ffC[2,4] = c₄ * lr
        ffC[3,1] = -c₁ * rf
        ffC[3,2] = c₂ * rf
        ffC[3,3] = c₃ * rr
        ffC[3,4] = -c₄ * rr
        # Create the stiffness matrix for the forcing function.
        ffK = zeros(Float64, 3, 4)
        ffK[1,1] = k₁
        ffK[1,2] = k₂
        ffK[1,3] = k₃
        ffK[1,4] = k₄
        ffK[2,1] = -k₁ * lf
        ffK[2,2] = -k₂ * lf
        ffK[2,3] = k₃ * lr
        ffK[2,4] = k₄ * lr
        ffK[3,1] = -k₁ * rf
        ffK[3,2] = k₂ * rf
        ffK[3,3] = k₃ * rr
        ffK[3,4] = -k₄ * rr
        
        new(M, C, K, ffF, ffC, ffK)
    end
end # FsaeRaceCar
    
# Create a roadway for the car to traverse.

struct Bump
    height::Float64    # verticle height of a bump in ft
    width::Float64     # horizontal width of a bump in ft
    top::Float64       # width of flat region at top of a bump in ft
    
    function Bump(height::Float64, width::Float64, top::Float64)
        new(height, abs(width), abs(top))
    end
end # Bump

# Characteristics of a speed bump.
bump = Bump(1.0/6.0, 2.0, 0.5)  # Bump(height, width, top)

function speedBump(position::Float64)::Float64
    if position < 0.0 || position > bump.width
        Ṙ = 0.0
    else
        if position < (bump.width - bump.top) / 2.0
            φ = 2π * position / (bump.width - bump.top)
        elseif position < (bump.width + bump.top) / 2.0
            φ = π
        else
            φ = 2π * (position - bump.top) / (bump.width - bump.top)
        end
        Ṙ = (π * bump.height / (bump.width - bump.top)) * sin(φ) * speed
    end
    return Ṙ
end # speedBump
    
function mogul(position::Float64)::Float64
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
    return speedBump(location)
end # mogul

function roadwayDF(car::FsaeRaceCar, time::Float64)::Float64
    position = speed * time
    Ṙ = mogul(position)
    return Ṙ
end # roadwayDF

function roadwayPF(car::FsaeRaceCar, time::Float64)::Float64
    offset = 0.5         # distance passenger side trails the driver's side
    
    position = speed * time
    position = position - offset
    Ṙ = mogul(position)
    return Ṙ
end # roadwayDF

function roadwayPR(car::FsaeRaceCar, time::Float64)::Float64
    lf = -car.ffK[2,1] / car.ffK[1,1]  # length CG to front axel
    lr =  car.ffK[2,3] / car.ffK[1,3]  # length CG to real axel
    wheelbase = lf + lr
    offset = 0.5         # distance passenger side trails the driver's side

    position = speed * time
    position = position - wheelbase - offset
    Ṙ = mogul(position)
    return Ṙ
end # roadwayDF

function roadwayDR(car::FsaeRaceCar, time::Float64)::Float64
    lf = -car.ffK[2,1] / car.ffK[1,1]  # length CG to front axel
    lr =  car.ffK[2,3] / car.ffK[1,3]  # length CG to real axel
    wheelbase = lf + lr

    position = speed * time
    position = position - wheelbase
    Ṙ = mogul(position)
    return Ṙ
end # roadwayDF

# Function roadway is the model's fₓ, viz., x′ = fₓ(model, t, x, y).
function roadway(car::FsaeRaceCar, t::Float64, x::Vector{Float64}, y::Vector{Float64})::Vector{Float64}
    # arguments
    #    car  the car's parameters.
    #    t    the independent variable, i.e., time.
    #    x    the control variables, i.e., displacement of the wheels.
    #    y    the response variables, viz, displacement of the vehicle.
    # returned
    #    x′   Vertical velocities of the four wheels caused by a roadway.
    Ṙ₁ = roadwayDF(car, t)  # driver front
    Ṙ₂ = roadwayPF(car, t)  # passenger front
    Ṙ₃ = roadwayPR(car, t)  # passenger rear
    Ṙ₄ = roadwayDR(car, t)  # driver rear
    
    Ṙ  = [Ṙ₁ , Ṙ₂, Ṙ₃, Ṙ₄]
    
    return Ṙ # which is x′
end # roadway

# Function vehicle is the model's fᵧ, viz., y″ = fᵧ(model, t, x, y, y′).
function vehicle(car::FsaeRaceCar, t::Float64, x::Vector{Float64}, y::Vector{Float64}, y′::Vector{Float64})::Vector{Float64}
    # arguments
    #    car  the car's parameters.
    #    t    the independent variable, i.e., time.
    #    x    the control variables, i.e., displacement of the wheels.
    #    y    the response variables, viz, displacement of the vehicle.
    #    y′   rate of response variables, viz, velocity of the vehicle.
    # returned
    #    y″   acceleration of the vehicle.
    
    # Forces caused by the control variables:
    # x  = {R₁, R₂, R₃, R₄}ᵀ
    # x′ = {Ṙ₁, Ṙ₂, Ṙ₃, Ṙ₄}ᵀ
    x′ = roadway(car, t, x, y) 
    ff = car.ffF - car.ffC*x′ - car.ffK*x

    # Forces caused by the response variables:
    # y  = {h, p, r}ᵀ           h = heave
    # y′ = {ḣ, ṗ, ṙ}ᵀ   where   p = pitch
    # y″ = {ḧ, p̈, r̈}ᵀ           r = roll
    y″ = car.M \ (ff - car.C*y′ - car.K*y)
    
    return y″
end # vehicle

function run()
    # Asign parameters that define a car.
    m  = 14.0     # mass in slugs
    w  = 450.0    # weight in lbs
    Jₓ = 20.0     # moment of inertia resisting roll  in ft.lbs/(rad/sec^2)
    Jᵧ = 45.0     # moment of inertia resisting pitch in ft.lbs/(rad/sec^2)
    lf = 3.2      # distance from CG to front axle in ft
    lr = 1.8      # distance from CG to rear  axle in ft
    rf = 2.1      # distance from CL to center tire patch at front axle in ft
    rr = 2.0      # distance from CL to center tire patch at rear  axle in ft
    c₁ = 120.0    # driver-front damping from shock absorber in lbs/(ft/sec)
    c₂ = 120.0    # passenger-front damping from shock absorber in lbs/(ft/sec)
    c₃ = 180.0    # passenger-rear damping from shock absorber in lbs/(ft/sec)
    c₄ = 180.0    # driver-rear damping from shock absorber in lbs/(ft/sec)
    k₁ = 1800.0   # driver-front spring stiffness in lbs/ft
    k₂ = 1800.0   # passenger-front spring stiffness in lbs/ft
    k₃ = 3600.0   # passenger-rear spring stiffness in lbs/ft
    k₄ = 3600.0   # driver-rear spring stiffness in lbs/ft
    car = FsaeRaceCar(m, w, Jₓ, Jᵧ, lf, lr, rf, rr, 
                      c₁, c₂, c₃, c₄, k₁, k₂, k₃, k₄)
                      
    # Properties for the integrator.
    tol = 0.0001   # upper bound on the local truncation error
    N = 500        # number of global steps
    T = 1.5        # time at the end of the run/analysis
    
    # Establish the initial conditions.
    t₀  = 0.0
    x₀  = zeros(Float64, 4)
    y₀  = car.K \ car.ffF
    x′₀ = roadway(car, t₀, x₀, y₀)
    y′₀ = zeros(Float64, 3)
    y″₀ = vehicle(car, t₀, x₀, y₀, y′₀)
    
    print("\nThe initial state is:\n")
    print("  z  = ", PF.toString(12y₀[1]), " in.\n")
    print("  z′ = ", PF.toString(12y′₀[1]), " in./sec.\n")
    print("  z″ = ", PF.toString(12y″₀[1]), " in./sec.²\n")
    print("  θ  = ", PF.toString(180y₀[2]/π), " deg.\n") 
    print("  θ′ = ", PF.toString(180y′₀[2]/π), " deg./sec.\n") 
    print("  θ″ = ", PF.toString(180y″₀[2]/π), " deg./sec.²\n") 
    print("  φ  = ", PF.toString(180y₀[3]/π), " deg.\n")
    print("  φ′ = ", PF.toString(180y′₀[3]/π), " deg./sec.\n")
    print("  φ″ = ", PF.toString(180y″₀[3]/π), " deg./sec.²\n")
    
    # Create the data vectors used to construct the figures.
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
    
    # Create the solver.
    pece = SecondOrderPECE(roadway, vehicle, car, N, t₀, T, x₀, y₀, y′₀, tol)
    
    # Begin populating the data vectors.
    t[1]  = PF.get(pece.t_prev)
    z[1]  = 12pece.y_prev[1]
    θ[1]  = 180pece.y_prev[2] / π
    φ[1]  = 180pece.y_prev[3] / π
    z′[1] = 12pece.y′_prev[1]
    θ′[1] = 180pece.y′_prev[2] / π
    φ′[1] = 180pece.y′_prev[3] / π
    z″[1] = 12pece.y″_prev[1]
    θ″[1] = 180pece.y″_prev[2] / π
    φ″[1] = 180pece.y″_prev[3] / π
    
    i = 1
    while pece.n < pece.N
        advance!(pece)
        if pece.atNode == true
            i = i + 1
            t[i]  = PF.get(pece.t_curr)
            ε[i]  = PF.get(pece.ε_curr)
            z[i]  = 12pece.y_curr[1]
            θ[i]  = 180pece.y_curr[2] / π
            φ[i]  = 180pece.y_curr[3] / π
            z′[i] = 12pece.y′_curr[1]
            θ′[i] = 180pece.y′_curr[2] / π
            φ′[i] = 180pece.y′_curr[3] / π
            z″[i] = 12pece.y″_curr[1]
            θ″[i] = 180pece.y″_curr[2] / π
            φ″[i] = 180pece.y″_curr[3] / π
        end
    end # while
    ε[1] = ε[2]
    print("\nThe FSAE race car analysis ran with statistics:\n")
    print("   ", PF.toString(pece.steps), " steps taken with ",
          PF.toString(pece.repeats), " steps repeated\n")
    print("   of which ", PF.toString(pece.doubled), " were doubled and ",
          PF.toString(pece.halved), " were halved.\n")
    
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
    ylabel!("rate of heave, dz/dt  (in/sec)")
    figName = string("FSAE_z′.png")
    figPath = string(dirPath, figName)
    savefig(figPath)
    
    plot(t, z″, label="d²z/dt²", linecolor=:black, linewidth=3)
    plot!(size=(809, 500)) # (500ϕ, 500), ϕ is golden ratio
    plot!(legend=:topright)
    plot!(left_margin=3mm, right_margin=3mm, top_margin=3mm, bottom_margin=3mm)
    xlabel!("time, t  (sec)")
    ylabel!("rate of heave, d²z/dt²  (in/sec²)")
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
