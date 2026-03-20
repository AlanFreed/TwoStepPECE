"""
    author:  Alan Freed
    date:    March 18, 2026
    
    To illustrate this class of problems, consider a vibration model for a car 
    in three degrees of freedom: bounce, pitch and roll, all measured at the 
    center of gravity of a car and its driver.  This example simulates an FSAE 
    race car.
    
    x = {b, p, r}^T   where  b = bounce,  p = pitch,  r = roll
    v = {db, dp, dr}^T
    a = {d2b, d2p, d2r}^T  and this is given by the equation:
    
    a = M^{-1} [fFn(t) - C*v - K*x]

    Bounce is in feet, while pitch and roll are in radians, per FSAE rules. 
    Bounce is positive downward (towards the ground).  Pitch is positive when 
    the nose is up and the tail is down.  Roll is positive when the driver is
    up and the passenger is down.
    
    The mass matrix  M  for this 3 degree-of-freedom (DOF) problem is
            / m   0  0  \                        / 1/m  0    0   \
        M = | 0  Jy  0  |    so that    M^{-1} = |  0  1/Jy      |
            \ 0   0  Jx /                        \  0   0   1/Jx /
    where  m  is the mass of the vehicle in slugs, while  Jx  and  Jy  are the
    moments of inertia in units of  ft.lbs/(rad/sec^2)  about the x and y axes, 
    per FSAE rules.  
    
    The symmetric damping matrix  C  for this 3 DOF car simulation is
            / c11 c12 c13 \
        C = | c12 c22 c23 |
            \ c13 c23 c33 /
    wherein
        c11 = c1 + c2 + c3 + c4
        c12 = −(c1 + c2) lf + (c3 + c4) lr
        c13 = −(c1 − c2) rf + (c3 − c4) rr 
        c22 = (c1 + c2) lf^2 + (c3 + c4) lr^2
        c23 = -(c1 − c2) lf rf + (c3 − c4) lr rr
        c33 = (c1 + c2) rf^2 + (c3 + c4) rr^2
    where  c1  is the damping of the driver front shock absorber,  c2  is the
    damping of the passenger front shock absorber,  c3  is the damping of the
    passenger rear shock absorber,  c4  is the damping of the driver rear 
    shock absorber, all of which have units of lb/(ft/sec).  Parameter  lf  
    is the distance from the center of gravity (CG) to the front axle,  lr 
    is the distance from the CG to the rear axle,  rf  is the radial distance 
    from the axial centerline (CL) to the center of the tire patch at the front 
    axle, and  rr  is the radial distance from the CL to the center of the tire 
    patch at the rear axle, with distances being in feet, per FSAE rules.
    
    The symmetric stiffness matrix  K  for this 3 DOF car simulation is
            / k11 k12 k13 \
        K = | k12 k22 k23 |
            \ k13 k23 k33 /
    wherein
        k11 = k1 + k2 + k3 + k4
        k12 = −(k1 + k2) lf + (k3 + k4) lr
        k13 = −(k1 − k2) rf + (k3 − k4) rr 
        k22 = (k1 + k2) lf^2 + (k3 + k4) lr^2
        k23 = -(k1 − k2) lf rf + (k3 − k4) lr rr
        k33 = (k1 + k2) rf^2 + (k3 + k4) rr^2
    where  k1  is the stiffness of the driver front spring,  k2  is the
    stiffness of the passenger front spring,  k3  is the stiffness of the
    passenger rear spring,  k4  is the stiffness of the driver rear spring,
    all of which have units of lb/ft, per FSAE rules.  The other parameters 
    are as defined for the damping matrix.
        
    The forcing function  fFn  for thie 3 DOF car simulation is
              / f1 \
        fFn = | f2 |
              \ f3 /
    wherein
       f1 = w − c1 Ṙ1 − c2 Ṙ2 − c3 Ṙ3 − c4 Ṙ4 
          −  k1 R1 − k2 R2 − k3 R3 − k4 R4
       f2 = (c1 Ṙ1 + c2 Ṙ2 + k1 R1 + k2 R2) lf 
          − (c3 Ṙ3 + c4 Ṙ4 + k3 R3 + k4 R4) lr
       f3 = (c1 Ṙ1 − c2 Ṙ2 + k1 R1 − k2 R2) rf 
          − (c3 Ṙ3 − c4 Ṙ4 + k3 R3 − k4 R4) rr
    where  w  is the weight of the car in pounds,  R1, R2, R3, R4  are the 
    upward displacements of the roadway, which are functions of time, and
    Ṙ1, Ṙ2, Ṙ3, Ṙ4  are their rates of change, which are also functions of 
    time.  Units are in ft and ft/sec, respectively.  The other parameters 
    are as defined for the damping and stiffness matrices.
    
    Representative parameters for a typical FSAE race car with driver are:
        m = 14       in slug = lbf.sec^s/ft
        w = 450      in lbf
        Jx = 20      in ft.lbf/(rad/sec^2)
        Jy = 45      in ft.lbf/(rad/sec^2)
        lf = 3.2     in ft
        lr = 1.8     in ft
        rf = 2.1     in ft
        rr = 2       in ft
        c1 = 120      in lbf/(in/sec)
        c2 = 120      in lbf/(in/sec)
        c3 = 180      in lbf/(in/sec)
        c4 = 180      in lbf/(in/sec)
        k1 = 1800     in lbf/ft
        k2 = 1800     in lbf/ft
        k3 = 3600     in lbf/ft
        k4 = 3600     in lbf/ft
$$
\begin{aligned}
    m & = 14       & & \text{slug} \\
    w & = 450      & & \text{lbf} \\
    J_x & = 20     & & \text{ft.lbf/(rad/sec^2)} \\
    J_y & = 45     & & \text{ft.lbf/(rad/sec^2)} \\
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
"""
module fsaeRaceCars

end # fsaeRaceCars