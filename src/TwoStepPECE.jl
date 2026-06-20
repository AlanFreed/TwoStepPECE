#=
-------------------------------------------------------------------------------
Created on Sun 22 Feb 2026
Updated on Tue 16 Jun 2026
=#

module TwoStepPECE

using
    LinearAlgebra,
    PhysicalFields
    
import
    LinearAlgebra  as LA,
    PhysicalFields as PF

export
    # abstract types
    Model,
    PECE,
    
    # concrete types
    FirstOrderPECE,
    SecondOrderPECE,
    
    # solver functions
    advance!
    
    # PECE solvers cannot be made persistent because their data structures
    # contain functions, which cannot be made persistent.
#=
-------------------------------------------------------------------------------
=#

abstract type Model end
abstract type PECE end

#=
-------------------------------------------------------------------------------
=#

include("FirstOrderPECE.jl")

include("SecondOrderPECE.jl")

end # module TwoStepPECE