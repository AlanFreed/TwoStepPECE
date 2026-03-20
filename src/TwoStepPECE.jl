#=
Created on Sun 22 Feb 2026
Updated on Fri 20 Mar 2026
=#

module TwoStepPECE

using
    JSON3,
    LinearAlgebra,
    PhysicalFields,
    StructTypes
    
import
    LinearAlgebra  as LA,
    PhysicalFields as PF

export
    # abstract type
    PECE,
    
    # concrete types
    FirstOrderPECE,
    SecondOrderPECE,
    
    # functions
    toFile,
    fromFile,
    
    # solver functions
    advance!
#=
-------------------------------------------------------------------------------
=#

abstract type PECE end

#=
-------------------------------------------------------------------------------
=#

include("FirstOrderPECE.jl")

include("SecondOrderPECE.jl")

end # module TwoStepPECE