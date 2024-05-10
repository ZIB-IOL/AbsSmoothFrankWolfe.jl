module AbsSmoothFW

using FrankWolfe
using LinearAlgebra
using SparseArrays
using JuMP
using HiGHS
using ADOLC

import MathOptInterface
const MOI = MathOptInterface

include("as_frank_wolfe.jl")
include("aasm.jl")
include("abs_linear.jl")
include("Abs_LMO.jl")

end
