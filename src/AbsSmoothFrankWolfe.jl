module AbsSmoothFrankWolfe

using FrankWolfe
using SparseArrays
using LinearAlgebra
using JuMP
using HiGHS
using ADOLC


import MathOptInterface
const MOI = MathOptInterface

include("aasm.jl")
include("as_frank_wolfe.jl")
include("abs_linear.jl")
include("abs_lmo.jl")

end
