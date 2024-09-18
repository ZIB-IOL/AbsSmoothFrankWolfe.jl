module AbsSmoothFrankWolfe

using FrankWolfe
using SparseArrays
using LinearAlgebra
using JuMP
using HiGHS
using ADOLC


import MathOptInterface
const MOI = MathOptInterface

export abs_linear, AbsSmoothLMO, as_frank_wolfe, signature_vec

include("aasm.jl")
include("as_frank_wolfe.jl")
include("abs_linear.jl")
include("abs_lmo.jl")

end
