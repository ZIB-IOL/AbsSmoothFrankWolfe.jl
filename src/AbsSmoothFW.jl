module AbsSmoothFW

using FrankWolfe
using LinearAlgebra
using SparseArrays
using JuMP
using HiGHS
using ADOLC

import MathOptInterface
const MOI = MathOptInterface

include("aasm.jl")
include("abs_linear.jl")

end
