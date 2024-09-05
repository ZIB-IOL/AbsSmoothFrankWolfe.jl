using FrankWolfe
using SparseArrays
using LinearAlgebra
using JuMP
using HiGHS
using ADOLC

import MathOptInterface
const MOI = MathOptInterface

include("../src/aasm.jl")
include("../src/as_frank_wolfe.jl")
include("../src/abs_linear.jl")
include("../src/abs_lmo.jl")

# DEM 
 function f(x)
 	return max(5*x[1]+x[2], -5*x[1]+x[2], x[1]^2+x[2]^2+4*x[2])
 end
 
 
# evaluation point x_base
x_base = [1.0,1.0]
n = length(x_base)
 
lb_x = [-5 for in in x_base] 
ub_x = [5 for in in x_base]

# call the abs-linear form of f
abs_normal_form = abs_linear(x_base,f)
#@show abs_normal_form

alf_a = abs_normal_form.Y
alf_b = abs_normal_form.J 
z = abs_normal_form.z  
s = abs_normal_form.num_switches

sigma_z = signature_vec(s,z)


# gradient formula in terms of abs-linearization
function grad!(storage, x)
    c = vcat(alf_a', alf_b'.* sigma_z)
    @. storage = c
end

# set bounds
o = Model(HiGHS.Optimizer)
MOI.set(o, MOI.Silent(), true)
@variable(o, lb_x[i] <= x[i=1:n] <= ub_x[i])

# initialise dual gap 
dualgap_asfw = Inf

# abs-smooth lmo
lmo_as = AbsSmoothLMO(o, x_base, f, n, s, lb_x, ub_x, dualgap_asfw)

# define termination criteria

# In case we want to stop the frank_wolfe algorithm prematurely after a certain condition is met,
# we can return a boolean stop criterion `false`.
# Here, we will implement a callback that terminates the algorithm if ASFW Dual gap < eps.
function make_termination_callback(state)
 return function callback(state,args...)
  return state.lmo.dualgap_asfw[1] > 1e-2
 end
end

callback = make_termination_callback(FrankWolfe.CallbackState)

# call abs-smooth-frank-wolfe
x, v, primal, dual_gap, traj_data = as_frank_wolfe(
    f, 
    grad!, 
    lmo_as, 
    x_base;
    gradient = ones(n+s),
    line_search = FrankWolfe.FixedStep(1.0),
    callback=callback,
    verbose=true,
    max_iteration=6078
)

@show x_base

