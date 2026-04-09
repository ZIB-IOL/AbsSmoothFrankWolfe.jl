using AbsSmoothFrankWolfe
using FrankWolfe
using LinearAlgebra
using JuMP
using HiGHS


import MathOptInterface
const MOI = MathOptInterface

# MAXQ

 # forgot to square x_x
 function f(x)
 	return max(x[1]^2,x[2]^2,x[3]^2,x[4]^2,x[5]^2,x[6]^2,x[7]^2,x[8]^2,x[9]^2,x[10]^2,x[11]^2,x[12]^2,x[13]^2,x[14]^2,x[15]^2,x[16]^2,x[17]^2,x[18]^2,x[19]^2,x[20]^2)
 end
 
 x_base = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,-11.0,-12.0,-13.0,-14.0,-15.0,-16.0,-17.0,-18.0,-19.0,-20]
n = length(x_base)
 
lb_x = [-20.0 for in in x_base] 
ub_x = [20.0 for in in x_base]

# call the abs-linear form of f
abs_normal_form = AbsSmoothFrankWolfe.abs_linear(x_base,f)
s = abs_normal_form.num_switches
# container for gradient - stores n+s instances
function grad!(storage, x)
    c = ones(n+s)
    @. storage = c
end

# set bounds
o = Model(HiGHS.Optimizer)
MOI.set(o, MOI.Silent(), true)
@variable(o, lb_x[i] <= x[i=1:n] <= ub_x[i])

# initialise dual gap 
dualgap_asfw = Inf

# abs-smooth lmo
lmo_as = AbsSmoothFrankWolfe.AbsSmoothLMO(o, x_base, f, n, s, lb_x, ub_x, dualgap_asfw)

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
x, v, primal, dual_gap, traj_data = AbsSmoothFrankWolfe.as_frank_wolfe(
    f, 
    grad!, 
    lmo_as, 
    x_base;
    gradient = ones(n+s),
    line_search = FrankWolfe.FixedStep(1.0),
    callback=callback,
    verbose=true,
    max_iteration=1e7
)
println("Grand Total Simplex Steps: ", lmo_as.total_simplex_count)
println("x_final = ", x_base)

