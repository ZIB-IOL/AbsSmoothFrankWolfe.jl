using AbsSmoothFrankWolfe

using FrankWolfe

using LinearAlgebra
using JuMP
using HiGHS

import MathOptInterface
const MOI = MathOptInterface

# Chained Mifflin 2  
 function f(x)
     n = length(x)
     l_1 = [-x[i]+2*(x[i]^2+x[i+1]^2-1)+1.75*abs(x[i]^2+x[i+1]^2-1) for i in 1:n-1]
     return sum(l_1[i] for i in 1:n-1)
 end
 
 n=200
 x_base = ones(n)

n = length(x_base)

lb_x = [-3.0 for in in x_base] 
ub_x = [3.0 for in in x_base]

# call the abs-linear form of f
abs_normal_form = AbsSmoothFrankWolfe.abs_linear(x_base,f)
s = abs_normal_form.num_switches
# gradient formula in terms of abs-linearization
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
