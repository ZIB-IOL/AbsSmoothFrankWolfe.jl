using AbsSmoothFrankWolfe
using FrankWolfe
using LinearAlgebra
using JuMP
using HiGHS

import MathOptInterface
const MOI = MathOptInterface

# Wong 2
 function f(x)
 	return max(x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45,
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(3*(x[1]-2)^2+4*(x[2]-3)^2+2*x[3]^2-7*x[4]-120),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(5*x[1]^2+8*x[2]+(x[3]-6)^2-2*x[4]-40),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(0.5*(x[1]-8)^2+2*(x[2]-4)^2+3*x[5]^2-x[6]-30),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(x[1]^2+2*(x[2]-2)^2-2*x[1]*x[2]+14*x[5]-6*x[6]),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(4*x[1]+5*x[2]-3*x[7]+9*x[8]-105),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(10*x[1]-8*x[2]-17*x[7]+2*x[8]),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(-3*x[1]+6*x[2]+12*(x[9]-8)^2-7*x[10]),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(-8*x[1]+2*x[2]+5*x[9]-2*x[10]-12))
 end
  
# evaluation point x_base
x_base = [2.0,3.0,5.0,5.0,1.0,2.0,7.0,3.0,6.0,10.0]
n = length(x_base)
 
lb_x = [-10.0 for in in x_base] 
ub_x = [10.0 for in in x_base]

# call the abs-linear form of f
abs_normal_form = abs_linear(x_base,f)  
s = abs_normal_form.num_switches

# gradient formula in terms of abs-linearization
function grad!(storage, x)
    c = ones(n+s)
    #c = vcat(alf_a', alf_b'.* sigma_z)
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
    max_iteration=1e7
)
@show f(x_base)
println("x_final = ", x_base)

