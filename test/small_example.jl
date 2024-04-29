include("../src/AbsSmooth_FW.jl")
include("../src/adolc_call.jl")

# CB3
 function f(x)
 	return max(x[1]^4+x[2]^2, (2-x[1])^2+(2-x[2])^2, 2*exp(x[2]-x[1]))
 end

# initialization
 x_base = [2.0,2.0]
 lb_x = [-5, -5]
 ub_x = [5, 5]

n = length(x_base)

# call the abs-linear form of f
abs_normal_form = call_adolc(x_base, f) 
alf_a = abs_normal_form.Y 
alf_b = reshape(abs_normal_form.J, size(abs_normal_form.J)[2], 1)
z = abs_normal_form.z  
s = abs_normal_form.num_switches
sigma_z = signature_vec(s,z)

# gradient formula interms of abs-linearization
function grad!(storage, x)
    c = vcat(alf_a', alf_b.* sigma_z)
    @. storage = c
end

# set bounds
o = Model(HiGHS.Optimizer)
MOI.set(o, MOI.Silent(), true)
@variable(o, lb_x[i] <= x[i=1:n] <= ub_x[i])

# abs-smooth lmo
lmo_as = AbsSmoothLMO(o, x_base, f, n, s, lb_x, ub_x)

# define termination criteria

# In case we want to stop the frank_wolfe algorithm prematurely after a certain condition is met,
# we can return a boolean stop criterion `false`.
# Here, we will implement a callback that terminates the algorithm if ||x_t+1 - x_t|| < eps.
function make_termination_callback(state)
 return function callback(state,args...)
  return norm(state.d) > 1e-3
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
    max_iteration=10
)

