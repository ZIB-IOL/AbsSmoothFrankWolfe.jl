include("../src/AbsSmooth_FW.jl")
using ADOLC

function abs_linearization()
    # Define the Mifflin 2 function
    function f(x)
        return -x[1] + 2*(x[1]^2 + x[2]^2 - 1) + 1.75*abs(x[1]^2 + x[2]^2 - 1)
    end
    
    # Define the derivative evaluation point x
    x_base = [-1.0,-1.0]
    n = length(x_base)
    
    # Initialize the AbsNormalForm object
    abs_normal_form = ADOLC.init_abs_normal_form(f, 1, n, x_base, tape_id=1)
    
    # Calculate the absolute normal form derivative
    derivative!(abs_normal_form, f, 1, n, x_base, :abs_normal, tape_id=abs_normal_form.tape_id, reuse_tape=true)
    
    println("AbsNormalForm at $x_base: ", abs_normal_form)  
    
    z = abs_normal_form.z
    Z = abs_normal_form.Z  
    L = abs_normal_form.L 
    alf_a = abs_normal_form.Y
    alf_b = reshape(abs_normal_form.J, size(abs_normal_form.J)[2], 1) 
    cz = abs_normal_form.cz 
    cy = abs_normal_form.cy 
    z = abs_normal_form.z  
    s = abs_normal_form.num_switches
    sigma_z = signature_vec(s,z)
    @show Z,L,alf_a,alf_b,cz,cy,z,s,sigma_z 
 return x_base,f,Z,L,alf_a,alf_b,cz,cy,z,s,sigma_z,n     
end

# bounds
lb_x = [-5, -5]
ub_x = [5, 5]

# call the abs-linear form of f
x_base,f,Z,L,alf_a,alf_b,cz,cy,z,s,sigma_z,n = abs_linearization()
@show x_base,f,Z,L,alf_a,alf_b,cz,cy,z,s,sigma_z,n

# gradient formula in terms of abs-linearization
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

