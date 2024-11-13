# LASSO problem

Let's say we want to minimize the [LASSO](https://www.jstor.org/stable/2346178?seq=1) problem: $\frac{1}{2}\|Ax - y\|_2^2 + \rho \|x\|_1$, subjected to simple box constraints. 
This is what the code looks like:

```julia
julia> using AbsSmoothFrankWolfe,FrankWolfe,LinearAlgebra,JuMP,HiGHS

julia> import MathOptInterface

julia> const MOI = MathOptInterface

julia> n = 5 # choose lenght(x)

julia> p = 3 # choose lenght(y)

julia> rho = 0.5

julia> A = rand(p,n) # randomly choose matrix A

julia> y = rand(p) # randomly choose y

#define the LASSO function
julia>  function f(x)
	
 	return 0.5*(norm(A*x - y))^2 + rho*norm(x)

 end

# evaluation point x_base
julia> x_base = ones(n)*1.0

# box constraints
julia> lb_x = [-5 for in in x_base]

julia> ub_x = [5 for in in x_base]

# call the abs-linear form of f
julia> abs_normal_form = AbsSmoothFrankWolfe.abs_linear(x_base,f)

# gradient formula in terms of abs-linearization
julia> alf_a = abs_normal_form.Y

julia> alf_b = abs_normal_form.J 

julia> z = abs_normal_form.z 

julia> s = abs_normal_form.num_switches

julia> sigma_z = AbsSmoothFrankWolfe.signature_vec(s,z)

julia> function grad!(storage, x)
    c = vcat(alf_a', alf_b'.* sigma_z)
    @. storage = c
end

# define the model using JuMP with HiGHS as inner solver
julia> o = Model(HiGHS.Optimizer)

julia> MOI.set(o, MOI.Silent(), true)

julia> @variable(o, lb_x[i] <= x[i=1:n] <= ub_x[i])

# initialise dual gap
julia> dualgap_asfw = Inf

# abs-smooth lmo
julia> lmo_as = AbsSmoothFrankWolfe.AbsSmoothLMO(o, x_base, f, n, s, lb_x, ub_x, dualgap_asfw)

# define termination criteria using Frank-Wolfe 'callback' function
julia> function make_termination_callback(state)
 return function callback(state,args...)
  return state.lmo.dualgap_asfw[1] > 1e-2
 end
end

julia> callback = make_termination_callback(FrankWolfe.CallbackState)

# call abs-smooth-frank-wolfe
julia> x, v, primal, dual_gap, traj_data = AbsSmoothFrankWolfe.as_frank_wolfe(
    f, 
    grad!, 
    lmo_as, 
    x_base;
    gradient = ones(n+s),
    line_search = FrankWolfe.FixedStep(1.0),
    callback=callback,
    verbose=true
)

```

Beyond those presented in the documentation, more test problems can be found in the `examples` folder.

