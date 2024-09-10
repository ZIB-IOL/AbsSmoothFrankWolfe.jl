# AbsSmoothFrankWolfe.jl

[![CI](https://github.com/ZIB-IOL/AbsSmoothFW.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ZIB-IOL/AbsSmoothFW.jl/actions/workflows/CI.yml)
[![DOI](https://zenodo.org/badge/793075266.svg)](https://zenodo.org/doi/10.5281/zenodo.11198550)

This package is a toolbox for Abs-Smooth Frank-Wolfe algorithm.

## Overview

Abs-Smooth Frank-Wolfe algorithms are designed to solve optimization problems of the form $\min\limits_{x\in C}$  $f(x)$ , for convex compact $C$ and an [abs-smooth](https://optimization-online.org/wp-content/uploads/2012/09/3597.pdf) function $f$.


## Installation

The most recent release is available on the main branch:

```julia
Pkg.add(url="https://github.com/ZIB-IOL/AbsSmoothFW.jl")
```

## Getting started

Let's say we want to minimize the [LASSO](https://www.jstor.org/stable/2346178?seq=1) problem: $\frac{1}{2}\|Ax - y\|_2^2 + \rho \|x\|_1$, subjected to simple box constraints. 
This is what the code looks like:

```julia
julia> using FrankWolfe,SparseArrays,LinearAlgebra,JuMP,HiGHS,ADOLC

julia> import MathOptInterface

julia> const MOI = MathOptInterface

julia> include("../src/aasm.jl")

julia> include("../src/as_frank_wolfe.jl")

julia> include("../src/abs_linear.jl")

julia> include("../src/abs_lmo.jl")

julia> n = 5 # lenght(x)

julia> p = 3 # lenght(y)

julia> rho = 0.5

julia> A = rand(p,n)

julia> y = rand(p)

julia>  function f(x)
	
 	return 0.5*(norm(A*x - y))^2 + rho*norm(x)

 end

# evaluation point x_base
julia> x_base = ones(n)*1.0

# box constraints
julia> lb_x = [-5 for in in x_base]

julia> ub_x = [5 for in in x_base]

# call the abs-linear form of f
julia> abs_normal_form = abs_linear(x_base,f)

# gradient formula in terms of abs-linearization
julia> alf_a = abs_normal_form.Y

julia> alf_b = abs_normal_form.J 

julia> z = abs_normal_form.z 

julia> s = abs_normal_form.num_switches

julia> sigma_z = signature_vec(s,z)

julia> function grad!(storage, x)
    c = vcat(alf_a', alf_b'.* sigma_z)
    @. storage = c
end

# define the model
julia> o = Model(HiGHS.Optimizer)

julia> MOI.set(o, MOI.Silent(), true)

julia> @variable(o, lb_x[i] <= x[i=1:n] <= ub_x[i])

# initialise dual gap
julia> dualgap_asfw = Inf

# abs-smooth lmo
julia> lmo_as = AbsSmoothLMO(o, x_base, f, n, s, lb_x, ub_x, dualgap_asfw)

# define termination criteria
julia> function make_termination_callback(state)
 return function callback(state,args...)
  return state.lmo.dualgap_asfw[1] > 1e-2
 end
end

julia> callback = make_termination_callback(FrankWolfe.CallbackState)

# call abs-smooth-frank-wolfe
julia> x, v, primal, dual_gap, traj_data = as_frank_wolfe(
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

