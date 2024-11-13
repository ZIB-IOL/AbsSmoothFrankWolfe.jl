# AbsSmoothFrankWolfe.jl

[![CI](https://github.com/ZIB-IOL/AbsSmoothFW.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ZIB-IOL/AbsSmoothFW.jl/actions/workflows/CI.yml)
[![DOI](https://zenodo.org/badge/793075266.svg)](https://zenodo.org/doi/10.5281/zenodo.11198550)

This package is a toolbox for non-smooth version of the Frank-Wolfe algorithm. 

## Overview 
Abs-Smooth Frank-Wolfe algorithms are designed to solve optimization problems of the form $\min\limits_{x\in C}$  $f(x)$ , for convex compact $C$ and an [abs-smooth](https://optimization-online.org/wp-content/uploads/2012/09/3597.pdf) function $f$. 

We solve such problems by using [ADOLC.jl](https://github.com/TimSiebert1/ADOLC.jl/tree/master) for the AD toolbox and using the [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) for conditional gradient methods.


## Installation

The most recent release is available via the julia package manager, e.g., with

```julia
using Pkg
Pkg.add("AbsSmoothFrankWolfe")
```

or the main branch:

```julia
Pkg.add(url="https://github.com/ZIB-IOL/AbsSmoothFrankWolfe.jl", rev="main")
```

## Getting started

Let us consider the minimization of the abs-smooth function $max(x_1^4+x_2^2, (2-x_1)^2+(2-x_2)^2, 2*e^{(x_2-x_1)})$ subjected to simple box constraints $-5\leq x_i \leq 5$. Here is what the code will look like:

```julia
julia> using AbsSmoothFrankWolfe

julia> using FrankWolfe

julia> using LinearAlgebra

julia> using JuMP

julia> using HiGHS

julia> import MathOptInterface

julia> const MOI = MathOptInterface

julia> function f(x)
 	return max(x[1]^4+x[2]^2, (2-x[1])^2+(2-x[2])^2, 2*exp(x[2]-x[1]))
 end
 
# evaluation point x_base
julia> x_base = [3.0,2.0]

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

Vanilla Abs-Smooth Frank-Wolfe Algorithm.
MEMORY_MODE: FrankWolfe.InplaceEmphasis() STEPSIZE: FixedStep EPSILON: 1.0e-7 MAXITERATION: 1.0e7 TYPE: Float64
MOMENTUM: nothing GRADIENTTYPE: Vector{Float64}
LMO: AbsSmoothLMO
[ Info: In memory_mode memory iterates are written back into x0!

-------------------------------------------------------------------------------------------------
  Type     Iteration         Primal    ||delta x||       Dual gap           Time         It/sec
-------------------------------------------------------------------------------------------------
     I             1   8.500000e+01   2.376593e+00   1.281206e+02   0.000000e+00            Inf
  Last             7   2.000080e+00   2.000000e-05   3.600149e-04   2.885519e+00   2.425907e+00
-------------------------------------------------------------------------------------------------
x_final = [1.00002, 1.0]
```

## Documentation
To explore the contents of this package, go to the [documentation](https://zib-iol.github.io/AbsSmoothFrankWolfe.jl/dev/).
Beyond those presented in the documentation, many more use cases are implemented in the `examples` folder.







