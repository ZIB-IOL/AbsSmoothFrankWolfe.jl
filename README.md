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
 using AbsSmoothFrankWolfe
using FrankWolfe
using LinearAlgebra
using JuMP
using HiGHS


import MathOptInterface
const MOI = MathOptInterface

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
    max_iteration=1e7
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

```

## Documentation
To explore the contents of this package, go to the [documentation](https://zib-iol.github.io/AbsSmoothFrankWolfe.jl/dev/).
Beyond those presented in the documentation, many more use cases are implemented in the `examples` folder.







