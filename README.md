# AbsSmoothFW.jl
Algorithm for nonsmooth optimization

Zenodo!!
## Pre-requisites
You need to export the ADOL-C libraries and the wrapper first by :
```julia
 ] add https://github.com/TimSiebert1/ADOLC_jll.jl.git#main
```
```julia
 ] add https://github.com/TimSiebert1/ADOLC.jl.git#master
```
## Overview
We are interested in optimization problems of the form $\min\limits_{x\in C}$  $f(x)$ , for convex compact $C$ and an [abs-smooth](https://optimization-online.org/wp-content/uploads/2012/09/3597.pdf) function $f$.

We use the AD toolbox from [ADOLC.jl](https://github.com/TimSiebert1/ADOLC.jl) package and [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) package for the Conditional Gradient Descent methods.

## Example
Look into examples folder to find 2 examples for small(n=2) and large(n=1000) test problems.

 
