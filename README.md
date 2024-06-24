[![DOI](https://zenodo.org/badge/793075266.svg)](https://zenodo.org/doi/10.5281/zenodo.11198550)

>[!IMPORTANT]
>ADOLC.jl has been changed, making changes here.
# AbsSmoothFW.jl
An algorithm for non-smooth optimization.

## Pre-requisites
You need to export the ADOL-C libraries and the wrapper first by :
```julia
 ] add https://github.com/TimSiebert1/ADOLC_jll.jl.git#main
```
```julia
 ] add https://github.com/TimSiebert1/ADOLC.jl.git#master
```
more information can be found at [ADOLC.jl](https://github.com/TimSiebert1/ADOLC.jl).
## Overview
We are interested in optimization problems of the form $\min\limits_{x\in C}$  $f(x)$ , for convex compact $C$ and an [abs-smooth](https://optimization-online.org/wp-content/uploads/2012/09/3597.pdf) function $f$.

## Example
Look into examples folder to find 2 examples for small(n=2) and large(n=1000) test problems.

 
