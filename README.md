# AbsSmoothFW.jl
more features to come soon..

## Overview
We are interested in optimization problems of the form $\min\limits_{x\in C}$  $f(x)$ , for convex compact $C$ and an [abs-smooth](https://optimization-online.org/wp-content/uploads/2012/09/3597.pdf) function $f$.
Using the AD toolbox from [ADOLC.jl](https://github.com/TimSiebert1/ADOLC.jl) package and [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) package for the Conditional Gradient Descent methods, we device an algortihm for nonsmooth optimization.

## Example
Look into examples folder to find 2 examples for small(n=2) and large(n=1000) test problems.

 
