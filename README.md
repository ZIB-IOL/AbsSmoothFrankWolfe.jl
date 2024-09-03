[![DOI](https://zenodo.org/badge/793075266.svg)](https://zenodo.org/doi/10.5281/zenodo.11198550)

[![Build Status](https://github.com/shtadinada/AbsSmoothFW.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/shtadinada/AbsSmoothFW.jl/actions/workflows/CI.yml?query=branch%3Amain)

# AbsSmoothFW.jl
An algorithm for non-smooth optimization.

## Overview
We are interested in optimization problems of the form $\min\limits_{x\in C}$  $f(x)$ , for convex compact $C$ and an [abs-smooth](https://optimization-online.org/wp-content/uploads/2012/09/3597.pdf) function $f$.
We solve such problems by using [ADOLC.jl](https://github.com/TimSiebert1/ADOLC.jl/tree/master) for the AD toolbox and using the [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) for conditional gradient methods.

## Example
Look into examples folder to find 2 examples for small(n=2) and large(n=1000) test problems.

 
