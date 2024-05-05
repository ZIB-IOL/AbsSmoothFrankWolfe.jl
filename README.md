# AbsSmooth_FW.jl
Algorithm for nonsmooth optimization

## How to use this package? 
First you have to download the [ADOLC.jl](https://github.com/TimSiebert1/ADOLC.jl) package to use the AD toolbox.

Then you need the [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) package for the Conditional Gradient Descent methods.

## What this package does?
We are interested in optimization problems of the form 

$min\limits_{x\in C}$  $f(x)$ 

for convex compact $C$ and an [abs-smooth](https://optimization-online.org/wp-content/uploads/2012/09/3597.pdf) function $f$.

## Example
...

 
