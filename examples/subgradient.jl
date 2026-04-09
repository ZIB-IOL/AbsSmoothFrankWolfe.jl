using FrankWolfe
using LinearAlgebra
using JuMP
using HiGHS
using ADOLC


import MathOptInterface
const MOI = MathOptInterface

# Wong 2
 function f(x)
 	return max(x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45,
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(3*(x[1]-2)^2+4*(x[2]-3)^2+2*x[3]^2-7*x[4]-120),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(5*x[1]^2+8*x[2]+(x[3]-6)^2-2*x[4]-40),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(0.5*(x[1]-8)^2+2*(x[2]-4)^2+3*x[5]^2-x[6]-30),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(x[1]^2+2*(x[2]-2)^2-2*x[1]*x[2]+14*x[5]-6*x[6]),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(4*x[1]+5*x[2]-3*x[7]+9*x[8]-105),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(10*x[1]-8*x[2]-17*x[7]+2*x[8]),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(-3*x[1]+6*x[2]+12*(x[9]-8)^2-7*x[10]),
 		 x[1]^2+x[2]^2+x[1]*x[2]-14*x[1]-16*x[2]+(x[3]-10)^2+4*(x[4]-5)^2+(x[5]-3)^2+2*(x[6]-1)^2+5*x[7]^2+7*(x[8]-11)^2+2*(x[9]-10)^2+(x[10]-7)^2+45+10*(-8*x[1]+2*x[2]+5*x[9]-2*x[10]-12))
 end
 
 x_base = [2.0,3.0,5.0,5.0,1.0,2.0,7.0,3.0,6.0,10.0]

n = length(x_base)
 
# call the abs-linear form of f
const tape_id = 1
ADOLC.derivative(f, x_base, :jac, tape_id=tape_id)
c = CxxVector(x_base)

function grad!(storage, x)
 ADOLC.gradient!(c, f, n, x, tape_id, true)   # derivate in the direction y
 @. storage = c
end

o = HiGHS.Optimizer()
x = MOI.add_variables(o, n)

# set bounds
for i in 1:n
    MOI.add_constraint(o, x[i], MOI.Interval(-10.0, 10.0))
end

# abs-smooth lmo
lmo = FrankWolfe.MathOptLMO(o)


# call abs-smooth-frank-wolfe
x, v, primal, dual_gap, traj_data = FrankWolfe.frank_wolfe(
    f, 
    grad!, 
    lmo, 
    x_base;
    line_search = FrankWolfe.Agnostic(),
    verbose=true,
    print_iter=1,
    max_iteration=10000,
    epsilon = 0.1
)

#println("x_final = ", x_base)

