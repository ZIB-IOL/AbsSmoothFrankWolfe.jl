# Adapted Active Signature Method 

###################################################################   

using LinearAlgebra
using SCIP, HiGHS, GLPK
import MathOptInterface
const MOI = MathOptInterface
using JuMP

# get signature vector up to accuracy of myeps

function signature_vec(s, z; myeps = 1.0e-12)         
  sigma_z=sign.(z)                                  
  
  for i=1:s
    if abs(z[i])<myeps
      sigma_z[i]=0;
    end
  end

  return sigma_z
end

# check optimality

###############################################################################
# Highs/GLPK uses L(x,z,lambda,mu) = f(x,z)-lambda^T eq_const - mu^T ineq_const
# this is different to the derivation in the ASM
#  => use -lambda in normal growth condition!!
###############################################################################

function check_normalGrowth(s, b, L, sigma_z, lambda, z; myeps=1.0e-12)   

  index = -1
  
  min_mu_val = 1.0/0.0
  
  for i=1:s
    mu = b[i]
    for j=1:s
      mu = mu-lambda[j]*L[i,j]                      
    end
    mu = mu-abs(lambda[i])
      
    if mu<0 && abs(sigma_z[i]==0) && abs(lambda[i])>=myeps
      if mu < min_mu_val
        index = i
        min_mu_val = mu
          
      end  
    end
  end
  return index
  
  end

###############################################################################

# adapted active signature method (AASM)

function aasm(x_base, alpha, f_eval, outer_iter; max_inner_iter=100, model="model", mps=false) 

#  x_base: point at which the piecewise local model is build
#  lb_x, ub_x: lower and upper bound on x
#  alpha: step size of FW method
#  remaining arguments: function pointer
   
# abs-linearization of f
  abs_normal_form = call_adolc(x_base, f_eval) 
  
  z = abs_normal_form.z
  Z = abs_normal_form.Z  
  L = abs_normal_form.L 
  alf_a = abs_normal_form.Y
  alf_b = reshape(abs_normal_form.J, size(abs_normal_form.J)[2], 1) 
  cz = abs_normal_form.cz 
  cy = abs_normal_form.cy 
  
# to store asm information
  lambdas = []
  solutions = []
  
  iter = 1
  while iter <= max_inner_iter  
    
    o = Model(HiGHS.Optimizer)
    MOI.set(o, MOI.Silent(), true)

# variables: xz = (x,z)  => n+s variables
    @variable(o, xz[1:n+s])
    
    set_start_value.(xz[1:n], x_base)
    set_start_value.(xz[n+1:end], z)
    
    sigma_z = signature_vec(s,z)
    
# define objective
    c = vcat(alf_a', alf_b .* sigma_z)       
    @objective(o, Min, dot(c, xz))
  
################## define bound constraints ############################## 
# bounds on x from original problem - linear shift due to linearization
# here we compute \Delta v hence
# v \in C => lb_v <= v <= ub_v
# lb_v - x_base <= v - x_base <= ub_v-x_base
# alpha(lb_v - x_base) <= alpha(v - x_base) <= alpha(ub_v-x_base)
# alpha(lb_v - x_base) <= x_delta <= alpha(ub_v-x_base)
##########################################################################
    
    @constraint(o, xz[1:n] <= min.(alpha*(ub_x-x_base),1.0e30))          
    @constraint(o, xz[1:n] >= max.(alpha*(lb_x-x_base),-1.0e30))         
  
    @constraint(o, sigma_z .* xz[n+1:end] .>= 0) 
    for i in eachindex(abs_normal_form.L)
     abs_normal_form.L[i]
    end  
    
    for i in eachindex(abs_normal_form.Z)
     abs_normal_form.Z[i]
    end      
    x_base
    A = [Z L.*sigma_z'-I]
    
    c1 = @constraint(o, A*xz .== -cz)
    optimize!(o)
    
    myxz = [value(var) for var in all_variables(o)]    
    myxz = round.(myxz, digits=5)
    push!(solutions, myxz)
  
    x_delta = myxz[1:n]
    z = myxz[n+1:n+s]
     
# new abs linearization
    fpl_new = cy + alf_a*x_delta + alf_b'*abs.(z) 
    
# dual-gap for abs-smooth case
   fabs = alf_a*x_delta + alf_b'*abs.(z)
 
    mylambda = round.(dual.(c1), digits=3)
    push!(lambdas, mylambda)
    
# local optimality criterion
    index_mu = check_normalGrowth(s, alf_b, L, sigma_z, mylambda, z) 
    
    if index_mu > -1
 # current point not optimal, change polyhedron 
      sigma_z[index_mu]=-mylambda[index_mu]/abs(mylambda[index_mu]);       
    else
 # local minimizer reached      
     return x_delta, fabs, solutions, lambdas, iter
    end
    iter = iter+1
  end

  return x_delta, fabs, solutions, lambdas, iter 
end
