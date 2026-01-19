"""
Adapted Active Signature Method - aasm( x, step_size, f, iterations;..)
 
x_base: point at which the piecewise local model is build
f: abs-smooth function
lb_x, ub_x: lower and upper bound on x
alpha: step size of FW method
remaining arguments: function pointer
"""
###################################################################   

# get signature vector up to accuracy of myeps
"""
Signature vector at the point 'x' upto certain accuracy
"""
function signature_vec(s, z; myeps = 1.0e-12)         
  sigma_z=sign.(z)                                  
  
  for i=1:s
    if abs(z[i])<myeps
      sigma_z[i]=1.0;
    end
  end

  return sigma_z
end

# check optimality

###############################################################################
# Highs/GLPK uses L(x,z,lambda,mu) = f(x,z) - lambda^T eq_const - mu^T ineq_const
# this is different to the derivation in the ASM
#  => use -lambda in normal growth condition!!
###############################################################################

function check_normalGrowth(s, b, L, lambda, z; myeps=1.0e-12)   

  index = -1
  
  min_mu_val = 1.0/0.0
  
  for i=1:s
    mu = b[i]
    for j=1:s
      mu = mu-lambda[j]*L[j,i]                      
    end
    mu = mu-abs(lambda[i])
 
      
    if mu<0 && abs(z[i])<= myeps && abs(lambda[i])>=myeps
  #  @show mu
      if mu < min_mu_val
        index = i
        min_mu_val = mu
        #  
      end  
    end
  end
  return index
  
  end

###############################################################################
simplex_count = []
lp_count = []

# adapted active signature method (AASM)
function aasm(x_base, alpha, f_eval, n, s, ub_x, lb_x, outer_iter; max_inner_iter=100, model="model", mps=false) 

# to store asm information
 lambdas = []
 solutions = []
   
 iter = 1
 x_delta = copy(x_base)
 gap = Inf
#     println("--------0----------")   
# call the abs-linear form of f
 abs_normal_form = abs_linear(x_base,f_eval)
#   println("--------1----------")   
# s = copy(abs_normal_form.num_switches)
#    println("--------2----------")   
 alf_a = copy(abs_normal_form.Y)
 alf_b = copy(abs_normal_form.J) 
#    println("--------3----------")   
 z_base = copy(abs_normal_form.z)
 sigma_z = signature_vec(s,z_base)
#    println("--------4----------")   
  Z = copy(abs_normal_form.Z)  
  L = copy(abs_normal_form.L)
#    println("--------5----------")  
    #@show Z,L 
 #cz = abs_normal_form.cz
 #cy = abs_normal_form.cy

 y = copy(abs_normal_form.y)
#    println("--------6----------")   
 cz = z_base - L * abs.(z_base)
 z_pl = copy(z_base)
 #@show z_pl,z_base
#  println("--------7----------")   
 while iter <= max_inner_iter  
	println("----- inner iteration:  $iter / $max_inner_iter -----")
	o = Model(HiGHS.Optimizer)
	MOI.set(o, MOI.Silent(), true)
#   println("--------8----------")   
	# variables: xz = (x,z)  => n+s variables
	@variable(o, xz[1:n+s])
   
	set_start_value.(xz[1:n], x_delta)
	set_start_value.(xz[n+1:end], z_pl)

	c = vcat(alf_a', alf_b' .* sigma_z)      
	@objective(o, Min, dot(c, xz))
#  println("--------9----------")   
################## define bound constraints ############################## 
# bounds on x from original problem - linear shift due to linearization
# here we compute \Delta v hence
# v \in C => lb_v <= v <= ub_v
# lb_v-x_base <= v-x_base <= ub_v-x_base
# alpha(lb_v-x_base) <= alpha(v-x_base) <= alpha(ub_v-x_base)
# alpha(lb_v-x_base) <= x_delta <= alpha(ub_v-x_base)
##########################################################################
    
	@constraint(o, xz[1:n] <= min.(alpha*(ub_x-x_base),1.0e30))          
	@constraint(o, xz[1:n] >= max.(alpha*(lb_x-x_base),-1.0e30))         
#	ub = alpha .* (ub_x .- x_base)
#	lb = alpha .* (lb_x .- x_base)

#	ub .= ifelse.(isfinite.(ub), ub, 1.0e30)
#	lb .= ifelse.(isfinite.(lb), lb, -1.0e30)

#	@constraint(o, xz[1:n] .<= ub)
#	@constraint(o, xz[1:n] .>= lb)

	@constraint(o, sigma_z .* xz[n+1:end] .>= 0.0)  
    
	A = [Z L.*sigma_z'-I]
#    println("--------10----------")
    #@show A   
	c1 = @constraint(o, A*xz .== -cz)
#    println("--------11----------")   
	optimize!(o)
#    println("--------11.1----------")

	# Retrieve the solver information to get the number of simplex iterations
	# Get the HiGHS backend
	optimizer = backend(o)  
	# Retrieve simplex iterations
	info = MOI.get(optimizer, MOI.SimplexIterations())  
#	num_lp_iters = MOI.get(optimizer, MOI.RawOptimizerAttribute("num_lp_solve_calls"))
#	info = HiGHS.get_info(o)
#	info.num_lp_iterations
      
	push!(simplex_count, info)
#	push!(lp_count, num_lp_iters)
	#@show simplex_count
     
	total_simplex_count = 0
#	total_lp_count = 0
     
	for value in simplex_count
	total_simplex_count += value
#	total_lp_count += value
	end
     
	@show total_simplex_count #total_lp_count
#	@dbg "total_simplex_count = ", total_simplex_count

#        println("simplex: ", total_simplex_count)
#        flush(stdout)
         
	myxz = [value(var) for var in all_variables(o)]    
	myxz = round.(myxz, digits=5)
	push!(solutions, myxz)
#     println("--------12----------")   
	x_delta = myxz[1:n]
	#@show x_delta + x_base
#	println("-------12.1-----------")   
	z_pl = myxz[n+1:n+s]
	#@show count(abs.(z_pl) .< 1.0e-10), length(z_pl)
	
#     	println("--------12.2----------")   
	# new abs linearization  
	fpl_new = y + alf_a*x_delta + alf_b*(abs.(z_pl) - abs.(z_base))  
	
#	println("--------13----------")    

	if (fpl_new[1] > f_eval(x_base)[1] + 1.0e-12) || abs(fpl_new[1] - f_eval(x_base)[1]) < 1.0e-12
	gap = [0.0]  #fpl_new = f_eval(x_base)		
	end

	mylambda = round.(dual.(c1), digits=3)
	push!(lambdas, mylambda)
#    println("--------14----------")   
	# dual-gap for abs-smooth case
	gap = (f_eval(x_base) .- fpl_new)/alpha
       
	# local optimality criterion
	index_mu = check_normalGrowth(s, alf_b, L, mylambda, z_pl) 
#	println("--------15----------")   
	if index_mu > -1
	# current point not optimal, change polyhedron 
	sigma_z[index_mu]=-mylambda[index_mu]/abs(mylambda[index_mu]);   
	else
	# local minimizer reached     
	return x_delta, gap, solutions, lambdas, iter
	end
    
	iter = iter+1
 end

 return x_delta, gap, solutions, lambdas, iter 
end

