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

# Highs/GLPK uses L(x,z,lambda,mu) = f(x,z) - lambda^T eq_const - mu^T ineq_const
# this is different to the derivation in the ASM
#  => use -lambda in normal growth condition!!

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
      if mu < min_mu_val
        index = i
        min_mu_val = mu
      end  
    end
  end
  return index
  
  end
  
function aasm(x_base, alpha, f_eval, n, s, ub_x, lb_x, outer_iter; max_inner_iter=10, model="model", mps=false) 

    # to store asm information
    lambdas = []
    solutions = []
    simplex_count = [] 
    
    iter = 1
    x_delta = copy(x_base)
    gap = Inf

    abs_normal_form = abs_linear(x_base,f_eval)

    alf_a = Float64.(abs_normal_form.Y)
    alf_b = Float64.(abs_normal_form.J) 
    z_base = Float64.(abs_normal_form.z)
    Z = Float64.(abs_normal_form.Z)  
    L = Float64.(abs_normal_form.L)
    y = Float64.(abs_normal_form.y)
    
    sigma_z = signature_vec(s, z_base)
    cz = Float64.(z_base - L * abs.(z_base))
    z_pl = copy(z_base)

    while iter <= max_inner_iter  
        println("----- inner iteration:  $iter / $max_inner_iter -----")

        o = Model(HiGHS.Optimizer)
        MOI.set(o, MOI.Silent(), true)
        @variable(o, xz[1:n+s])

        set_start_value.(xz[1:n], x_delta)
        set_start_value.(xz[n+1:end], z_pl)
        
        c_vec = Vector{Float64}(vec(vcat(alf_a', alf_b' .* sigma_z)))
        A_mat = Matrix{Float64}([Z L.*sigma_z'-I])
        cz_vec = Vector{Float64}(vec(cz))
        
        @objective(o, Min, sum(c_vec[i] * xz[i] for i in 1:n+s))

        lb_vals = max.(alpha .* (lb_x .- x_base), -1.0e30)
        ub_vals = min.(alpha .* (ub_x .- x_base),  1.0e30)

        for i in 1:n
            set_lower_bound(xz[i], lb_vals[i])
            set_upper_bound(xz[i], ub_vals[i])
        end

        # sigma_z[i] is always +1 or -1. This safely enforces xz >= 0 or xz <= 0
        for i in 1:s
            if sigma_z[i] > 0
                set_lower_bound(xz[n+i], 0.0)
            else
                set_upper_bound(xz[n+i], 0.0)
            end
        end
        
        c1 = @constraint(o, [i=1:s], sum(A_mat[i,j] * xz[j] for j in 1:n+s) == -cz_vec[i])
        
        optimize!(o)

        optimizer = backend(o)  
        info = MOI.get(optimizer, MOI.SimplexIterations())  

        push!(simplex_count, info)
        
        total_simplex_count = sum(simplex_count) 
          
        myxz = [value(var) for var in all_variables(o)]    
        myxz = round.(myxz, digits=5)
        push!(solutions, myxz)
        
        x_delta = myxz[1:n]
        z_pl = myxz[n+1:n+s]
            
        fpl_new = y + alf_a*x_delta + alf_b*(abs.(z_pl) - abs.(z_base))  
        
        if (fpl_new[1] > f_eval(x_base)[1] + 1.0e-12) || abs(fpl_new[1] - f_eval(x_base)[1]) < 1.0e-12
            gap = [0.0]        
        end

        mylambda = round.(dual.(c1), digits=5)
        push!(lambdas, mylambda)
        
        gap = (f_eval(x_base) .- fpl_new)/alpha
         
        index_mu = check_normalGrowth(s, alf_b, L, mylambda, z_pl) 
        
        if index_mu > -1
            sigma_z[index_mu] = -mylambda[index_mu]/abs(mylambda[index_mu])    
        else
            println("Local minimizer reached at iteration $iter. Exiting early.")
            total_inner_simplex = sum(simplex_count)
            return x_delta, gap, solutions, lambdas, iter, total_inner_simplex
        end
        
        iter = iter+1
    end
    total_inner_simplex = sum(simplex_count)
    return x_delta, gap, solutions, lambdas, iter, total_inner_simplex
end
