"""
 Custom LMO stores the abs-linearization of the function 
 compute_extreme_point() finds its minima at each iteration using AASM
"""
# abs-smooth custom lmo 

mutable struct AbsSmoothLMO <: FrankWolfe.LinearMinimizationOracle
    o
    x_base
    f
    n
    s
    lb_x 
    ub_x
    dualgap_asfw
    fw_iteration_counter::Int # for FW
    iteration_counter::Int # for inner iterations
end

AbsSmoothLMO(o, x_base, f, n, s, lb_x, ub_x, dualgap_asfw) = AbsSmoothLMO(o, x_base, f, n, s, lb_x, ub_x, dualgap_asfw, 0, 0) 

# abs-smooth custom extreme point -- AASM

function FrankWolfe.compute_extreme_point(lmo::AbsSmoothLMO, direction; kwargs...)
    f = lmo.f
    n = lmo.n
    s = lmo.s
    lb_x = lmo.lb_x
    ub_x = lmo.ub_x
    x_base = lmo.x_base
    inner_iteration = lmo.iteration_counter
    t = lmo.fw_iteration_counter
    alpha = 2 /(t+2)    
    x_delta, gap, _, _, iter = aasm(x_base, alpha, f, n, ub_x, lb_x, inner_iteration)
    
    # update abs-smooth dual gap
    lmo.dualgap_asfw = gap
    #@show lmo.dualgap_asfw
    
    # update inner_iteration
     lmo.fw_iteration_counter += 1
     lmo.iteration_counter += iter
     lmo.x_base = x_base + x_delta
    return x_base + x_delta        
end
