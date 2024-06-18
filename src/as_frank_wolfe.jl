"""
    frank_wolfe(f, grad!, lmo, x0; ...)

Simplest form of the Frank-Wolfe algorithm.
Returns a tuple `(x, v, primal, dual_gap, traj_data)` with:
- `x` final iterate
- `v` last vertex from the LMO
- `primal` primal value `f(x)`
- `dual_gap` final Frank-Wolfe gap
- `traj_data` vector of trajectory information.
"""
function as_frank_wolfe(
    f,
    grad!,
    lmo,
    x0;
    line_search::FrankWolfe.LineSearchMethod=FrankWolfe.Adaptive(),
    momentum=nothing,
    epsilon=1.0e-7,
    max_iteration=1.0e6,
    print_iter=1.0,
    trajectory=false,
    verbose=false,
    memory_mode::FrankWolfe.MemoryEmphasis=FrankWolfe.InplaceEmphasis(),
    gradient=nothing,
    callback=nothing,
    traj_data=[],
    timeout=3600,
    linesearch_workspace=nothing,
    dual_gap_compute_frequency=1.0,
)

    # header and format string for output of the algorithm
    headers = ["Type", "Iteration", "Primal", "||delta x||", "Primal gap", "Time", "It/sec"]
    format_string = "%6s %13s %14e %14e %14e %14e %14e\n"
    function format_state(state)
        rep = (
            FrankWolfe.st[Symbol(state.tt)],
            string(state.t),
            Float64(state.primal),
            Float64(norm(v-x0)),
            Float64(norm(f(v)-f(x0))),
            state.time,
            state.t / state.time,
        )
        return rep
    end

    t = 0
    dual_gap = Inf
    primal = Inf
    v = []
    x = x0
    tt = FrankWolfe.regular

    if trajectory
        callback = FrankWolfe.make_trajectory_callback(callback, traj_data)
    end

    if verbose
        callback = FrankWolfe.make_print_callback(callback, print_iter, headers, format_string, format_state)
    end

    time_start = time_ns()

    if (momentum !== nothing && line_search isa Union{FrankWolfe.Shortstep, FrankWolfe.Adaptive, FrankWolfe.Backtracking})
        @warn("Momentum-averaged gradients should usually be used with agnostic stepsize rules.",)
    end

    if verbose
        println("\nVanilla Frank-Wolfe Algorithm.")
        NumType = eltype(x0)
        println(
            "MEMORY_MODE: $memory_mode STEPSIZE: $line_search EPSILON: $epsilon MAXITERATION: $max_iteration TYPE: $NumType",
        )
        grad_type = typeof(gradient)
        println("MOMENTUM: $momentum GRADIENTTYPE: $grad_type")
        println("LMO: $(typeof(lmo))")
        if memory_mode isa FrankWolfe.InplaceEmphasis
            @info("In memory_mode memory iterates are written back into x0!")
        end
    end
    if memory_mode isa FrankWolfe.InplaceEmphasis && !isa(x, Union{Array,SparseArrays.AbstractSparseArray})
        # if integer, convert element type to most appropriate float
        if eltype(x) <: Integer
            x = copyto!(similar(x, float(eltype(x))), x)
        else
            x = copyto!(similar(x), x)
        end
    end

    # instanciating container for gradient
    if gradient === nothing
        gradient = collect(x)
    end

    first_iter = true
    if linesearch_workspace === nothing
        linesearch_workspace = FrankWolfe.build_linesearch_workspace(line_search, x, gradient)
    end

    # container for direction
    d = similar(x)
    gtemp = momentum === nothing ? d : similar(x)

    while t <= max_iteration 

        #####################
        # managing time and Ctrl-C
        #####################
        time_at_loop = time_ns()
        if t == 0
            time_start = time_at_loop
        end
        # time is measured at beginning of loop for consistency throughout all algorithms
        tot_time = (time_at_loop - time_start) / 1e9

        if timeout < Inf
            if tot_time ≥ timeout
                if verbose
                    @info "Time limit reached"
                end
                break
            end
        end

        #####################


        if momentum === nothing || first_iter
            grad!(gradient, x)
            if momentum !== nothing
                gtemp .= gradient
            end
        else
            grad!(gtemp, x)
            FrankWolfe.@memory_mode(memory_mode, gradient = (momentum * gradient) + (1 - momentum) * gtemp)
        end

        v = if first_iter
            compute_extreme_point(lmo, gradient)
        else
            compute_extreme_point(lmo, gradient, v=v)
        end   

        first_iter = false
        # go easy on runtime - only compute primal and dual if needed
        compute_iter = (
            (mod(t, print_iter) == 0 && verbose) ||
            callback !== nothing ||
            line_search isa FrankWolfe.Shortstep
        )
        if compute_iter
            primal = f(x)
        end
         
        d = FrankWolfe.muladd_memory_mode(memory_mode, d, x, v) 

        gamma = FrankWolfe.perform_line_search(
            line_search,
            t,
            f,
            grad!,
            gradient,
            x,
            d,
            1.0,
            linesearch_workspace,
            memory_mode,
        )
        t = t + 1
        if callback !== nothing
            state = FrankWolfe.CallbackState(
                t,
                primal,
                primal - dual_gap,
                dual_gap,
                tot_time,
                x,
                v,
                d,
                gamma,
                f,
                grad!,
                lmo,
                gradient,
                tt,
            )
            if callback(state) === false
                break
            end
        end  
               
        x = FrankWolfe.muladd_memory_mode(memory_mode, x, gamma, d)  
           
    end    

    return x, v, primal, dual_gap, traj_data, x0
end
