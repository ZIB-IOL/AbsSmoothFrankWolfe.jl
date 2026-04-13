using AbsSmoothFrankWolfe
using FrankWolfe
using LinearAlgebra
using JuMP
using HiGHS
using SparseArrays
using Plots
using Statistics 

import MathOptInterface
const MOI = MathOptInterface

# Path to your local diabetes dataset
file_path = "dataset"  

# Initialize containers
A_rows = Vector{Vector{Float64}}()
y = Float64[]

# Read file line by line
for line in eachline(file_path)
    line = strip(line)
    # Skip lines starting with @ (like @DATA)
    if startswith(line, "@") || isempty(line)
        continue
    end
    # Split by comma
    values = split(line, ',')
    # Last column is target
    push!(y, parse(Float64, values[end]))
    # Features: first 10 columns
    push!(A_rows, parse.(Float64, values[1:end-1]))
end

# Convert A to proper matrix (rows × features)
A = hcat(A_rows...)' 

println("Feature matrix A size: ", size(A))  
println("Target vector y length: ", length(y)) 

p, n = size(A)

y_mean = mean(y)
A_mean = mean(A, dims=1)
A_std = std(A, dims=1)

y_centered = y .- y_mean
A_scaled = (A .- A_mean) ./ A_std # Z-score scaling

rho = 5
@show cond(A_scaled)


function f(x)
    return 0.5 * dot(A_scaled*x - y_centered, A_scaled*x - y_centered) + rho * sum(abs.(x))
end
 
x_base = zeros(n)
 
lb_x = fill(-20, n) 
ub_x = fill(500, n)

# call the abs-linear form of f
abs_normal_form = abs_linear(x_base,f)
s = abs_normal_form.num_switches

# gradient formula in terms of abs-linearization
function grad!(storage, x)
    c = ones(n+s)
    @. storage = c
end

# set bounds
o = Model(HiGHS.Optimizer)
MOI.set(o, MOI.Silent(), true)
@variable(o, lb_x[i] <= x[i=1:n] <= ub_x[i])

# initialise dual gap 
dualgap_asfw = Inf

# abs-smooth lmo
lmo_as = AbsSmoothLMO(o, x_base, f, n, s, lb_x, ub_x, dualgap_asfw)

# define termination criteria
function make_termination_callback(state)
    return function callback(state,args...)
        return state.lmo.dualgap_asfw[1] > 1e3
    end
end

callback = make_termination_callback(FrankWolfe.CallbackState)

x_final, v, primal, dual_gap, traj_data, x0 = as_frank_wolfe(
    f, 
    grad!, 
    lmo_as, 
    x_base;
    gradient = ones(n+s),
    line_search =  FrankWolfe.FixedStep(1.0),
    callback=callback,
    verbose=true,
    max_iteration=45000
)


x_original = x_final ./ vec(A_std)

intercept = y_mean - dot(vec(A_mean), x_original)
println("Intercept = ", intercept)

y_hat = intercept .+ A * x_original
mse = sum((y .- y_hat).^2) / length(y)
println("Mean Squared Error: ", mse)

scatter(
    y, y_hat,
    xlabel = "Actual Disease Progression (y)",
    ylabel = "Predicted Disease Progression (ŷ)",
    title = "Predicted vs. Actual Responses",
    label = "Patients",
    marker = (:circle, 8),
    legend = :topright
)

println("Grand Total Simplex Steps: ", lmo_as.total_simplex_count)
println("x_final = ", x_base)

# Add diagonal line for reference (perfect prediction)
plot!([minimum(y), maximum(y)], [minimum(y), maximum(y)], 
      line=:dash, color=:red, label="Perfect Prediction")

savefig("predictions.pdf")
