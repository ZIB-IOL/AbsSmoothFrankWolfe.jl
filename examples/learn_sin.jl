using AbsSmoothFrankWolfe
using FrankWolfe
using LinearAlgebra
using JuMP
using HiGHS
using SparseArrays
using Random

using Plots

import MathOptInterface
const MOI = MathOptInterface

# Data
w = 4
xs = range(0, 2*pi, length=50)
X = reshape(collect(xs), :, 1)        
Y = reshape(sin.(w .* xs), :, 1)            

delta = step(xs)/2  
xs_test = xs .+ delta 

X_test = reshape(collect(xs_test), :, 1)
Y_test = reshape(sin.(w .* xs_test), :, 1)

input_size  = size(X, 1)   
output_size = size(Y, 1)   
hidden_size = 64         


# Weights & Biases
Random.seed!(1234)  
W1 = randn(hidden_size, input_size) * sqrt(2 / input_size)
b1 = zeros(hidden_size)                       
W2 = randn(output_size, hidden_size) * sqrt(2 / hidden_size)
b2 = zeros(output_size)                       

# Flatten lengths
const len_W1 = hidden_size * input_size
const len_b1 = hidden_size
const len_W2 = output_size * hidden_size
const len_b2 = output_size

# Flatten parameters
function vec_params(W1, b1, W2, b2)
    return vcat(vec(W1), b1, vec(W2), b2)
end

# Unpack parameters
function unpack_params(x)
    W1 = reshape(x[1:len_W1], hidden_size, input_size)
    b1 = x[len_W1+1 : len_W1+len_b1]

    W2_start = len_W1 + len_b1 + 1
    W2_end   = W2_start + len_W2 - 1
    W2 = reshape(x[W2_start:W2_end], output_size, hidden_size)

    b2_start = W2_end + 1
    b2_end   = b2_start + len_b2 - 1
    b2 = x[b2_start:b2_end]

    return W1, b1, W2, b2
end

# Activation
relu(x) = max.(0.0, x)


# Loss function (MSE)
function f(x)
    W1, b1, W2, b2 = unpack_params(x)
    a1 = relu.(W1 * X .+ b1)   
    y_pred = W2 * a1 .+ b2
    return sum((y_pred .- Y).^2)
end


# Set up ASFW
x_base = vec_params(W1, b1, W2, b2)
n = length(x_base)

lb_x = fill(-10.0, n)
ub_x = fill(10.0, n)

abs_normal_form = abs_linear(x_base,f)
alf_a = abs_normal_form.Y
alf_b = abs_normal_form.J 
z = abs_normal_form.z  
s = abs_normal_form.num_switches
sigma_z = signature_vec(s,z)

function grad!(storage, x)
    c = vcat(alf_a', alf_b'.* sigma_z)
    @. storage = c
end

o = Model(HiGHS.Optimizer)
MOI.set(o, MOI.Silent(), true)
@variable(o, lb_x[i] <= x[i=1:n] <= ub_x[i])

dualgap_asfw = Inf
lmo_as = AbsSmoothLMO(o, x_base, f, n, lb_x, ub_x, dualgap_asfw)

function make_termination_callback(state)
    return function callback(state,args...)
        return state.lmo.dualgap_asfw[1] > 1e0
    end
end

callback = make_termination_callback(FrankWolfe.CallbackState)

x, v, primal, dual_gap, traj_data = as_frank_wolfe(
    f, 
    grad!, 
    lmo_as, 
    x_base;
    gradient = ones(n+s),
    line_search = FrankWolfe.FixedStep(1.0),
    callback=callback,
    verbose=true,
    max_iteration=5e3
)


#--------TEST--------

function predict(x_vec, X_in)
    W1, b1, W2, b2 = unpack_params(x_vec)
    a1 = relu.(W1 * X_in .+ b1)         
    y_pred = W2 * a1 .+ b2             
    return y_pred
end

Y_pred = predict(x, X)           
Y_pred_test  = predict(x, X_test)      

train_loss = sum((Y_pred .- Y).^2) / length(Y)
test_loss  = sum((Y_pred_test .- Y_test).^2) / length(Y_test)

println("\nTraining MSE: ", train_loss)
println("Test MSE: ", test_loss)


xs_dense = range(0, 2*pi, length=200)  
plot(xs_dense, sin.(w .* xs_dense), label="True sin(w⋅x)", lw=1)
scatter!(xs, Y, label="Train points", ms=2, color=:blue)
scatter!(xs_test, vec(Y_pred_test), label="Predicted test", ms=2, color=:orange)
savefig("sin_w$(w).pdf")

