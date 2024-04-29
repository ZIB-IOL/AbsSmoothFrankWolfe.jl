using ADOLC

function call_adolc(x, f)
# call the abs_normal driver
n = length(x)
y = Vector{Float64}(undef, 1)
m = length(y)


# Initialize the AbsNormalForm object
 abs_normal_form = ADOLC.AbsNormalForm()

# Calculate the absolute normal form derivative
 derivative!(abs_normal_form, f, m, n, x, :abs_normal)


    return abs_normal_form
end

function call_adolc_reuse(x, f)
# call the abs_normal driver
n = length(x)
y = Vector{Float64}(undef, 1)
m = length(y)


# Initialize the AbsNormalForm object
 abs_normal_form = ADOLC.AbsNormalForm()

# Calculate the absolute normal form derivative
 derivative!(abs_normal_form, f, m, n, x, :abs_normal, tape_id = abs_normal_form.tape_id, reuse_tape =true)


    return abs_normal_form
end
