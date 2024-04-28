using ADOLC
using ADOLC.TbadoubleModule
using ADOLC.array_types

function call_adolc(x, f)
# call the abs_normal driver
n = length(x)
y = Vector{Float64}(undef, 1)
m = length(y)
a = [Adouble{TbAlloc}() for _ in 1:length(x)]
b = [Adouble{TbAlloc}() for _ in 1:length(y)]

tape_num = 0
trace_on(tape_num, 1)
a << x
b[1] = f(a)
b >> y
trace_off(0)

m = length(b)

abs_normal_problem = AbsNormalProblem{Float64}(tape_num, m, n, x, y)


abs_normal!(abs_normal_problem, tape_num)

    return abs_normal_problem
end
