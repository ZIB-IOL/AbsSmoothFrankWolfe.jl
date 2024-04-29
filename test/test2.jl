using ADOLC, ADOLC.array_types

function demo()
 function f(x)
        return (max(-x[1]-x[2], -x[1]-x[2]+x[1]^2+x[2]^2-1) + max(-x[2]-x[3], -x[2]-x[3]+x[2]^2+x[3]^2-1))
    end 
    
    x = [-1.5, -1.5, -1.5]

    abs_normal_problem = ADOLC.AbsNormalForm()
    derivative!(abs_normal_problem, f, 1, 3, x, :abs_normal)
    y = f(x)

    @show abs_normal_problem.y[1]

    x = [-0.5, -0.5, -0.5]
    # reuse abs_normal_problem with same id and without retaping
    derivative!(abs_normal_problem, f, 1, 3, x, :abs_normal, tape_id=abs_normal_problem.tape_id, reuse_tape=true)
    y = f(x)

     @show abs_normal_problem.Y[1, 1]
end
 demo()

