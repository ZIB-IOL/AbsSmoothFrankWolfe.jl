using ADOLC, ADOLC.array_types

function demo()
    # Define the Mifflin 2 function
    function f(x)
        return -x[1] + 2*(x[1]^2 + x[2]^2 - 1) + 1.75*abs(x[1]^2 + x[2]^2 - 1)
    end
    
    # Define the derivative evaluation point x
    x = [-1.0, -1.0]
    
    # Initialize the AbsNormalForm object
    abs_normal_form = ADOLC.AbsNormalForm()
    
    # Calculate the absolute normal form derivative
    derivative!(abs_normal_form, f, 1, 2, x, :abs_normal)
    @show abs_normal_form.Z
    @show abs_normal_form.Y
    
end

demo()
