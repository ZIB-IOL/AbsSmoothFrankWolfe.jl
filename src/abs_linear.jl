"""
  one can get the abs-linearization of given f by abs_linear(x,f)
x: base point
f: abs-smooth func
"""

function abs_linear(x,f)
    n = length(x)

    # Initialize the AbsNormalForm object
    abs_normal_form = ADOLC.init_abs_normal_form(f, 1, n, x, tape_id=1)
    
    # Calculate the absolute normal form derivative
    derivative!(abs_normal_form, f, 1, n, x, :abs_normal, tape_id=abs_normal_form.tape_id, reuse_tape=true)
       
    return abs_normal_form
end
