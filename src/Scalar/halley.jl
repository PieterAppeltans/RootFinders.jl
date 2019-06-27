######
# Halley's method for root finding
# Iteratively computes root of linear-over-linear Padé approximation to function
# Constraints on f? Twice continuous differentiable
# Convergence? Cubic
#
######


function halleys_method(x0::Number,f_df_ddf,accuracy, Kmax::Integer = 100)
    x = x0
    fx,dfx,ddfx = f_df_ddf(x)
    k = 1
    while (fx>accuracy && k < Kmax )
        x = x - (2*fx*dfx)/(2*dfx**2-fx*ddfx)
        fx,dfx = f_df(x)
        k+=1
    end
    return x
end
######
# If its cheaper to compute fx/dfx and dfx/ddfx then individually
#
######
function halleys_method(x0::Number,f_df_ddf,accuracy, Kmax::Integer = 100)
    x = x0
    fx_dfx,ddfx_dfx = f_df_ddf(x)
    k = 1
    while (fx>accuracy && k < Kmax )
        x = x - fx_dfx/(1-fx_dfx*ddfx_dfx/2)
        fx,dfx = f_df(x)
        k+=1
    end
    return x
end
