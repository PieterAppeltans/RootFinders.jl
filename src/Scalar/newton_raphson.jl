function newton_raphson(x0::Number,f,df,accuracy,Kmax::Integer = 100)


end
###
# Constraints on f? differentiable
#
# Convergence: quadratic
###
function newton_raphson(x0::Number,f_df,accuracy, Kmax::Integer = 100)
    x = x0
    fx,dfx = f_df(x)
    k = 1
    while (fx>accuracy && k < Kmax )
        x = x - fx/dfx
        fx,dfx = f_df(x)
        k+=1
    end
    return x
end
function newton_raphson_(x0::Number,f,df,accuracy,Kmax::Integer = 100)

end
function newton_raphson_(x0::Number,f_df,accuracy, Kmax::Integer = 100)
    x = x0
    fx,dfx = f_df(x)

    k = 1
    x_vec = Vector{AbstractFloat}(nothing,Kmax)
    fx_vec = Vector{AbstractFloat}(nothing,Kmax)
    dfx_vec = Vector{AbstractFloat}(nothing,Kmax)
    x_vec[k] = x; fx_vec[k] = fx; dfx_vec[k] = dfx;
    while (fx>accuracy && k < Kmax )
        x = x - fx/dfx
        fx,dfx = f_df(x)
        k+=1
        x_vec[k] = x; fx_vec[k] = fx; dfx_vec[k] = dfx;
    end
    return newton_raphson_result(x,x_vec,fx_vec,dfx_vex)
end

struct newton_raphson_result
    x::Number
    x_vec::AbstractVector{Number}
    fx_vec::AbstractVector{Number}
    dfx_vec::AbstractVector{Number}
end
