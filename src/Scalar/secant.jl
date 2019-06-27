### secant
# Constraints on f?
# Starting values must be "close" to root.
# Convergence: superlineair (convergence-order = (1+√5)/2 (golden ratio))
####
function secant(f,x0,x1,ftol,Kmax = 100,fx0=nothing,fx1=nothing)
    if (Kmax<2) # TODO include initial funcion evaluations
        return nothing,-3;
    end
    if (fx0 == nothing)
        fx0 = f(x0);
    end
    if (fx1 == nothing)
        fx1 = f(x1)
    end
    k = 2
    while (abs(fx1)>ftol)
        if (k>Kmax) # Maximum number of function evaluations reached
            return nothing,-3;
        end
        xnext = x1 - fx1*(x1-x0)/(fx1-fx0)
        if (isinf(xnext) || isnan(xnext))
            return nothing,-4;
        end
        x0,x1 = x1,xnext; fx0,fx1 = fx1,f(xnext);
        k += 1
    end
    return x1,0;
end

# Working principle
# Seek root of 1 order Taylor serie (linear) approximate of function where
# f'(x) is approximated using finite differences.
#
#  f(xk) + (f(xk)-f(xk_1))/(xk-xk_1)*(x-xk) = 0
#  xnext = xk - f(xk)*(xk-xk_1)/(f(xk)-f(xk_1))
