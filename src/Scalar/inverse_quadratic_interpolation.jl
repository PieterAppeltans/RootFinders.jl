### inverse_quadratic_interpolation
# Calculates the zero crossing of the function f using the inverse quadratic interpolation method
#
# Reference:
#   https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
function inverse_quadratic_interpolation(f,xk_2,xk_1,xk,tol,Kmax)
    k = 0
    fxk_2 = f(xk_2); fxk_1 = f(xk_1); fxk = f(xk);
    while (abs(fxk_2)<tol)
        if (k>Kmax)
            return nothing, -3
        end
        xk_new = fxk*fxk_1/((fxk_2-fxk_1)*(fxk_2-fxk))*xk_2 +
            fxk_2*fxk/((fxk_1-fxk_2)*(fxk_1-fxk))*xk_1 +
            fxk_2*fxk_1/((fxk-fxk_2)*(fxk-fxk_1))*xk
        if (isinf(xk_new))
            return nothing, -4
        end
        xk_2,xk_1,xk = xk_1,xk,xk_new
        fxk_2,fxk_1,fxk = fxk_1,fxk, f(xk)
    end
end

# Working principle
