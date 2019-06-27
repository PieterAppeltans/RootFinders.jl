### muller
# Calculates a zero-crossing of f(x) using Muller's method
# Constraints on f?
#   - real-valued
#   -
# Input:
#   - f f(x) (!x may be complex!)
# Reference:
#  Muller, David E., "A Method for Solving Algebraic Equations Using an Automatic Computer," Mathematical Tables and Other Aids to Computation, 10 (1956), 208-215.
#  https://en.wikipedia.org/wiki/Muller%27s_method

# TODO reuse divided differences
function muller(f,xk::Complex{T0},xk_1::Complex{T1},xk_2::Complex{T2},ftol,Kmax=100) where T0,T1,T2

    fxk = f(xk); fxk_1 = f(xk_1); fxk_2 = f(xk_2);
    while (abs(fxk)>ftol)
        w = (fxk_1-fxk)/(xk_1-xk)+(fxk_2-fxk)/(xk_2-xk)-
            (fxk_2-fxk_1)/(xk_2-xk_1)
        div_dif = 1/(x2-x0)*((fxk_2-fxk_1)/(xk_2-xk_1)-(fxk_1-fxk)/(xk_1-xk))
        t = sqrt(w**2-4*fxk*div_dif)
        den = max(abs2(w+t),abs2(w-t))
        xnext = xk - (2*fxk)/den
        if (isinf(xnext) || isnan(xnext))
            return nothing,-4
        end
        xk_2,xk_1,xk = xk_1,xk,xnext;
        fxk_2,fxk_1,fxk = fxk_1,fxk,f(xk)
    end
    return fxk,0
end
function muller(f,xk<:Real,xk_1<:Real,xk_2<:Real,ftol,Kmax=100)
    muller(f,comlex(xk),complex(xk_1),complex(xk_2),ftol,Kmax=100)
end
