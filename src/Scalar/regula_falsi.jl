### regula_falsi
# Calculates zero-crossing of f in the interval [lower,upper] using regula falsi
#
#
#
#
function regula_falsi(f,lower,upper,xtol,ftol=0,Kmax=100,flower,fupper)
    if (lower > upper)
        return nothing,-1
    end
    flower = f(lower);
    fupper = f(upper);
    k = 0
    while (upper-lower>xtol)
        if (k>Kmax)
            return nothing,-4
        end
        xnext = lower - (upper-lower)/(fupper-flower)*flower
        fnext = f(xnext)
        if (fnext*flower>0)
            lower,flower = xnext,fnext
        else
            upper,fupper = xnext,fupper
        end
        k += 1
    end
    return lower - (upper-lower)/(fupper-flower)*flower,0
end

# Working principle
# Similar to bisection method. But next iterate is x intercept of line connecting
# f(lower) and f(upper) instead of middle
# f(lower) + (x-lower)*(f(upper)-f(lower))/(upper-lower) = 0
# xnext = lower-(upper-lower)/(f(upper)-f(lower))*f(lower)
