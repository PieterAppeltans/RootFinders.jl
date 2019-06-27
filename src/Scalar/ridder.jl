### ridder
# Calculates (unique) root crossing of f in interval [lower,upper] using Ridder's method
# Constraints on f?
#
#
# Convergence: superlineair (convergence-order = √2)
# Reference:
#
###
function ridder(f,lower,upper,xtol,ftol,Kmax = 100)
    if (lower>upper)
        return nothing,-1
    end
    k = 0;
    flower = f(lower); fupper = f(upper);
    if (flower == 0)
        return lower,0
    elseif (fb == 0)
        return upper,0
    elseif (flower*fupper>0)
        return nothing,-1
    end
    while ((upper-lower)/2>tolerance)
        if (k>Kmax)
            return nothing,-3
        end
        beta = 0.5*(a+b)
        fbeta = f(beta)
        if (fbeta == 0)
            return beta
        end
        x = beta + 0.5*(b-a)*(sign(fa)*fbeta)/(sqrt(fbeta**2-fa*fb))
        fx = f(x)
        if (fx == 0)
            return x
        end
        if(fbeta*fx<0)
            if (sign(fa)>0)
                a = beta; b =x;
                fa = fbeta; fb = fx;
            else
                a = x; b = beta;
                fa = fx; fb = fbeta;
            end
        elseif (fa*fx<0)
            b = x; fb = fx;
        else
            a = x; fa = fx;
        end
        k += 1;
end
