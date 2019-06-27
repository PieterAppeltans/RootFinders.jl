### dekker
# Calculate zero-crossing of f using Brent Dekker method in interval [lower,upper].
# Combines root bracketing (to avoid divergence), secant method and inverse quadratic interpolation.
# Also known as Wijngaarden-Brent-Dekker method
# Assumptions on f?
#   Continuous? yes, on [lower,upper]
#   Smooth?
# Convergence?
# Input:
#   - f: f(x) must return function value for x
#   - lower: lowerbound searchinterval
#   - upper: upperbound searchinterval
#   - tol: bound on bound on forward error |root-exact|<=tol
#   - Kmax: maximum number of function evaluations
# Output:
#   - root: a zero-crossing of f in interval [lower,upper]
#   - info:
#        0 - Successful
#       -1 - Lower bound is larger than upper bound
#       -2 - Search interval does not contain a root, root of even multiplicity, or multiple roots
#       -3 - Maximum number of function reached
# Reference:
#   - Dekker, T. J. (1969), "Finding a zero by means of successive linear interpolation", in Dejon, B.; Henrici, P., Constructive Aspects of the Fundamental Theorem of Algebra, London: Wiley-Interscience, ISBN 978-0-471-20300-1
###
function dekker(f,lower,upper,tol,Kmax)

    if (lower>upper)
        return nothing,-1
    end
    if (Kmax<2)
        return nothing,-3
    end
    flower = f(lower)
    fupper = f(upper)
    k = 2
    if (flower == 0)
        return lower,0
    elseif (fupper == 0)
        return upper, 0
    elseif (flower*fupper>0)
        # Function has same sign in both lower and upper boundary, f has possible
        # no root, multiple roots or a root of even multiplicity in the interval
        return nothing,-2
    end
    if (abs(flower)<abs(fupper))
        bk = lower; fbk = flower;
        ak = upper; fak = fupper;
    else
        bk = upper; fbk = fupper;
        ak = lower; fak = flower;
    end
    bk_1 = ak; fbk_1 = fak;
    while (abs(bk-ak)>tol)
        if (k>Kmax)
            return nothing,-3
        end
        # Calculate secant estimate
        estimate = bk - (bk-bk_1)/(fbk-fbk_1)
        # Calculate middle point interval
        middle = (ak+bk)/2
        bk_1 = bk; fbk_1 = fbk;
        if ((estimate-middle)*(bk-estimate)>=0)
            bk = estimate; fbk = f(estimate);
        else
            bk = middle; fbk = f(middle);
        end
        if (fbk*fak>0) #f(b_k+1) and f(a_k) have same sign
            ak = bk_1; fak = fbk_1;
        end
        if (abs(fak)<abs(fbk))
            ak,bk = bk,ak; fak,fbk = fbk,fak;
        end
        k = k+1
    end
    return bk
end

# Working principle
# Combines secant method with bisection method
