###
# brent
#
# Reference
#   Brent, R. P. (1973), "Chapter 4: An Algorithm with Guaranteed Convergence for Finding a Zero of a Function", Algorithms for Minimization without Derivatives, Englewood Cliffs, NJ: Prentice-Hall, ISBN 0-13-022335-2
###

# use inverse quadratic interpolation from __helper.jl
# TODO: tolerance interpolation
# http://mathfaculty.fullerton.edu/mathews/n2003/BrentMethodMod.html
# numpy
function brent(f,lower,upper,tol,Kmax)
    k = 0
    if (lower > upper)
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
    ck = ak; fck = fck;
    ck_1 = ck; fck_1 = fck;
    interpolation_k_1 = false
    while (abs(bk-ak)>tol)
        if (k>Kmax)
            return nothing,-3
        end
        if (abs(fak-fck)>1e-5 && abs(fbk-fck)>1e-5)
            estimate = (fbk*fck)/((fak-fbk)*(fak-fck)) +
                        (fak*fck)/((fbk-fak)*(fbk-fck)) +
                        (fak*fbk)/((fck-fak)*(fck-fbk))
        else
            # Calculate secant estimate
            estimate = bk - (bk-bk_1)/(fbk-fbk_1)
        end
        # Calculate middle point interval
        middle = (ak+bk)/2
        intermediate = (3*ak+bk)/4

        if (estimate-intermediate)*(bk-estimate)>=0 &&
            (interpolation_k_1 || (abs(bk-ck)>tol_interpolation && abs(estimate-bk)<0.5*abs(bk-ck))) &&
            (!interpolation_k_1 || (abs(ck-ck_1)>tol_interpolation && abs(estimate-bk)<0.5*abs(ck-ck_1) ))

            interpolation_k_1 = true;
            s = estimate; fs = f(estimate);
        else
            interpolation_k_1 = false;
            s = middle; fs = f(middle);
        end
        ck_1,ck = ck,bk;fck_1,fck = fck,fbk;
        bk = s; fbk = fs;
        if (fs*fak>0) #f(b_k+1) and f(a_k) have same sign
            ak = ck; fak = fck;
        end
        if (abs(fak)<abs(fbk))
            ak,bk = bk,ak; fak,fbk = fbk,fak;
        end
        k = k+1
    end
    return bk
end
# Working principle: https://en.wikipedia.org/wiki/Brent%27s_method
# Tries to use inverse quadratic interpolation
# (instead of secant method in Dekker's method) if possible.
# If not possible uses Secant method
# If estimate outside bracket or following conditions not met
#  -> bisection method
# inverse quadratic Lagrange interpolation of degree 2
#
