###
# bisection
# Calculate (unique) zero-crossing of scalar function f in interval [lower,upper].
# Assumptions on f?
#   scalar
#   continuous: no
#   smooth: no, derivatives may not exist or be continuous
#   other: searchinterval only contains one root of uneven multiplicity
# Convergence:
#   Each iterationstep search interval halfs.
#   Number of function evaluations: 1+ceil(log_2(|upper-lower|/tol))
# Input:
#   - f: f(x) must return function value for x
#   - lower: lowerbound searchinterval
#   - upper: upperbound searchinterval
#   - xtol: bound on forward error |root-exact|<=tol
#   - ftol: bound on backward error |f(root)| <= tol (TODO)
# Output:
#   - root: (unique) zero-crossing of f in interval [lower,upper]
#   - info:
#        0 - Successful
#       -1 - Lower bound is larger than upper bound
#       -2 - Search interval does not contain a root, root of even multiplicity, or multiple roots
function bisection(f,lower,upper,xtol,ftol)
    if (lower>upper)
        # Lower bound is larger than upper bound
        return nothing,-1
    end
    flower = f(lower)
    fupper = f(upper)
    if (flower == 0)
        return lower,0
    elseif (fupper == 0)
        return upper, 0
    elseif (flower*fupper>0)
        # Function has same sign in both lower and upper boundary, f has possible
        # no root, multiple roots or a root of even multiplicity in the interval
        return nothing,-2
    end
    while ((upper-lower)>xtol*2)
        middle = (lower+upper)/2
        fmiddle = f(middle)
        if (fmiddle == 0)
            return middle,0
        elseif (fmiddle*flower>0)
            lower = middle
            flower = fmiddle
        else
            upper = middle
            fupper = fmiddle
        end
    end
    return (upper+lower)/2, 0
end



function bisection_(f,lower,upper,xtol,ftol)
  if (upper<lower)
    return nothing
  end
  flower = f(lower)
  fupper = f(upper)
  if(flower*fupper>0)
    # Function has same sign in both lower and upper boundary, f has possible
    # no root, multiple roots or a root of even multiplicity in the interval
    return bisection_result(nothing,f,[lower],[upper])
  elseif (flower == 0)
    return bisection_result(lower,f,[lower],[upper])
  elseif (fupper == 0)
    return bisection_result(upper,f,[lower],[upper])
  end

  n = ceil(log2((upper-lower)/tolerance))-1
  lower_vec = Vector{AbstractFloat}(nothing,n+1)
  upper_vec = Vector{AbstractFloat}(nothing,n+1)
  lower_vec[1] = lower
  upper_vec[1] = upper
  for i = 1:n
    middle = (upper+lower)/2
    fmiddle = f(middle)
    if (fmiddle == 0)
      return bisection_resul(middle,f,lower_vec[1:i],upper_vec[1:i])
    elseif (fmiddle*flower>0)
      lower = middle

      flower = fmiddle
    else # (fmiddle*fupper>0)
      upper = middle
      fupper = fmiddle
    end
    lower_vec[i+1] = lower
    upper_vec[i+1] = upper
  end
  return bisection_result((lower+upper)/2,f,lower_vec,upper_vec)
end
struct bisection_result
  x::AbstractFloat
  f
  lower_vec::Base.AbstractVector{AbstractFloat}
  upper_vec::Base.AbstractVector{AbstractFloat}
end
