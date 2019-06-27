### sidi
# Calculate a root of f using Sidi's generalised secand method
# Constraints on f?
#
# Input:
#   - f
#   - initial_points: [xk_n,..,xk]
# Reference
#  Sidi, Avram, "Generalization Of The Secant Method For Nonlinear Equations", Applied Mathematics E-notes 8 (2008), 115–123

# TODO start with lower order methods (eg. only one initial point known)
# TODO if fxk_n,..,fxk already known pass as argument
function sidi(f,initial_points::Array{T,1},ftol,Kmax = 100) where T
    n = length(initial_points)
    points = initial_points
    fpoints = zeros{T}(n)
    for i=1:n
        fpoints[i] = f(points[i])
    end
    while (fpoints[n]>ftol)
        nip = NewtonInterpolatingPolynomial(points,fpoints)
        xnext = points[n]-fpoints[n]/dfeval(nip,points[n])
        points = [points[2:n],xnext];
        fpoints = [fpoints[2:n],f(xnext)];
    end
    return points[n],0
end
