function inverse_quadratic_interpolation(a,b,c,fa,fb,fc)
    # TODO improved implementation: http://mathworld.wolfram.com/BrentsMethod.html
    return fb*fc/((fa-fb)*(fa-fc))*a + fa*fc/((fb-fa)*(fb-fc))*b + fa*fb/((fc-fa)*(fc-fb))
end

struct NewtonInterpolatingPolynomial
    points<:Array
    fval<:Array
end

# TODO efficient evaluation
function feval(p::NewtonInterpolatingPolynomial,x)

end
# TODO check and test
# TODO efficient evaluation
function dfeval(p::NewtonInterpolatingPolynomial,x)
    n = length(p.points)
    DD = zeros{typeof(p[1])}(n,n)
    DD[:,1] = fval
    for i=2:n
        for j=i:n
            DD[j,i] = (DD[j,i-1]-DD[j-1,i-1])/(points[j]-points[j-i+1])
        end
    end
    f = DD[n,n]
    df = 0
    for k = n-1:-1:1
        df = df*(x-point[k])+p
        f = f*point[k]+DD[k,k]
    end
    return df
end
#TODO
Base.show(io::IO,p::NewtonInterpolatingPolynomial) = print(io,'NIP')
