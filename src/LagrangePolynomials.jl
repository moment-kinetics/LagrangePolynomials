"""
Lagrange polynomials can be useful for finite element methods on any set of basis points,
as they give a representation of the interpolating function within an element whose
coefficients are just the function values at the grid points.
"""
module LagrangePolynomials

export lagrange_poly,
    lagrange_poly_derivative,
    LagrangePolyData


struct LagrangePolyDataBase
    # the set of points {x_i, i /= j} for each point x_j in {x_j}
    other_nodes::Array{Float64,1}
    # the set of points {x_i, i /= j,k} for each point x_j in {x_j}
    other_nodes_derivative::Array{Float64,2}
    # One over the denominators of the Lagrange polynomials
    # 1/ prod ((x_j - x_i) for i /= j)
    one_over_denominator::Float64
end

struct LagrangePolyData
    # the set of points {x_i} used to construct the interpolant
    x_nodes::Array{Float64,1}
    # precomputed data
    lpoly_data::Array{LagrangePolyDataBase,1}
    """
    Internal constructor function that takes only `x_nodes` as an argument.
    """
    function LagrangePolyData(x_nodes::AbstractArray{Float64,1})
        # Avoid using iterators below for maximum readibility
        ngrid = size(x_nodes,1)
        if ngrid < 2
            error("ERROR: LagrangePolynomials requires ngrid = size(x_nodes,1) >= 2")
        end
        lpoly_data = Array{LagrangePolyDataBase,1}(undef,ngrid)
        for j in 1:ngrid
            other_nodes = Array{Float64,1}(undef,ngrid-1)
            # collect nodes x_i, i /= j
            for i in 1:j-1
                other_nodes[i] = x_nodes[i]
            end
            for i in j+1:ngrid
                other_nodes[i-1] = x_nodes[i]
            end
            # collect nodes x_i i /= j,k
            other_nodes_derivative = Array{Float64,2}(undef,ngrid-2,ngrid-1)
            for i in 1:ngrid-1
                for k in 1:i-1
                    other_nodes_derivative[k,i] = other_nodes[k]
                end
                for k in i+1:ngrid-1
                    other_nodes_derivative[k-1,i] = other_nodes[k]
                end
            end
            # form 1 / Prod ( (x_j - x_i) for j /= i )
            one_over_denominator = 1.0
            for i in 1:j-1
                one_over_denominator /= (x_nodes[j]-x_nodes[i])
            end
            for i in j+1:ngrid
                one_over_denominator /= (x_nodes[j]-x_nodes[i])
            end
            lpoly_data[j] = LagrangePolyDataBase(other_nodes,
                                                    other_nodes_derivative,
                                                    one_over_denominator)
        end
        # return deepycopy of x_nodes to ensure that this data
        # cannot be changed unintentionally
        return new(deepcopy(x_nodes),
                            lpoly_data)
    end
end

"""
    lagrange_poly(jth_lpoly_data, x)

Lagrange polynomial calculation, making use of pre-calculated quantities.

`jth_lpoly_data` contains `other_nodes`, a vector of the grid points in this element where this Lagrange
polynomial is zero (the other nodes than the one where it is 1), and
`one_over_denominator = 1/prod(x0 - n for n in other_nodes)` where `x0` is the grid
point where this Lagrange polynomial is 1.

`x` is the point to evaluate the Lagrange polynomial at.
"""
function lagrange_poly(jth_lpoly_data::LagrangePolyDataBase,
                                 x::Float64)
    other_nodes = jth_lpoly_data.other_nodes
    one_over_denominator = jth_lpoly_data.one_over_denominator
    n = size(other_nodes,1)
    poly = 1.0
    for i in 1:n
        poly *= (x - other_nodes[i])
    end
    poly *= one_over_denominator
    return poly
end

"""
    lagrange_poly_derivative(jth_lpoly_data::LagrangePolyDataBase,
                            x::Float64)

Calculation of the first derivative of a Lagrange polynomial, making use of
pre-calculated quantities.

`x` is the point to evaluate the Lagrange polynomial at.
"""
function lagrange_poly_derivative(jth_lpoly_data::LagrangePolyDataBase,
                                            x::Float64)
    other_nodes = jth_lpoly_data.other_nodes_derivative
    one_over_denominator = jth_lpoly_data.one_over_denominator
    result = 0.0
    n = size(other_nodes,1)
    m = size(other_nodes,2)
    for i in 1:m
        poly = 1.0
        for j in 1:n
            poly *= (x - other_nodes[j,i])
        end
        result += poly
    end
    result *= one_over_denominator
    return result
end

end
