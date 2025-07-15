# LagrangePolynomials

This package provides an implementation of the Lagrange interpolating polynomials
and their first derivatives. This package is based on code originally implemented in
[moment kinetics](https://github.com/mabarnes/moment_kinetics).
To provide speed when evaluating for the Lagrange polynomials within loops, this package pre-computes quantities to attempt to make the Lagrange polynomial evaluation as rapid as possible.

# Definitions
We define the $`j^{\rm th}`$ interpolating polynomial to be
```math
    l_j(x) = \Pi_{k \neq j} \frac{x - x_k}{x_j - x_k},
```
where the $`x_i`$ are in the set of nodes used to define the interpolation.
The interpolant $`F(x)`$ of $`f(x)`$ is then constructed by
```math
    F(x) = \sum_j f(x_j) l_j(x).
```

The derivative of the $`j^{\rm th}`$ interpolating polynomial is calculated by
```math
    \frac{d l_j}{d x} = \sum_{i} \frac{1}{x_j - x_i} \Pi_{k \neq j,} \frac{x - x_k}{x_j - x_k},
```
and the derivative of the interpolant $`dF/dx`$ by
```math
    \frac{dF}{dx} = \sum_j f(x_j) \frac{d l_j}{d x}.
```

# Example code

To construct the interpolant of a function $`f(x)`$ and the derivative $`f^\prime(x)`$
we can use the following code.
```
using LagrangePolynomials:  lagrange_poly,
                            lagrange_poly_derivative,
                            lagrange_poly_data
using FastGaussQuadrature: gausslobatto

function test(ngrid)
    # construct a set of nodes x on which to interpolate
    # here we choose a Gauss-Legendre-Lobatto grid from FastGaussQuadrature
    x, w = gausslobatto(ngrid)

    # precompute data for interpolation
    lpoly_data = lagrange_poly_data(x)

    # initialise some function to interpolate, here a sine wave
    f = Array{Float64,1}(undef,ngrid)
    for i in 1:ngrid
        f[i] = sin(x[i])
    end

    # choose where to interpolate the data in x
    x_interp = 0.675

    # construct the interpolant
    f_interp = 0.0
    for j in 1:ngrid
        jth_lpoly_data = lpoly_data.lpoly_data[j]
        f_interp += f[j]*lagrange_poly(jth_lpoly_data,x_interp)
    end
    println("f(x) ",f_interp," ", sin(x_interp))
    # construct the interpolant derivative
    f_prime_interp = 0.0
    for j in 1:ngrid
        jth_lpoly_data = lpoly_data.lpoly_data[j]
        f_prime_interp += f[j]*lagrange_poly_derivative(jth_lpoly_data,x_interp)
    end
    println("f'(x) ",f_prime_interp," ", cos(x_interp))
end
test(20)
```
The above code creates the following output.
```
f(x) 0.6248973167277001 0.6248973167276999
f'(x) 0.7807069511324476 0.7807069511324468
```