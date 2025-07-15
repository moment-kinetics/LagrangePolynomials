module LagrangePolynomialsTests

using Test: @testset, @test
using LagrangePolynomials
using FastGaussQuadrature

function test_polynomial(x;n=2)
    return (x-1)^n
end

function test_polynomial_prime(x;n=2)
    return n*(x-1)^(n-1)
end


function test_lagrange(;ngrid=20::Int64,
                        func::Function=sin,
                        atol::Float64=1.0e-14,
                        label="trig")
    @testset "lagrange_poly" begin
        println("  -- test lagrange_poly ngrid=$(ngrid) type=$(label)")
        x, w = gausslobatto(ngrid)
        lpoly_data = lagrange_poly_data(x)
        f = Array{Float64,1}(undef,ngrid)
        for i in 1:ngrid
            f[i] = func(x[i])
        end
        for jmid in 1:ngrid-1
            x_interp = 0.5*(x[jmid] + x[jmid + 1])
            f_interp_exact = func(x_interp)
            f_interp = 0.0
            for j in 1:ngrid
                jth_lpoly_data = lpoly_data.lpoly_data[j]
                f_interp += f[j]*lagrange_poly(jth_lpoly_data,x_interp)
            end
            interp_err = abs(f_interp - f_interp_exact)
            @test interp_err < atol
        end
    end
end

function test_lagrange_derivative(;ngrid=20::Int64,
                                func::Function=sin,
                                dfunc::Function=cos,
                                atol::Float64=3.0e-14,
                                label="trig")
    @testset "lagrange_poly_derivative" begin
        println("  -- test lagrange_poly_derivative ngrid=$(ngrid) type=$(label)")
        x, w = gausslobatto(ngrid)
        lpoly_data = lagrange_poly_data(x)
        # test function
        f = Array{Float64,1}(undef,ngrid)
        for i in 1:ngrid
            f[i] = func(x[i])
        end
        for jmid in 1:ngrid-1
            x_interp = 0.5*(x[jmid] + x[jmid + 1])
            f_interp_exact = dfunc(x_interp)
            f_interp = 0.0
            for j in 1:ngrid
                jth_lpoly_data = lpoly_data.lpoly_data[j]
                f_interp += f[j]*lagrange_poly_derivative(jth_lpoly_data,x_interp)
            end
            interp_err = abs(f_interp - f_interp_exact)
            @test interp_err < atol
        end
    end
end

function runtests()
    @testset "LagrangePolynomials" begin
        println("test LagrangePolynomials")
        # test trig function
        test_lagrange()
        # test polynomial functions
        for n in 2:10
            function polynomial(x)
                return test_polynomial(x,n=n-1)
            end
            test_lagrange(ngrid=n,
                        func=polynomial,
                        atol=2.0e-13,
                        label="polynomial")
        end
        # test trig function
        test_lagrange_derivative()
        # test polynomial functions
        for n in 2:10
            function polynomial(x)
                return test_polynomial(x,n=n-1)
            end
            function polynomial_prime(x)
                return test_polynomial_prime(x,n=n-1)
            end
            test_lagrange_derivative(ngrid=n,
                                    func=polynomial,
                                    dfunc=polynomial_prime,
                                    atol=3.0e-12,
                                    label="polynomial")
        end
    end
end

end

using .LagrangePolynomialsTests
LagrangePolynomialsTests.runtests()