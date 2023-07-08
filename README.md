# Extended Laurent Polynomials

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.juliahub.com/ExtendedLaurentPolynomials/I3uj7/)

Basic arithmetic, evaluation and substitution for Laurent Polynomials.

## Installation

```julia
(v1.8) pkg> add ExtendedLaurentPolynomials
```

This packege supports Julia 1.8 and later.

## Usage

```julia
julia> using ExtendedLaurentPolynomials
```

## Construction

Construct polynomial from a single number.
```julia
julia> Polynomial(8)
8
```

Construct polynomial from vectors of coefficients and exponents.
```julia
julia> Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2
```

You can also choose another variable for the polynomial.
```julia
julia> Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2]; variable=:x)
x^-2 + 2*x^-1 + 3 + 4*x + 5*x^2
```

Instead of chopping really small numbers using default constant $\epsilon = 10^{-6}$, you can set $\epsilon$ manually.
```julia
julia> Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2]; ϵ=10^(-3))
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2
```

Create a new polynomial using an exsting one with another variable.
```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> Polynomial(p, :y)
y^-2 + 2*y^-1 + 3 + 4*y + 5*y^2
```

## Evaluating and substituting

Basic syntax for evaluating is the following.
```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> p(2.3)
39.7086011342155
```

The $\to$ (***\to***) syntax also can be used for evaluation.
```julia
julia> 2.3 → p
39.7086011342155
```

Substitution is like evaluation, but except of using numbers monomials represented with Julia Exprs are used.
```julia
julia> substitute(p, :(√x))
x^-1 + 2*x^⁻¹/₂ + 3 + 4*x^¹/₂ + 5*x
```

The $\to$ (***\to***) syntax also can be used for substitution.
```julia
julia> :(√x) → p
x^-1 + 2*x^⁻¹/₂ + 3 + 4*x^¹/₂ + 5*x 
```

## Arithmetic

Operaions with two polynomials.
```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> q = Polynomial([1, 2, 3], [-1, 0, 1])
z^-1 + 2 + 3*z 

julia> p + q
z^-2 + 3*z^-1 + 5 + 7*z + 5*z^2 

julia> p - q
z^-2 + z^-1 +  + z + 5*z^2 

julia> p * q
z^-3 + 4*z^-2 + 10*z^-1 + 16 + 22*z + 22*z^2 + 15*z^3 

julia> p ÷ q
(0.2962962962962963z^-1 + 0.22222222222222218 + 1.6666666666666667*z , 0.7037037037037037z^-2 + 1.1851851851851853*z^-1 )

julia> p // q
0.2962962962962963z^-1 + 0.22222222222222218 + 1.6666666666666667*z 

julia> p % q
0.7037037037037037z^-2 + 1.1851851851851853*z^-1 

julia> -p
-z^-2 - 2*z^-1 - 3 - 4*z - 5*z^2 

julia> p / q # evaluates as `p * inv(q)`
ERROR: "z^-1 + 2 + 3*z  is Uninvertible"
```

Operations with a polynomial and a number.
```julia
julia> p * 2
2z^-2 + 4*z^-1 + 6 + 8*z + 10*z^2 

julia> p / 2
0.5z^-2 + z^-1 + 1.5 + 2.0*z + 2.5*z^2 

julia> p ^ 2
z^-4 + 4*z^-3 + 10*z^-2 + 20*z^-1 + 35 + 44*z + 46*z^2 + 40*z^3 + 25*z^4
```

Inverse polynomial can be found only for polynomials of 0 degree
```julia
julia> r = Polynomial([2], [3])
2z^3 

julia> degree.([p, q, r])
3-element Vector{Int64}:
 4
 2
 0
 
julia> inv(r)
0.5z^-3

julia> p / r
0.5z^-5 + z^-4 + 1.5*z^-3 + 2.0*z^-2 + 2.5*z^-1
```

Using $\div$ operator the GCD of two polynomials can be found. As a bonus it returns quotients from the Euler algorithm.
```julia
julia> gcd(p, q)
(0.8701171874999999z^-1 , Any[0.2962962962962963z^-1 + 0.22222222222222218 + 1.6666666666666667*z , 0.18457031250000017z + 2.5312499999999996*z^2 , 0.8087458951656483z^-1 + 1.362098349752671 ])
```

## Some other stuff, I wrote for Wavelet Analysis course

```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> LaurentLeft(p)
z^-1 + 2 + 3*z + 4*z^2 + 5*z^3 

julia> LaurentRight(p)
z^-3 + 2*z^-2 + 3*z^-1 + 4 + 5*z 

julia> LaurentUp(p)
z^-4 + 2*z^-2 + 3 + 4*z^2 + 5*z^4 

julia> LaurentDown(p)
z^-1 + 3.0 + 5.0*z 

julia> InvZ(p)
5-element Vector{Int64}:
 1
 2
 3
 4
 5
```