module ExtendedLaurentPolynomials

import LinearAlgebra: gcd
import Base: (==), (+), (-), (*), (^), (÷), (%), (/), (\), (//)

function unicode_exponent(io, j)
    a = ("⁻", "", "", "⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
    for i in string(j)
        print(io, a[Int(i)-44])
    end
end
function unicode_subscript(io, j)
    a = ("₋", "", "", "₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉")
    for i in string(j)
        print(io, a[Int(i)-44])
    end
end
function Base.show(io::IO, r::Rational)
    unicode_exponent(io, r.num)
    print(io, '/')
    unicode_subscript(io, r.den)
end

"""
    Polynomial{T}(coeffs, exponent; variable)

Creates a Polynomial with `coeffs` of type `T`. Supports `Rational` exponents and substitution.

#Examples

```julia
julia> Polynomial(8)
8

julia> Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2

julia> Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2]; variable=:x)
x^-2 + 2*x^-1 + 3 + 4*x + 5*x^2

julia> Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2]; ϵ=10^(-3))
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2
```
"""
struct Polynomial{T,X}
    coeffs::Dict{Number,T}

    function Polynomial(
        coeffs::AbstractVector{T},
        exponent::AbstractVector;
        variable::Symbol=:z,
        ϵ::Number=10^(-6)
    ) where {T}
        coeffs = collect(coeffs)
        exponent = collect(exponent)
        length(coeffs) != length(exponent) && throw("lengths of `coeffs` and `exponent` must be the same")

        if T <: Number
            coeffs_zeros = abs.(coeffs) .<= ϵ
            deleteat!(coeffs, coeffs_zeros)
            deleteat!(exponent, coeffs_zeros)
        end

        new{T,variable}(
            Dict(
                exponent .=> coeffs
            )
        )
    end

    function Polynomial(p::Polynomial{T,X}, variable::Symbol) where {T,X}
        new{T,variable}(p.coeffs)
    end
end
Polynomial(a₀::T; variable::Symbol=:z) where {T} = Polynomial([a₀], [0]; variable=variable)
Base.copy(p::Polynomial{T,X}) where {T,X} = Polynomial(p, X)

Base.zero(p::Polynomial{T,X}) where {T,X} = Polynomial(zero(T); variable=X)
Base.zero(::Type{Polynomial{T,X}}) where {T,X} = Polynomial(zero(T); variable=X)

Base.one(p::Polynomial{T,X}) where {T,X} = Polynomial(one(T); variable=X)
Base.one(::Type{Polynomial{T,X}}) where {T,X} = Polynomial(one(T); variable=X)

"""
    exponent(p::Polynomial)

Returns exponents of Polynomial `p`.
"""
Base.exponent(p::Polynomial) = p.coeffs |> keys |> collect .|> (r -> typeof(r) <: Rational ? r.den == 1 ? r.num : r : r)

"""
    degree(p::Polynomial)

Computes degree of `p` as difference between maximal and minimal exponent of it.
"""
degree(p::Polynomial) = length(p.coeffs) == 0 ? 0 : begin
    ext = extrema(exponent(p))
    ext[2] - ext[1]
end

(==)(left::Polynomial{T,X}, right::Polynomial{T,X}) where {T,X} = left.coeffs == right.coeffs # length(left.coeffs)==length(right.coeffs)==0 ? false : 

function Base.show(io::IO, p::Polynomial{T,X}) where {T,X}
    if length(p.coeffs) == 0
        print(io, T <: Number ? zero(T) : 0.0)
        return
    end

    exponents = exponent(p) |> sort

    if T <: Number # T is Number
        # first
        if p.coeffs[exponents[1]] != zero(T) && abs(p.coeffs[exponents[1]]) != one(T)
            print(io, p.coeffs[exponents[1]])
        elseif exponents[1] == 0
            print(io, p.coeffs[exponents[1]])
        end
        if p.coeffs[exponents[1]] != zero(T) && exponents[1] != 0
            if p.coeffs[exponents[1]] == -one(T)
                print(io, '-')
            end
            print(io, X)

            if exponents[1] != 1
                print(io, '^', exponents[1])
            end
        end
        print(io, ' ')

        # rest
        for i in 2:length(exponents)
            if i != 2
                if p.coeffs[exponents[i]] > 0
                    print(io, "+ ")
                else
                    print(io, "- ")
                end
            elseif p.coeffs[exponents[1]] != zero(T)
                if p.coeffs[exponents[i]] > 0
                    print(io, "+ ")
                else
                    print(io, "- ")
                end
            end

            if p.coeffs[exponents[i]] != zero(T) && abs(p.coeffs[exponents[i]]) != one(T)
                print(io, p.coeffs[exponents[i]] |> abs)
            end
            if p.coeffs[exponents[i]] != zero(T) && exponents[i] != 0
                if abs(p.coeffs[exponents[i]]) != one(T)
                    print(io, '*')
                end
                print(io, X)
                if exponents[i] != 1
                    print(io, '^', exponents[i])
                end
            end
            print(io, ' ')
        end
    else # T is Symbol (probably)
        # first
        print(io, p.coeffs[exponents[1]])
        if exponents[1] != 0
            print(io, '*', X)
            if exponents[1] != 1
                print(io, '^', exponents[1])
            end
        end
        print(io, ' ')

        # rest
        for i in 2:length(exponents)
            print(io, "+ ", p.coeffs[exponents[i]])
            if exponents[i] != 0
                print(io, '*', X)
                if exponents[i] != 1
                    print(io, '^', exponents[i])
                end
            end
            print(io, ' ')
        end
    end
end

"""
    substitute(p::Polynomial, expr::Expr)
    (→)(p::Polynomial, expr::Expr)

Substitutes `expr` for `p.variable`. `expr` can be only `√x`, `a*x^b` or `√(a*x^b)`.

# Examples

```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2

julia> substitute(p, :(√x))
x^-1 + 2*x^⁻¹/₂ + 3 + 4*x^¹/₂ + 5*x

julia> :(√x) → p
x^-1 + 2*x^⁻¹/₂ + 3 + 4*x^¹/₂ + 5*x
```
"""
(→)(variable::Symbol, p::Polynomial) = Polynomial(p, variable)
function substitute(p::Polynomial{T,X}, expr::Union{Expr,Symbol}) where {T<:Number,X}
    if typeof(expr) == Symbol
        return expr → p
    elseif expr.args[1] == :√
        coeffs = p.coeffs |> values |> collect
        exponents = p.coeffs |> keys |> collect .|> (e -> e // 2)
        exponents = Rational <: eltype(exponents) ? exponents .|> (r -> r.den == 1 ? r.num : r) : exponents
        return substitute(
            Polynomial(coeffs, exponents; variable=X),
            expr.args[2]
        )
    elseif expr.args[1] == :*
        exponents = p.coeffs |> keys |> collect
        coeffs = p.coeffs |> values |> collect .|> Float64
        coeffs .*= [expr.args[2]^Rational(e) for e in exponents]
        return substitute(
            Polynomial(coeffs, exponents; variable=X),
            expr.args[3]
        )
    elseif expr.args[1] == :^
        coeffs = p.coeffs |> values |> collect
        exponents = p.coeffs |> keys |> collect .|> (e -> e * expr.args[3])
        exponents = Rational <: eltype(exponents) ? exponents .|> (r -> r.den == 1 ? r.num : r) : exponents
        return substitute(
            Polynomial(coeffs, exponents; variable=X),
            expr.args[2]
        )
    end

    return p
end
(→)(expr::Expr, p::Polynomial) = substitute(p, expr)

"""
    (p::Polynomial)(x::Number)
    (→)(x::Number, p::Polynomial)

Evaluates `p` of `x`.

# Examples

```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> p(2.3)
39.7086011342155

julia> 2.3 → p
39.7086011342155
```
"""
(p::Polynomial)(x::Number) =
    let exponents = exponent(p), result = zero(x)
        for e in exponents
            result += p.coeffs[e] * x^Rational(e)
        end
        return result
    end
(→)(x::Number, p::Polynomial) = p(x)

function +(left::Polynomial{T1,X}, right::Polynomial{T2,X}) where {T1<:Number,T2<:Number,X}
    left_exp = exponent(left)
    right_exp = exponent(right)
    exponents = convert.(Rational, union(left_exp, right_exp))

    coeffs = Dict(exponents .=> zeros(promote_type(T1, T2), length(exponents)))
    for le in left_exp
        coeffs[le] += left.coeffs[le]
    end
    for re in right_exp
        coeffs[re] += right.coeffs[re]
    end

    return Polynomial(
        coeffs |> values |> collect,
        coeffs |> keys |> collect;
        variable=X
    )
end

function -(left::Polynomial{T1,X}, right::Polynomial{T2,X}) where {T1<:Number,T2<:Number,X}
    left_exp = exponent(left)
    right_exp = exponent(right)
    exponents = convert.(Rational, union(left_exp, right_exp))

    coeffs = Dict(exponents .=> zeros(promote_type(T1, T2), length(exponents)))
    for le in left_exp
        coeffs[le] += left.coeffs[le]
    end
    for re in right_exp
        coeffs[re] -= right.coeffs[re]
    end

    return Polynomial(
        coeffs |> values |> collect,
        coeffs |> keys |> collect;
        variable=X
    )
end

function *(left::Polynomial{T1,X}, right::Polynomial{T2,X}) where {T1,T2,X}
    left_exp = exponent(left)
    right_exp = exponent(right)
    exponents = convert.(Rational, unique(vcat([
        right_exp .+ le
        for le in left_exp
    ]...)))

    T = promote_type(T1, T2)
    coeffs = Dict(exponents .=> zeros(T == Any ? Float64 : T, length(exponents)))

    for le in left_exp, re in right_exp
        coeffs[le+re] += left.coeffs[le] * right.coeffs[re]
    end

    return length(coeffs) != 0 ? Polynomial(
        coeffs |> values |> collect,
        coeffs |> keys |> collect;
        variable=X
    ) : Polynomial([0], [0]; variable=X)
end
function *(left::Polynomial{T1,X}, right::T2) where {T1<:Number,T2<:Number,X}
    return Polynomial(
        (left.coeffs |> values |> collect) * right,
        left.coeffs |> keys |> collect;
        variable=X
    )
end
(*)(left::Number, right::Polynomial) = right * left
(-)(p::Polynomial{T,X}) where {T,X} = -1 * p

(^)(left::Polynomial{T,X}, right::Int) where {T,X} = Base.power_by_squaring(left, right)

function ÷(left::Polynomial{T1,X}, right::Polynomial{T2,X}) where {T1<:Number,T2<:Number,X}
    div_degree = degree(left) - degree(right)
    div = Polynomial([zero(promote_type(T1, T2))], [0]; variable=X)

    d_right = maximum(exponent(right))
    a_right = right.coeffs[d_right]

    if div_degree == 0
        d_left = maximum(exponent(left))
        a_left = left.coeffs[d_left]
        buf = Polynomial([a_left / a_right], [d_left - d_right]; variable=X)

        div += buf
        left -= right * buf
    end

    while degree(div) < div_degree
        d_left = maximum(exponent(left))
        a_left = left.coeffs[d_left]
        buf = Polynomial([a_left / a_right], [d_left - d_right]; variable=X)

        div += buf
        left -= right * buf
    end
    return div, left
end

(//)(left::Polynomial, right::Polynomial) = ÷(left, right)[1]
(%)(left::Polynomial, right::Polynomial) = ÷(left, right)[2]

function /(left::Polynomial{T1,X}, right::T2) where {T1<:Number,T2<:Number,X}
    return Polynomial(
        (left.coeffs |> values |> collect) / right,
        left.coeffs |> keys |> collect;
        variable=X
    )
end
(\)(left::Number, right::Polynomial) = right / left

function Base.inv(p::Polynomial{T,X}) where {T,X}
    (p == zero(p) || length(p.coeffs) != 1) && throw("$p is Uninvertible")

    coeff = p.coeffs |> collect |> first
    return Polynomial(
        [1 / last(coeff)],
        [-first(coeff)];
        variable=X
    )
end
(/)(left::Number, right::Polynomial) = left * inv(right)
(/)(left::Polynomial{T1,X}, right::Polynomial{T2,X}) where {T1,T2,X} = left * inv(right)

function gcd(a::Polynomial{T1,X}, b::Polynomial{T2,X}) where {T1<:Number,T2<:Number,X}
    Q = []
    while b != zero(b)
        q, r = a ÷ b
        push!(Q, q)
        a, b = b, r
    end
    return a, Q
end

"""
	LaurentLeft(p::Polynomial)

Shifts polynomial `p` to the left.

# Examples

```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> LaurentLeft(p)
z^-1 + 2 + 3*z + 4*z^2 + 5*z^3 
```
"""
LaurentLeft(p::Polynomial{T,X}) where {T,X} = p * Polynomial([1], [1]; variable=X)

"""
	LaurentRight(p::Polynomial)

Shifts polynomial `p` to the right.

# Examples

```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> LaurentRight(p)
z^-3 + 2*z^-2 + 3*z^-1 + 4 + 5*z
```
"""
LaurentRight(p::Polynomial{T,X}) where {T,X} = p * Polynomial([1], [-1]; variable=X)

"""
	LaurentUp(p::Polynomial)

Upscales p.

# Examples

```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> LaurentUp(p)
z^-4 + 2*z^-2 + 3 + 4*z^2 + 5*z^4
```
"""
LaurentUp(p::Polynomial{T,X}) where {T,X} = :($X^2) → p

"""
	LaurentDown(p::Polynomial)

Downscales p.

# Examples

```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> LaurentDown(p)
z^-1 + 3.0 + 5.0*z
```
"""
LaurentDown(p::Polynomial{T,X}) where {T,X} = ((:(√$X) → p) + (:(-1√$X) → p)) / 2

"""
	InvZ(p::Polynomial)

Performs inverse z-transform on polynomial `p` that has only integer exponents.

# Examples

```julia
julia> p = Polynomial([1, 2, 3, 4, 5], [-2, -1, 0, 1, 2])
z^-2 + 2*z^-1 + 3 + 4*z + 5*z^2 

julia> InvZ(p)
5-element Vector{Int64}:
 1
 2
 3
 4
 5
```
"""
function InvZ(p::Polynomial{T,X}, exponents::AbstractVector{Int}=sort(exponent(p))) where {T,X}
    eltype(exponents) == Int && map(
        e -> e in keys(p.coeffs) ? p.coeffs[e] : zero(T),
        exponents
    )
end

export Polynomial, degree, →, substitute, LaurentDown, LaurentUp, LaurentLeft, LaurentRight, InvZ

end # module ExtLaurentPolynomials