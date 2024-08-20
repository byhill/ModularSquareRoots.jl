module ModularSquareRoots

using Primes

export sqrtmod
export sqrtmodprime


sqrtmod(n::Integer, m::Integer) = sqrtmod(promote(n, m)...)
sqrtmodprimepower(n::Integer, p::Integer, k::Integer) = sqrtmodprimepower(promote(n, p)..., k)
sqrtmodprime(n::Integer, p::Integer) = sqrtmodprime(promote(n, p)...)


"""
    sqrtmod(n::Integer, m::Integer)

Finds all integers `0 ≤ x < m` that solve the congruence ``x^2 ≡ n \\pmod m``.

Returns an unsorted list.

!!! note
    Calculating the square root modulo a composite number
    is computationally equivalent to integer factorization.
    If you know that `m` is prime, consider using [`sqrtmodprime`](@ref),
    which uses a polynomial time algorithm.

See also [`sqrtmodprime`](@ref).

# Examples
```julia-repl
julia> sqrtmod(4, 5)
2-element Vector{Int64}:
 3
 2

julia> all(powermod(x, 2, 5) == 4 for x in sqrtmod(4, 5))
true

julia> sqrtmod(1240, 289032)
8-element Vector{Int64}:
 107056
 251572
  10712
 155228
 278320
 133804
 181976
  37460

julia> all(powermod(x, 2, 289032) == 1240 for x in sqrtmod(1240, 289032))
true

julia> sqrtmod(23, 200)
Int64[]

julia> !any(powermod(x, 2, 200) == 23 for x in sqrtmod(23, 200))
true
```
"""
function sqrtmod(n::T, m::T) where {T<:Integer}
    m ≤ 0 && throw(DomainError(m, "The modulus `m` must be a positive integer"))
    m == 1 && return T[0]

    roots = T[0]
    newroots = T[]
    mm = one(T)
    for (p, e) in eachfactor(m)
        for a1 in sqrtmodprimepower(n, T(p), e)
            for a2 in roots
                push!(newroots, crt([a1, a2], [p^e, mm], p^e * mm))
            end
        end

        mm *= p^e
        roots, newroots = newroots, roots
        empty!(newroots)
    end

    return roots
end


"""
    sqrtmodprimepower(n::Integer, p::Integer, k::Integer)

For prime `p` and `q = p^k`,
finds all integers `0 ≤ x < q` that solve the congruence ``x^2 ≡ n \\pmod q``.

Returns an unsorted list.

!!! warning
    The behaviour of `sqrtmodprimepower(n, p, k)` is undefined when `p` is not prime.
    This function assumes `p` is a prime number
    and there are no checks to ensure `p` is prime.
    If you can not gurantee that `p` is a prime number, use [`sqrtmod`](@ref) instead.

See also [`sqrtmod`](@ref), [`sqrtmodprime`](@ref).
"""
function sqrtmodprimepower(n::T, p::T, k::Integer) where {T<:Integer}
    k == 1 && return sqrtmodprime(n, p)

    q = p^k
    n = mod(n, q)
    powermod(n, (p - 1) >> 1, p) == -1 && return T[]  # Euler's criterion

    # Use Hensel's lifting lemma
    roots = T[]
    for r in sqrtmodprimepower(n, p, k - 1)
        if mod(2r, p) != 0
            s = mod(r - (r^2 - n) * invmod(2r, p), q)
            push!(roots, s)
        elseif mod(r^2 - n, q) == 0
            for t = 0:p-1
                push!(roots, r + t * p^(k - 1))
            end
        end
    end

    return roots
end


"""
    sqrtmodprime(n::Integer, p::Integer)

For prime `p`,
finds all integers `0 ≤ x < p` that solve the congruence ``x^2 ≡ n \\pmod p``.

Returns an unsorted list.

!!! warning
    The behaviour of `sqrtmodprime(n, p)` is undefined when `p` is not prime.
    This function assumes `p` is a prime number
    and there are no checks to ensure `p` is prime.
    If you can not gurantee that `p` is a prime number, use [`sqrtmod`](@ref) instead.

See also [`sqrtmod`](@ref).

# Examples
```julia-repl
julia> sqrtmodprime(16, 101)
2-element Vector{Int64}:
 97
  4

julia> sqrtmodprime(15, 101)
Int64[]

julia> sqrtmodprime(0, 101)
1-element Vector{Int64}:
 0
```
"""
function sqrtmodprime(n::T, p::T) where {T<:Integer}
    n = mod(n, p)

    # If n is zero, then the only solution is r = 0
    iszero(n) && return [zero(T)]

    # If p = 2 and n is not zero, then the only solution is r = 1
    p == 2 && return [one(T)]

    # If (n, p) fails Euler's criterion, then n is a non-residue mod p
    powermod(n, (p - 1) >> 1, p) == p - 1 && return T[]

    # If p = 3 (mod 4) (and passes Euler's criterion),
    # then sqrt(n) = ±n^((p + 1)/4) (mod p)
    if mod(p, 4) == 3
        r = powermod(n, p >> 2 + 1, p)
        return T[r, p-r]
    end

    # Otherwise, use the Tonelli-shanks algorithm
    (Q, S) = (p - 1, 0)
    while iseven(Q)
        Q >>= 1
        S += 1
    end
    z = first(z for z in 2:p-1 if powermod(z, (p - 1) >> 1, p) == p - 1)

    M = S
    c = powermod(z, Q, p)
    t = powermod(n, Q, p)
    R = powermod(n, (Q + 1) >> 1, p)

    while true
        isone(t) && return [R, p - R]
        i = first(i for i in 1:M-1 if powermod(t, 2^i, p) == 1)

        b = powermod(c, 2^(M - i - 1), p)
        M = i
        c = powermod(b, 2, p)

        (tc, flag) = Base.mul_with_overflow(t, c)
        t = flag ? T(mod(widemul(t, c), p)) : mod(tc, p)

        (Rb, flag) = Base.mul_with_overflow(R, b)
        R = flag ? T(mod(widemul(R, b), p)) : mod(Rb, p)
    end
end


# Chinese Remainder Theorem
function crt(a, n, m::T) where {T<:Integer}
    ans = zero(T)
    for (ai, ni) in zip(a, n)
        u = m ÷ ni
        x = invmod(u, ni) * u

        (y, flag) = Base.mul_with_overflow(ai, x)
        x = flag ? T(mod(widemul(ai, x), m)) : mod(y, m)

        ans = ans < m - x ? ans + x : -m + ans + x  # i.e., ans = mod(ans + x, m) with overflow protection
    end

    return ans
end

end  # module ModularSquareRoots
