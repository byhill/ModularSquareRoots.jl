module ModularSquareRoots

using Primes

export sqrtmod
export sqrtmodprime


sqrtmod(n::Integer, m::Integer) = sqrtmod(promote(n, m)...)
sqrtmodprime(n::Integer, p::Integer) = sqrtmodprime(promote(n, p)...)


"""
    sqrtmod(n::Integer, m::Integer)

Finds all `0 ≤ x < m` that solve the congruence ``x^2 ≡ n (mod m)``.

Returns an unsorted list.

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

!!! note
    Calculating the square root of a composite number
    is computationally equivalent to integer factorization.
    If you know that `m` is prime,
    consider using `sqrtmodp(n, m)`,
    which uses a polynomial time algorithm.
"""
function sqrtmod(n::T, m::T) where {T<:Integer}
    m ≤ 0 && throw(DomainError(m, "The modulus `m` must be a positive integer"))
    m == 1 && return T[0]

    roots = T[0]
    newroots = T[]
    mm = one(T)
    for (p, e) in eachfactor(m)
        for a1 in _sqrtmodq(n, T(p), e)
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
    _sqrtmodq(n::Integer, p::Integer, e::Integer)

Let `q = p^k` be a prime power.
Returns an unsorted list of all `0 ≤ x < q` such that `x^2 ≡ n (mod q)`.
Assumes `p` is prime.
"""
function _sqrtmodq(n::T, p::T, k::Integer) where {T<:Integer}
    k == 1 && return sqrtmodprime(n, p)

    q = p^k
    n = mod(n, q)
    powermod(n, (p - 1) >> 1, p) == -1 && return T[]  # Euler's criterion

    # Use Hensel's lifting lemma
    roots = T[]
    for r in _sqrtmodq(n, p, k - 1)
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

Returns an unsorted list of all `0 ≤ x < p` such that `x^2 ≡ n (mod p)`.
Assumes `p` is prime.

!!! warning
    This function assumes that `p` is a prime number
    and does not check to ensure that `p` is indeed prime.
    The behaviour of `sqrtmodprime(n, p)` is undefined when `p` is not prime.
    Only use this function if you know that `p` is prime.
    If there is a chance that `p` is not prime, use `sqrtmod(n, p)` instead.
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
