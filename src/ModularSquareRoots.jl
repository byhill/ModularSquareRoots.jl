module ModularSquareRoots

using Primes: factor

export sqrtmod
export sqrtmodp


sqrtmod(n::Integer, m::Integer) = sqrtmod(promote(n, m)...)
sqrtmodp(n::Integer, p::Integer) = sqrtmodp(promote(n, p)...)


# Chinese Remainder Theorem
function crt(a, n, m::T) where {T<:Integer}
    ans = zero(T)
    for (ai, ni) in zip(a, n)
        u = m ÷ ni
        x = invmod(u, ni) * u

        (y, flag) = Base.mul_with_overflow(ai, x)
        x = flag ? T(mod(widemul(ai, x), m)) : mod(y, m)

        (y, flag) = Base.add_with_overflow(ans, x)
        ans = flag ? T(mod(widen(ans) + widen(x), m)) : mod(y, m)
    end

    return ans
end


"""
    sqrtmod(n::T, m::T) where {T<:Integer}

Returns an unsorted list of all `0 ≤ x < m` such that ``x^2 ≡ n (mod m)``.

!!! note
    Calculating the square root of a composite number
    is computationally equivalent to integer factorization.
    If you know that `m` is prime,
    then consider using `sqrtmodp(n, m)`,
    which uses a polynomial time algorithm.
"""
function sqrtmod(n::T, m::T) where {T<:Integer}
    m ≤ 0 && throw(DomainError(m, "The modulus `m` must be a positive integer"))
    m == 1 && return T[0]

    factorm = factor(m)
    N = (p^k for (p, k) in factorm)

    roots = T[]
    for A in Iterators.product((sqrtmodq(n, p, T(k)) for (p, k) in factorm)...)
        push!(roots, crt(A, N, m))
    end

    return roots
end


"""
    sqrtmodq(n::Integer, p::Integer, e::Integer)

Let `q = p^k` be a prime power.
Returns an unsorted list of all `0 ≤ x < q` such that `x^2 ≡ n (mod q)`.
Assumes `p` is prime.
"""
function sqrtmodq(n::T, p::T, k::T) where {T<:Integer}
    k == 1 && return sqrtmodp(n, p)

    q = p^k
    n = mod(n, q)
    powermod(n, (p - 1) >> 1, p) == -1 && return T[]  # Euler's criterion

    # Use Hensel's lifting lemma
    roots = T[]
    for r in sqrtmodq(n, p, k - 1)
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
    sqrtmodp(n::Integer, p::Integer)

Returns an unsorted list of all `0 ≤ x < p` such that `x^2 ≡ n (mod p)`.
Assumes `p` is prime.

!!! warning
    This function assumes that `p` is a prime number
    and does not check to ensure that `p` is indeed prime.
    The behaviour of `sqrtmodp(n, p)` is undefined when `p` is not prime.
    Only use this function if you know that `p` is prime.
    If there is a chance that `p` is not prime, use `sqrtmod(n, p)` instead.
"""
function sqrtmodp(n::T, p::T) where {T<:Integer}
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

end  # module ModularSquareRoots
