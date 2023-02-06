module ModularSquareRoots

using Primes: factor

export modularsquareroots
export sqrtmod
export sqrtmodq
export sqrtmodp

function modularsquareroots(n::Integer, m::Integer)
    solns = []
    n = mod(n, m)
    for x = 0:m÷2
        powermod(x, 2, m) == n && (push!(solns, x), push!(solns, mod(-x, m)))
    end
    return solns
end

# Chinese Remainder Theorem
function crt(a, n, m)
    return mod(sum(ai * invmod(m ÷ ni, ni) * (m ÷ ni) for (ni, ai) in zip(n, a)), m)
end


"""
    sqrtmod(n::Integer, m::Integer)

Returns an unsorted list of all `0 ≤ x < m` such that `x^2 ≡ n (mod p)`.
"""
function sqrtmod(n::T, m::T) where {T<:Integer}
    factorm = factor(m)
    N = (p^k for (p, k) in factorm)

    roots = T[]
    for A in Iterators.product((sqrtmodq(n, p, k) for (p, k) in factorm)...)
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
        r = powermod(n, (p + 1) >> 2, p)
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
