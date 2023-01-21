module ModularSquareRoots

export sqrtmodp


function modularsquareroots(n::Integer, m::Integer)
    solns = []
    n = mod(n, m)
    for x = 0:m÷2
        powermod(x, 2, m) == n && (push!(solns, x), push!(solns, -x))
    end
    return solns
end


"""
    sqrtmodp(n::Integer, p::Integer)

Returns a list of all `0 ≤ x < p` such that `x^2 ≡ n (mod p)`.
Assumes `p` is prime and `0 ≤ n < p`.
"""
function sqrtmodp(n::Integer, p::Integer)
    T = typeof(p)

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
        return [r, p - r]
    end

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
        iszero(t) && error("n is a non-residue (mod p) that wasn't caught")
        isone(t) && return [R, p - R]
        i = first(i for i in 1:M-1 if powermod(t, 2^i, p) == 1)
        b = powermod(c, 2^(M - i - 1), p)
        M = i
        c = powermod(b, 2, p)
        t = mod(t * c, p)  # TODO figure out type-instability warning
        R = mod(R * b, p)
    end
end

end  # module ModularSquareRoots
