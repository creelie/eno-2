# Exact checks for the Oslo rice pile part of the paper, in rational arithmetic.
# Run: julia oslo.jl
using Printf

const Conf = NTuple{N,Int8} where N
checks = 0; fails = 0
function check(ok, msg)
    global checks += 1
    ok || (global fails += 1; println("FAIL: ", msg))
end

# law of the final configuration of a complete relaxation of (sigma, eta),
# processing the smallest site with a pending particle (Lemma: order is irrelevant)
function relax(s::Vector{Int8}, e::Vector{Int}, p::Rational{BigInt}, memo)
    key = (Tuple(s), Tuple(e))
    haskey(memo, key) && return memo[key]
    L = length(s); x = findfirst(>(0), e)
    res = Dict{NTuple,Rational{BigInt}}()
    if x === nothing
        res[Tuple(s)] = 1
    else
        e2 = copy(e); e2[x] -= 1
        function toppled(snew)
            ee = copy(e2); x > 1 && (ee[x-1] += 1); ee[min(x+1, L)] += 1
            relax(snew, ee, p, memo)
        end
        branches = Tuple{Rational{BigInt},Dict}[]
        if s[x] == 0
            s2 = copy(s); s2[x] = 1; push!(branches, (1//1, relax(s2, e2, p, memo)))
        elseif s[x] == 1
            s2 = copy(s); s2[x] = 2; push!(branches, (p, relax(s2, e2, p, memo)))
            s3 = copy(s); s3[x] = 0; push!(branches, (1 - p, toppled(s3)))
        else
            s3 = copy(s); s3[x] = 1; push!(branches, (1//1, toppled(s3)))
        end
        for (w, d) in branches, (c, q) in d
            res[c] = get(res, c, 0//1) + w*q
        end
    end
    memo[key] = res
end

g(s) = sum(i*s[i] for i in eachindex(s))
inR(s) = all(s[x] != 0 || (j = findfirst(y -> s[y] != 1, x+1:length(s)); j !== nothing && s[x+j] == 2) for x in eachindex(s))
K(s) = (L = length(s); L*(L+1) - g(s) + (s[end] == 2))

for L in 1:6, p in (1//3, 2//5, 3//4)
    L == 6 && p != 1//3 && continue
    pb = Rational{BigInt}(p)
    memo = Dict()
    confs = [Int8.(collect(c)) for c in Iterators.product((0:2 for _ in 1:L)...)]
    e1 = [1; zeros(Int, L-1)]
    T = Dict(Tuple(c) => relax(c, e1, pb, memo) for c in confs)
    Wv(v) = (out = Dict{NTuple,Rational{BigInt}}();
             for (c, q) in v, (c2, q2) in T[c]; out[c2] = get(out, c2, 0//1) + q*q2; end;
             filter!(kv -> kv[2] != 0, out); out)
    two = Tuple(fill(Int8(2), L)); one_ = Tuple(fill(Int8(1), L))
    Psi = Wv(Dict(two => 1//1))
    check(Wv(Psi) == Psi, "Psi stationary L=$L")
    rec = Set(keys(Psi))
    check(rec == Set(Tuple(c) for c in confs if inR(c)), "recurrent set = R, L=$L")
    # Theorem A: W^k|sigma> = Psi iff k >= K(sigma)
    for c in confs
        v = Dict{NTuple,Rational{BigInt}}(Tuple(c) => 1//1); k = 0
        while v != Psi; v = Wv(v); k += 1; k > L*(L+1) + 2 && break; end
        check(k == K(c), "K L=$L sigma=$c")
    end
    # Corollary B(2): max K on recurrent configurations is L(L+1)/2, attained at 1
    check(maximum(K(c) for c in rec) == L*(L+1)÷2 && K(collect(one_)) == L*(L+1)÷2, "rec max K L=$L")
    p == 1//3 && @printf("  L=%d: %4d configurations, %4d recurrent, Psi(1,...,1) = %s\n", L, 3^L, length(rec), string(Psi[one_])); flush(stdout)
end

println("Number of configurations in R against Fibonacci numbers F_{2L+1}")
function fib(n)
    a, b = big(0), big(1)
    for _ in 1:n; a, b = b, a + b; end
    a
end
for L in 1:12
    nR = count(c -> inR(collect(c)), Iterators.product((0:2 for _ in 1:L)...))
    check(nR == fib(2L+1), "count L=$L")
    @printf("  L=%2d  |R| = %7d  F_%d = %7d\n", L, nR, 2L+1, fib(2L+1)); flush(stdout)
end

println("D^{-1} = (min(x,y)) and the weight identity")
for L in 1:12
    D = zeros(Int, L, L)
    for x in 1:L
        D[x, x] = x < L ? 2 : 1
        x > 1 && (D[x-1, x] = -1); x < L && (D[x+1, x] = -1)
    end
    M = [min(x, y) for x in 1:L, y in 1:L]
    check(D*M == [Int(x == y) for x in 1:L, y in 1:L], "Dinv L=$L")
    check([sum(x*D[x, y] for x in 1:L) for y in 1:L] == [Int(y == L) for y in 1:L], "weight L=$L")
end
@printf("oslo.jl: %d checks, %d failures\n", checks, fails)
fails == 0 || exit(1)
