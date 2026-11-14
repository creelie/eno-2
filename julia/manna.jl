# Exhaustive and constructive checks for the Manna sandpile part of the paper.
# Run: julia manna.jl
using Printf

# ---------- boxes: sites, neighbours, configurations as bit masks ----------
struct Geo
    dims::Vector{Int}
    n::Int
    nb::Vector{Vector{Tuple{Int,Int}}}   # nb[s][axis] = (left, right), 0 = outside
end
function geo(dims; periodic=Int[])
    n = prod(dims)
    idx(c) = 1 + sum((c[i]-1)*prod(dims[1:i-1]) for i in eachindex(dims))
    coords = vec(collect(Iterators.product((1:k for k in dims)...)))
    nb = [Tuple{Int,Int}[] for _ in 1:n]
    for c in coords, a in eachindex(dims)
        lr = Int[]
        for e in (-1, 1)
            t = collect(c); t[a] += e
            if a in periodic; t[a] = mod1(t[a], dims[a]); end
            push!(lr, 1 <= t[a] <= dims[a] ? idx(t) : 0)
        end
        push!(nb[idx(collect(c))], (lr[1], lr[2]))
    end
    Geo(collect(dims), n, nb)
end

# all stable outcomes of relaxing an unstable height vector, toppling the smallest
# unstable site along every axis in turn (the support of one step of the chain)
function outcomes!(out::Set{UInt64}, g::Geo, h0::Vector{Int8})
    seen = Set{Vector{Int8}}([h0]); stack = [h0]
    while !isempty(stack)
        h = pop!(stack)
        s = findfirst(>=(2), h)
        if s === nothing
            m = UInt64(0); for i in 1:g.n; h[i] == 1 && (m |= UInt64(1) << (i-1)); end
            push!(out, m); continue
        end
        for (l, r) in g.nb[s]
            h2 = copy(h); h2[s] -= 2
            l > 0 && (h2[l] += 1); r > 0 && (h2[r] += 1)
            h2 in seen || (push!(seen, h2); push!(stack, h2))
        end
    end
end

function recurrent_set(g::Geo)
    full = (UInt64(1) << g.n) - 1
    seen = Set{UInt64}([full]); frontier = [full]
    while !isempty(frontier)
        nxt = UInt64[]
        for m in frontier, s in 1:g.n
            h = Int8[(m >> (i-1)) & 1 for i in 1:g.n]; h[s] += 1
            out = Set{UInt64}(); outcomes!(out, g, h)
            for t in out
                if !(t in seen); push!(seen, t); push!(nxt, t); end
            end
        end
        frontier = nxt
    end
    seen
end

G(ls) = prod(l ÷ 2 for l in ls) + prod((l - 1) ÷ 2 for l in ls)

checks = 0; fails = 0
function check(ok, msg)
    global checks += 1
    ok || (global fails += 1; println("FAIL: ", msg))
end

println("Exhaustive recurrent sets on boxes")
for dims in ([2,2],[2,3],[2,4],[2,5],[3,3],[3,4],[3,5],[4,4],[2,2,2],[3,2,2],[4,2,2])
    g = geo(dims); R = recurrent_set(g)
    mn = minimum(count_ones, R)
    @printf("  %-10s recurrent %6d of %6d   min grains %d   G = %d\n", join(dims, "x"), length(R), 2^g.n, mn, G(dims)); flush(stdout)
    check(mn == G(dims), "box $dims")
end

println("Exhaustive recurrent sets on cylinders [l] x Z/mZ")
for (l, m) in ((2,3),(2,4),(2,5),(2,6),(3,3),(3,4),(3,5),(4,3),(5,3),(4,4))
    g = geo([l, m]; periodic=[2]); R = recurrent_set(g)
    mn = minimum(count_ones, R)
    @printf("  %dx%d  recurrent %6d of %6d   min grains %d   bound %d\n", l, m, length(R), 2^g.n, mn, cld((l-1)*m, 2)); flush(stdout)
    check(mn == cld((l-1)*m, 2), "cylinder $l x $m")
end

# ---------- the constructions of the upper bound, executed with legality checks ----------
mutable struct Pile
    dims::Vector{Int}
    z::Dict{Vector{Int},Int}
end
Pile(dims) = Pile(dims, Dict(collect(c) => 1 for c in Iterators.product((1:k for k in dims)...)))
add!(P, s) = (P.z[s] += 1)
function topple!(P, s, a)
    P.z[s] >= 2 || error("illegal toppling at $s")
    P.z[s] -= 2
    for e in (-1, 1)
        t = copy(s); t[a] += e
        haskey(P.z, t) && (P.z[t] += 1)
    end
end
# region = vector of ranges; returns list of sites
sites(reg) = [collect(c) for c in Iterators.product(reg...)]

# bring the sub-box `reg` (currently full on reg) to a sparsest configuration,
# toppling only along the axes listed in `axes` (the axes of the sub-box)
function sparsest!(P, reg, axes)
    long = [a for a in axes if length(reg[a]) >= 3]
    if any(length(reg[a]) == 1 for a in axes)
        a = first(a for a in axes if length(reg[a]) == 1)
        for s in sites(reg); add!(P, s); topple!(P, s, a); end      # Lemma (thin)
        return
    end
    if length(axes) == 1
        a = axes[1]; r = reg[a]
        for s in sites(reg); s[a] == last(r) && add!(P, s); end    # Lemma (chain)
        for x in reverse(r), s in sites(reg)
            s[a] == x && topple!(P, s, a)
        end
        return
    end
    if isempty(long)                                               # all sides two: Gray cycle
        k = length(axes)
        gray = [i ⊻ (i >> 1) for i in 0:2^k-1]
        for base in sites([a in axes ? (first(reg[a]):first(reg[a])) : reg[a] for a in eachindex(reg)])
            v(j) = (t = copy(base); for (b, a) in enumerate(axes); t[a] += (gray[j+1] >> (b-1)) & 1; end; t)
            add!(P, v(0))
            for j in 0:2^k-1
                d = gray[j+1] ⊻ gray[mod(j+1, 2^k)+1]
                a = axes[trailing_zeros(d)+1]
                topple!(P, v(j), a)
            end
        end
        return
    end
    a = long[1]; r = reg[a]
    # Proposition (layers): first length(r)-2 layers, then layer r[end] emptied,
    # layer r[end-1] relaxed along the other axes, refilled, and made sparsest
    sub = copy(reg); sub[a] = first(r):last(r)-2
    sparsest!(P, sub, axes)
    for s in sites(reg); s[a] == last(r) && (add!(P, s); topple!(P, s, a)); end
    lay = copy(reg); lay[a] = last(r)-1:last(r)-1
    others = [b for b in axes if b != a]
    while true
        bad = [s for s in sites(lay) if P.z[s] >= 2]
        isempty(bad) && break
        topple!(P, bad[1], others[1])
    end
    for s in sites(lay); P.z[s] == 0 && add!(P, s); end
    sparsest!(P, lay, others)
end

println("Constructions on all boxes with d <= 4, sides <= 7 and at most 600 sites")
ncons = 0
for d in 1:4, ls in Iterators.product((1:7 for _ in 1:d)...)
    dims = collect(ls); issorted(dims) || continue; prod(dims) <= 600 || continue
    P = Pile(dims); sparsest!(P, [1:k for k in dims], collect(1:d))
    stable = all(v <= 1 for v in values(P.z))
    check(stable && sum(values(P.z)) == G(dims), "construction $dims")
    global ncons += 1
end
println("  $ncons boxes, each history legal and ending stable with G grains")

println("Same-parity family of cells")
for d in 1:4, ls in Iterators.product((1:6 for _ in 1:d)...)
    dims = collect(ls); issorted(dims) || continue
    cells = [collect(c) for c in Iterators.product((1:k-1 for k in dims)...)]
    F = [c for c in cells if all(isodd, c) || all(iseven, c)]
    check(length(F) == G(dims), "family size $dims")
    for i in eachindex(F), j in i+1:length(F)
        common = prod(length(intersect(F[i][a]:F[i][a]+1, F[j][a]:F[j][a]+1)) for a in 1:d)
        check(common <= 1, "family overlap $dims")
    end
end

println("G satisfies the layer recursion G(l+2, l') = G(l, l') + G(l')")
for d in 2:5, ls in Iterators.product((1:8 for _ in 1:d)...)
    dims = collect(ls)
    check(G([dims[1]+2; dims[2:end]]) == G(dims) + G(dims[2:end]), "recursion $dims")
end
println("Two-dimensional formula G(l,m) = ceil((l-1)(m-1)/2)")
for l in 1:60, m in 1:60
    check(G([l, m]) == cld((l-1)*(m-1), 2), "2D $l $m")
end
@printf("manna.jl: %d checks, %d failures\n", checks, fails)
fails == 0 || exit(1)
