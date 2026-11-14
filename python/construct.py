"""Literal execution of the layer construction (Proposition layers) with legality checks."""
import itertools

def G(ls):
    a = 1; b = 1
    for l in ls: a *= l // 2; b *= (l - 1) // 2
    return a + b

class Pile:
    def __init__(self, dims):
        self.dims = list(dims)
        self.z = {c: 1 for c in itertools.product(*[range(1, k + 1) for k in dims])}
    def add(self, s): self.z[s] += 1
    def topple(self, s, a):
        if self.z[s] < 2: raise RuntimeError(f"illegal toppling at {s}")
        self.z[s] -= 2
        for e in (-1, 1):
            t = list(s); t[a] += e; t = tuple(t)
            if t in self.z: self.z[t] += 1

def sites(reg): return list(itertools.product(*reg))

def sparsest(P, reg, axes):
    long = [a for a in axes if len(reg[a]) >= 3]
    one = [a for a in axes if len(reg[a]) == 1]
    if one:
        a = one[0]
        for s in sites(reg): P.add(s); P.topple(s, a)
        return
    if len(axes) == 1:
        a = axes[0]; r = reg[a]
        for s in sites(reg):
            if s[a] == r[-1]: P.add(s)
        for x in reversed(r):
            for s in sites(reg):
                if s[a] == x: P.topple(s, a)
        return
    if not long:
        k = len(axes); gray = [i ^ (i >> 1) for i in range(2 ** k)]
        basereg = [range(reg[a][0], reg[a][0] + 1) if a in axes else reg[a] for a in range(len(reg))]
        for base in sites(basereg):
            def v(j):
                t = list(base)
                for b, a in enumerate(axes): t[a] += (gray[j] >> b) & 1
                return tuple(t)
            P.add(v(0))
            for j in range(2 ** k):
                d = gray[j] ^ gray[(j + 1) % 2 ** k]
                a = axes[(d & -d).bit_length() - 1]
                P.topple(v(j), a)
        return
    a = long[0]; r = reg[a]
    sub = list(reg); sub[a] = range(r[0], r[-1] - 1)
    sparsest(P, sub, axes)
    for s in sites(reg):
        if s[a] == r[-1]: P.add(s); P.topple(s, a)
    lay = list(reg); lay[a] = range(r[-1] - 1, r[-1])
    others = [b for b in axes if b != a]
    while True:
        bad = [s for s in sites(lay) if P.z[s] >= 2]
        if not bad: break
        P.topple(bad[0], others[0])
    for s in sites(lay):
        if P.z[s] == 0: P.add(s)
    sparsest(P, lay, others)

def build(dims):
    P = Pile(dims); sparsest(P, [range(1, k + 1) for k in dims], list(range(len(dims))))
    assert all(v <= 1 for v in P.z.values()) and sum(P.z.values()) == G(dims), dims
    return sorted(s for s, v in P.z.items() if v == 1)

if __name__ == "__main__":
    for d in ([9, 8], [5, 5, 5], [4, 4, 4], [4, 4, 5]):
        g = build(d); print(d, len(g), G(d), g)
