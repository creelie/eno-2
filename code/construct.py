"""Explicit legal sequence from the all-ones configuration to a sparsest
recurrent configuration of the Manna sandpile on a box (Section 5 of the paper).
Every toppling is checked for legality; the final configuration is checked to
be stable, and its number of grains is compared with the lower bound."""
import itertools
from math import prod

class Box:
    def __init__(self, dims):
        self.dims = tuple(dims)
        self.z = {s: 1 for s in itertools.product(*[range(n) for n in dims])}
    def add(self, s):
        self.z[s] += 1
    def topple(self, s, axis):
        assert self.z[s] >= 2, ("illegal toppling", s)
        self.z[s] -= 2
        for e in (-1, 1):
            t = list(s); t[axis] += e; t = tuple(t)
            if t in self.z:
                self.z[t] += 1
    def stable(self):
        return all(v <= 1 for v in self.z.values())
    def grains(self, sites=None):
        return sum(self.z[s] for s in (sites if sites is not None else self.z))

def sites_of(region):
    return list(itertools.product(*region))

def relax(box, region, axis):
    """Topple unstable sites of the region along one fixed axis until stable."""
    S = sites_of(region)
    while True:
        bad = [s for s in S if box.z[s] >= 2]
        if not bad:
            return
        box.topple(min(bad), axis)

def fill(box, region):
    for s in sites_of(region):
        if box.z[s] == 0:
            box.add(s)

def chain(box, region, a):
    """1D model along axis a on each line of the region, from all-ones:
    add a grain at the last site of each line and relax (Lemma 5.1)."""
    for s in sites_of(region):
        if s[a] == region[a][-1]:
            box.add(s)
    relax(box, region, a)

def two_by_two(box, region, a, b):
    """The explicit sequence of Lemma 5.3 on a 2x2 region in axes a, b."""
    def site(i, j):
        s = [r[0] for r in region]; s[a] += i; s[b] += j; return tuple(s)
    box.add(site(0, 0)); box.topple(site(0, 0), a)
    box.topple(site(1, 0), b); box.topple(site(1, 1), a); box.topple(site(0, 1), b)

def gray_cube(box, region, axes):
    """Lemma 6.2: a region in which every axis has length two, brought from
    all-ones to a single grain along a Hamiltonian cycle (reflected Gray code)."""
    k = len(axes)
    origin = [r[0] for r in region]
    def site(bits):
        s = list(origin)
        for i, a in enumerate(axes):
            s[a] += bits[i]
        return tuple(s)
    code = [[(g >> i) & 1 for i in range(k)] for g in (j ^ (j >> 1) for j in range(2 ** k))]
    box.add(site(code[0]))
    for j in range(2 ** k):
        u, v = code[j], code[(j + 1) % 2 ** k]
        i = next(i for i in range(k) if u[i] != v[i]) if k > 1 else 0
        box.topple(site(u), axes[i])

def sparsest(box, region, axes):
    """Region is at all-ones; bring it to a sparse recurrent configuration of
    the model using `axes`: Lemmas 5.1-5.3, 6.2 and Proposition 5.4.  An odd
    axis is removed from base 1 (Lemma 5.2); with no odd axis, an even axis of
    length at least 4 is removed from base 2; a cube of side 2 is Lemma 6.2."""
    if len(axes) == 1:
        chain(box, region, axes[0]); return
    if all(len(region[a]) == 2 for a in axes):
        gray_cube(box, region, axes); return
    odd = [a for a in axes if len(region[a]) % 2 == 1]
    if odd:
        a = odd[0]; base = 1
    else:
        a = next(a for a in axes if len(region[a]) >= 4); base = 2
    rest = [b for b in axes if b != a]
    lo, hi = region[a][0], region[a][-1]
    base_region = list(region); base_region[a] = range(lo, lo + base)
    if base == 1:
        for s in sites_of(base_region):       # Lemma 5.2: dump the single layer
            box.add(s); box.topple(s, a)
    else:                                     # same dimension, side a equal to 2
        sparsest(box, base_region, axes)
    c = lo + base
    while c <= hi:                            # Proposition 5.4: add two layers
        l1 = list(region); l1[a] = range(c, c + 1)
        l2 = list(region); l2[a] = range(c + 1, c + 2)
        for s in sites_of(l2):
            box.add(s); box.topple(s, a)
        relax(box, l1, rest[0])               # stabilise layer c with axes != a
        fill(box, l1)                         # back to all-ones on layer c
        sparsest(box, l1, rest)
        c += 2

def run(dims):
    box = Box(dims)
    sparsest(box, [range(n) for n in dims], list(range(len(dims))))
    assert box.stable()
    return box

def lower(dims):
    """Theorem 4.1, and Theorem 6.1 when a side equals 2."""
    if min(dims) == 1:
        return 0
    P = prod(n - 1 for n in dims)
    big = sum(1 for n in dims if n >= 3)
    den = 2 ** (len(dims) - 1) if big == len(dims) else 2 ** big
    return -(-P // den)

def settled(dims):
    """Boxes where the paper proves that the construction is optimal."""
    d = len(dims); odd = sum(n % 2 for n in dims)
    if min(dims) == 1 or d <= 2 or odd >= d - 1:
        return True                                   # Lemma 5.2, Theorems A, B
    if 3 in dims:                                     # a side equal to 3 (Section 6.1)
        rest = list(dims); rest.remove(3)
        if settled(tuple(rest)) and lower(tuple(rest)) == lower(dims):
            return True
    long_ = [n for n in dims if n >= 3]
    if len(long_) < d:                                # a side equal to 2: Theorem 6.1
        ev = [n for n in long_ if n % 2 == 0]
        return len(long_) <= 1 or not ev or (len(ev) == 1 and all(n == 3 for n in long_ if n % 2))
    return False

if __name__ == "__main__":
    cases = [(1,), (5,), (1, 1), (1, 6), (2, 2), (3, 3), (4, 4), (3, 5), (4, 5), (5, 6), (6, 6),
             (7, 6), (9, 8), (10, 12), (2, 2, 2), (4, 2, 2), (5, 2, 2), (2, 2, 6), (2, 3, 4),
             (2, 5, 3), (2, 3, 6), (2, 2, 2, 2), (2, 2, 2, 5), (3, 3, 3), (5, 5, 3), (3, 4, 4),
             (3, 3, 3, 3), (3, 5, 3, 4), (2, 4, 4), (4, 4, 5), (4, 4, 4)]
    ok = True
    for dims in cases:
        b = run(dims)
        g, lb = b.grains(), lower(dims)
        if settled(dims):
            ok &= (g == lb)
            print(dims, "grains", g, "lower bound", lb, "OK" if g == lb else "MISMATCH")
        else:
            print(dims, "grains", g, "lower bound", lb, "open")
    print("all settled cases equal" if ok else "FAIL")
