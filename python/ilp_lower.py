"""Integer programme of Appendix A: the fewest grains compatible with the conclusions of
Lemmas last and charge (each grain is delivered along at most one edge at its site; every
cell contains a delivered grain on one of its edges or has all corners untouched).
Usage: python3 ilp_lower.py "[(4,4,4),(4,4,5)]"   (needs numpy and scipy >= 1.9)"""
import itertools, sys, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
def lb(dims, useU=True):
    d=len(dims)
    sites=list(itertools.product(*[range(n) for n in dims])); S=set(sites)
    edges=[]
    for s in sites:
        for i in range(d):
            t=list(s); t[i]+=1; t=tuple(t)
            if t in S: edges.append((s,t))
    cells=list(itertools.product(*[range(n-1) for n in dims]))
    var=[]  # (type, data)
    yidx={}
    for k,(a,b) in enumerate(edges):
        for h in (a,b): yidx[(h,k)]=len(var); var.append(('y',h,k))
    uidx={}; cidx={}
    if useU:
        for s in sites: uidx[s]=len(var); var.append(('u',s))
        for Q in cells: cidx[Q]=len(var); var.append(('c',Q))
    n=len(var)
    cost=np.array([0.0 if v[0]=='c' else 1.0 for v in var])
    rows=[]; lo=[]; hi=[]
    def add(coef,l,h): rows.append(coef); lo.append(l); hi.append(h)
    for s in sites:
        coef={i:1 for (h,k),i in yidx.items() if h==s}
        if useU: coef[uidx[s]]=1
        add(coef,-np.inf,1)
    for Q in cells:
        cs=set(tuple(Q[i]+o[i] for i in range(d)) for o in itertools.product((0,1),repeat=d))
        coef={i:1 for (h,k),i in yidx.items() if edges[k][0] in cs and edges[k][1] in cs}
        if useU:
            coef[cidx[Q]]=1
            for s in cs: add({cidx[Q]:1,uidx[s]:-1},-np.inf,0)
        add(coef,1,np.inf)
    A=lil_matrix((len(rows),n))
    for r,coef in enumerate(rows):
        for i,v in coef.items(): A[r,i]=v
    res=milp(cost,constraints=LinearConstraint(A.tocsr(),lo,hi),integrality=np.ones(n),bounds=Bounds(0,1),options={'time_limit':900})
    return round(res.fun) if res.x is not None else None, res.status
for dims in eval(sys.argv[1]):
    print(dims, lb(dims), flush=True)
