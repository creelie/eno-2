"""Support of the Oslo dynamics: the configurations reachable from (2,...,2) by dropping
grains, compared with the set R (every 0 followed, after a run of 1's, by a 2).
Usage: python3 oslo_support.py   (prints L, |reachable|, |R|, equal?)"""
import itertools,sys
def outcomes(sig,L):
    res=set()
    stack=[(tuple(sig),(1,)+(0,)*(L-1))]
    seen=set()
    while stack:
        s,e=stack.pop()
        if (s,e) in seen: continue
        seen.add((s,e))
        xs=[i for i in range(L) if e[i]>0]
        if not xs: res.add(s); continue
        x=xs[0]; s=list(s); e=list(e); e[x]-=1
        def topple(snew):
            ee=e[:]
            if x>0: ee[x-1]+=1
            if x+1<L: ee[x+1]+=1
            else: ee[L-1]+=1
            return ee
        if s[x]==0:
            s2=s[:]; s2[x]=1; stack.append((tuple(s2),tuple(e)))
        elif s[x]==1:
            s2=s[:]; s2[x]=2; stack.append((tuple(s2),tuple(e)))
            s3=s[:]; s3[x]=0; stack.append((tuple(s3),tuple(topple(s3))))
        else:
            s3=s[:]; s3[x]=1; stack.append((tuple(s3),tuple(topple(s3))))
    return res
def inR(s):
    L=len(s)
    for x in range(L):
        if s[x]==0:
            j=x+1
            while j<L and s[j]==1: j+=1
            if j>=L or s[j]!=2: return False
    return True
for L in range(1,9):
    start=(2,)*L
    R=set([start]); fr=[start]
    while fr:
        n=[]
        for s in fr:
            for t in outcomes(s,L):
                if t not in R: R.add(t); n.append(t)
        fr=n
    Rset=set(s for s in itertools.product(range(3),repeat=L) if inR(s))
    print(L,len(R),len(Rset),R==Rset)
