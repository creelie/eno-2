// The same search on a box of dimension up to four.  Usage: ./bfsd l1 l2 [l3 [l4]]
// "lower" is the bound of Theorem 4.1 (Theorem 6.1 is sharper when a side equals 2).
#include <bits/stdc++.h>
using namespace std;
int D, Ls[4], N, stride[4]; typedef vector<uint8_t> C;
unordered_set<uint32_t> outs; set<C> seen;
void relax(C &c){int k=-1; for(int i=0;i<N;i++) if(c[i]>=2){k=i;break;}
  if(k<0){uint32_t m=0; for(int i=0;i<N;i++) if(c[i]) m|=1u<<i; outs.insert(m); return;}
  if(seen.count(c)) return; seen.insert(c);
  for(int a=0;a<D;a++){ C d=c; d[k]-=2; int co=(k/stride[a])%Ls[a];
    if(co>0) d[k-stride[a]]++; if(co<Ls[a]-1) d[k+stride[a]]++; relax(d);} }
int main(int argc,char**argv){ D=argc-1; N=1; for(int a=0;a<D;a++){Ls[a]=atoi(argv[a+1]); stride[a]=N; N*=Ls[a];}
  vector<char> vis(1u<<N,0); uint32_t full=(1u<<N)-1; vector<uint32_t> q{full}; vis[full]=1;
  for(size_t h=0;h<q.size();h++){uint32_t s=q[h]; for(int i=0;i<N;i++){ C c(N); for(int j=0;j<N;j++) c[j]=(s>>j)&1; c[i]++; outs.clear(); seen.clear(); relax(c);
      for(auto o:outs) if(!vis[o]){vis[o]=1;q.push_back(o);} }}
  int mn=99; for(auto s:q) mn=min(mn,__builtin_popcount(s)); long P=1; for(int a=0;a<D;a++) P*=Ls[a]-1; long den=1L<<(D-1);
  printf("dims"); for(int a=0;a<D;a++) printf(" %d",Ls[a]); printf(": recurrent=%zu Nmin=%d lower=%ld\n",q.size(),mn,(P+den-1)/den);
}
