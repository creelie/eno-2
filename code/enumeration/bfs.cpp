// Recurrent configurations of the Manna model on the L x M rectangle, by breadth-first
// search from the all-ones configuration (Appendix A).  Usage: ./bfs L M [print]
#include <bits/stdc++.h>
using namespace std;
int L,M,N; // sites i=x*M+y, x in [0,L), y in [0,M)
typedef vector<uint8_t> C;
unordered_set<uint32_t> outs;
set<C> seen;
void relax(C &c){
  int k=-1; for(int i=0;i<N;i++) if(c[i]>=2){k=i;break;}
  if(k<0){uint32_t m=0; for(int i=0;i<N;i++) if(c[i]) m|=1u<<i; outs.insert(m); return;}
  if(seen.count(c)) return; seen.insert(c);
  int x=k/M,y=k%M;
  for(int dir=0;dir<2;dir++){ C d=c; d[k]-=2;
    if(dir==0){ if(x>0)d[k-M]++; if(x<L-1)d[k+M]++; } else { if(y>0)d[k-1]++; if(y<M-1)d[k+1]++; }
    relax(d);}
}
int main(int argc,char**argv){L=atoi(argv[1]);M=atoi(argv[2]);N=L*M;
  vector<char> vis(1u<<N,0); uint32_t full=(N==32)?~0u:((1u<<N)-1); vector<uint32_t> q{full}; vis[full]=1;
  for(size_t h=0;h<q.size();h++){uint32_t s=q[h];
    for(int i=0;i<N;i++){ C c(N); for(int j=0;j<N;j++) c[j]=(s>>j)&1; c[i]++; outs.clear(); seen.clear(); relax(c);
      for(auto o:outs) if(!vis[o]){vis[o]=1;q.push_back(o);} }}
  int mn=99; map<int,int> hist; for(auto s:q){int p=__builtin_popcount(s); mn=min(mn,p); hist[p]++;}
  int F=(L-1)*(M-1); printf("%dx%d recurrent=%zu of %u  Nmin=%d  ceil(F/2)=%d\n",L,M,q.size(),1u<<N,mn,(F+1)/2);
  if(argc>3){ for(auto s:q) if(__builtin_popcount(s)==mn){ for(int y=M-1;y>=0;y--){for(int x=0;x<L;x++)printf("%d",(s>>(x*M+y))&1);printf("\n");} printf("\n"); } }
}
