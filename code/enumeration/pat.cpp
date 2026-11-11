// Dhar's patterns of Figures 1a-1g and 2 (Physica A 270 (1999) 69-81): whether each
// occurs in a recurrent configuration of its bounding box, and the modified burning test.
#include <bits/stdc++.h>
using namespace std;
int L,M,N; typedef vector<uint8_t> C;
unordered_set<uint32_t> outs; set<C> seen;
void relax(C &c){int k=-1; for(int i=0;i<N;i++) if(c[i]>=2){k=i;break;}
  if(k<0){uint32_t m=0; for(int i=0;i<N;i++) if(c[i]) m|=1u<<i; outs.insert(m); return;}
  if(seen.count(c)) return; seen.insert(c); int x=k/M,y=k%M;
  for(int dir=0;dir<2;dir++){ C d=c; d[k]-=2;
    if(dir==0){ if(x>0)d[k-M]++; if(x<L-1)d[k+M]++; } else { if(y>0)d[k-1]++; if(y<M-1)d[k+1]++; } relax(d);} }
vector<uint32_t> rec(){ vector<char> vis(1u<<N,0); uint32_t full=(1u<<N)-1; vector<uint32_t> q{full}; vis[full]=1;
  for(size_t h=0;h<q.size();h++){uint32_t s=q[h]; for(int i=0;i<N;i++){ C c(N); for(int j=0;j<N;j++) c[j]=(s>>j)&1; c[i]++; outs.clear(); seen.clear(); relax(c);
      for(auto o:outs) if(!vis[o]){vis[o]=1;q.push_back(o);} }} return q;}
int main(){ vector<pair<string,vector<string>>> P={{"1a",{"00","00"}},{"1b",{".00","010","00."}},{"1c",{"0000","0110",".00."}},
 {"1d",{"1010","0000","0101"}},{"1e",{"0000.","01110",".0000"}},{"1f",{"..00",".010","010.","00.."}},{"1g",{"..00","0010","0110","0000"}},{"2",{"000.","0110",".000"}}};
 map<pair<int,int>,vector<uint32_t>> cache;
 for(auto&[name,rows]:P){ int H=rows.size(), W=rows[0].size(); L=W;M=H;N=L*M;
   if(!cache.count({L,M})) cache[{L,M}]=rec(); auto&R=cache[{L,M}];
   uint32_t care=0,val=0; for(int r=0;r<H;r++)for(int x=0;x<W;x++){char ch=rows[r][x]; int y=H-1-r; int i=x*M+y; if(ch!='.'){care|=1u<<i; if(ch=='1')val|=1u<<i;}}
   auto occurs=[&](uint32_t cr,uint32_t v){for(auto s:R) if((s&cr)==v) return true; return false;};
   bool occ=occurs(care,val);
   // minimality: drop each specified site
   string minimal="yes"; for(int i=0;i<N;i++) if(care>>i&1){ if(!occurs(care&~(1u<<i),val&~(1u<<i))) {minimal="NO(site "+to_string(i/M)+","+to_string(i%M)+")";} }
   // modified burning on pattern sites
   vector<int> burnt(N,0); bool ch=true; while(ch){ch=false; for(int i=0;i<N;i++) if((care>>i&1)&&!burnt[i]){int x=i/M,y=i%M; int h=0,v=0;
      auto un=[&](int xx,int yy){ if(xx<0||xx>=L||yy<0||yy>=M) return 0; int j=xx*M+yy; return (int)((care>>j&1)&&!burnt[j]);};
      h=un(x-1,y)+un(x+1,y); v=un(x,y-1)+un(x,y+1); int z=val>>i&1; if(z>=min(h,v)){burnt[i]=1;ch=true;}}}
   bool allb=true; for(int i=0;i<N;i++) if((care>>i&1)&&!burnt[i]) allb=false;
   printf("Fig %s (%dx%d): %s in box; minimal-in-box: %s; modified burning: %s\n",name.c_str(),W,H, occ?"OCCURS (not forbidden)":"forbidden", minimal.c_str(), allb?"burns":"does not burn");
 }}
