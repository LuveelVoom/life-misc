// rulesearch: constrained search over range-1 Moore rules (totalistic / isotropic / MAP).
// Build: g++ -O2 -march=native -fopenmp -std=c++17 -o rulesearch rulesearch.cpp   (or: make)
//        clang++ -O2 -std=c++17 also works (single-threaded unless OpenMP is available)
//
// Conventions (checked against Golly gollybase/liferules.cpp):
//   512-entry table index: row-major from NW, NW = bit 8 (256) ... SE = bit 0, centre = bit 4 (16).
//   MAP: 86 base64 chars, 6 bits per char MSB first, bit k of the stream = table[k].
//   Hensel letters/neighbourhoods: Golly's entry tables; counts 5..8 are complements of 3..0.
//   Patterns: x to the right, y down (RLE order).
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <functional>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_thread_num(){ return 0; } static int omp_get_max_threads(){ return 1; }
#endif
#if defined(__ARM_NEON)
#include <arm_neon.h>
#define HAVE_NEON 1
#endif
using namespace std;
typedef uint64_t u64;
typedef vector<pair<int,int>> Cells;

static u64 splitmix(u64 &x){ u64 z=(x+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
[[noreturn]] static void die(const string&m){ fprintf(stderr,"error: %s\n",m.c_str()); exit(1); }

// ===================== rule notation =====================
static const char* LET[5]={"","ce","ceaikn","ceaiknjqry","ceaiknjqrytwz"};
static const int ENT1[2]={1,2}, ENT2[6]={5,10,3,40,33,68}, ENT3[10]={69,42,11,7,98,13,14,70,41,97},
                 ENT4[13]={325,170,15,45,99,71,106,102,43,101,105,78,108};
static const int* ENT[5]={nullptr,ENT1,ENT2,ENT3,ENT4};
static int flipB(int x){ return ((x&0x07)<<6)|((x&0x1c0)>>6)|(x&0x38); }
static int rotB(int x){ return ((x&0x4)<<6)|((x&0x20)<<2)|((x&0x100)>>2)|((x&0x2)<<4)|(x&0x10)|((x&0x80)>>4)|((x&0x1)<<2)|((x&0x8)>>2)|((x&0x40)>>6); }
static vector<int> d4orbit(int x){ set<int> s; int y=x; for(int i=0;i<4;i++){s.insert(y);y=rotB(y);} y=flipB(y); for(int i=0;i<4;i++){s.insert(y);y=rotB(y);} return vector<int>(s.begin(),s.end()); }

// Hensel class of every 512 index: hclass[i] = 0..101 (B: 0..50, S: 51..101), hname[c] e.g. "B2a"
static int hclass[512]; static string hname[102]; static int hcount[102]; static char hlet[102];
static void initHensel(){
    for(int i=0;i<512;i++) hclass[i]=-1;
    int id=0;
    for(int s=0;s<2;s++) for(int n=0;n<=8;n++){
        int m=min(n,8-n); int nl = (n==0||n==8)?1:(int)strlen(LET[m]);
        for(int li=0;li<nl;li++){
            hname[id]=string(s?"S":"B")+char('0'+n)+((n==0||n==8)?string():string(1,LET[m][li]));
            hcount[id]=n; hlet[id]=(n==0||n==8)?0:LET[m][li];
            if(n==0||n==8){ for(int x=0;x<512;x++) if(((x>>4)&1)==s && __builtin_popcount(x&~0x10)==n) hclass[x]=id; }
            else { int x = n<=4 ? ENT[n][li] : (ENT[8-n][li]^0x1ef); x=(x&~0x10)|(s<<4);
                   for(int y:d4orbit(x)){ if(hclass[y]!=-1) die("hensel table overlap"); hclass[y]=id; } }
            id++;
        }
    }
    for(int i=0;i<512;i++) if(hclass[i]<0) die("hensel table incomplete");
    if(id!=102) die("hensel count");
}
static const char* B64="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
struct Table { uint8_t t[512]; bool operator==(const Table&o)const{return !memcmp(t,o.t,512);} };

static Table parseRule(string r){
    Table T; memset(T.t,0,512);
    string s; for(char c:r) if(!isspace((unsigned char)c)) s+=c;
    if(s.size()>=3 && toupper(s[0])=='M'&&toupper(s[1])=='A'&&toupper(s[2])=='P'){
        string b=s.substr(3); while(!b.empty()&&b.back()=='=') b.pop_back();
        if(b.size()!=86) die("MAP rule must have 86 base64 characters (Moore): "+r);
        int k=0;
        for(size_t i=0;i<b.size()&&k<512;i++){ const char*p=strchr(B64,b[i]); if(!p) die("bad MAP char"); int c=p-B64;
            for(int bit=5;bit>=0&&k<512;bit--) T.t[k++]=(c>>bit)&1; }
        return T;
    }
    // B.../S... (Hensel, possibly totalistic). Accept "B3/S23", "b3s23", "B3a/S2-a3".
    size_t sp=string::npos; for(size_t i=1;i<s.size();i++) if(s[i]=='S'||s[i]=='s'){sp=i;break;}
    if(toupper(s[0])!='B'||sp==string::npos) die("rule must be B.../S... or MAP...: "+r);
    string parts[2]={s.substr(1,sp-1),s.substr(sp+1)}; if(!parts[0].empty()&&parts[0].back()=='/') parts[0].pop_back();
    for(int half=0;half<2;half++){
        const string&p=parts[half]; size_t i=0;
        while(i<p.size()){
            if(!isdigit((unsigned char)p[i])) die("bad rule near '"+p.substr(i)+"' in "+r);
            int n=p[i++]-'0'; if(n>8) die("count >8");
            bool neg=false; string lets;
            if(i<p.size()&&p[i]=='-'){neg=true;i++;}
            while(i<p.size()&&isalpha((unsigned char)p[i])) lets+=(char)tolower(p[i++]);
            for(int c=0;c<102;c++){
                if(c/51!=half||hcount[c]!=n) continue;
                bool in = lets.empty() ? true : ((lets.find(hlet[c])!=string::npos)!=neg);
                if(!lets.empty()){ for(char L:lets){ int m=min(n,8-n); if(!strchr(LET[m],L)) die(string("letter ")+L+" invalid for count "+char('0'+n)); } }
                if(in) for(int x=0;x<512;x++) if(hclass[x]==c) T.t[x]=1;
            }
        }
    }
    return T;
}
static bool isIsotropic(const Table&T){ for(int x=0;x<512;x++) for(int y:d4orbit(x)) if(T.t[x]!=T.t[y]) return false; return true; }
static string mapString(const Table&T){ string o="MAP"; for(int k=0;k<512;k+=6){ int c=0; for(int b=0;b<6;b++){ int i=k+b; c=(c<<1)|((i<512)?T.t[i]:0);} o+=B64[c]; } return o; }
static string ruleString(const Table&T){
    if(!isIsotropic(T)) return mapString(T);
    string out;
    for(int half=0;half<2;half++){
        out+= half?"/S":"B";
        for(int n=0;n<=8;n++){
            vector<int> all,have; for(int c=half*51;c<half*51+51;c++) if(hcount[c]==n){ all.push_back(c); int rep=-1; for(int x=0;x<512;x++) if(hclass[x]==c){rep=x;break;} if(T.t[rep]) have.push_back(c);}
            if(have.empty()) continue;
            if(have.size()==all.size()){ out+=char('0'+n); continue; }
            string pos,neg; for(int c:all){ bool h=find(have.begin(),have.end(),c)!=have.end(); (h?pos:neg)+=hlet[c]; }
            out+=char('0'+n); out+= (pos.size()<=neg.size()+1)? pos : "-"+neg;
        }
    }
    return out;
}

// ===================== symmetry classes (the search variables) =====================
static int SYM=1; // 0 totalistic, 1 isotropic, 2 none
static int NC; static int cls[512]; static vector<vector<int>> members; static vector<string> cname;
static void initClasses(){
    for(int i=0;i<512;i++) cls[i]=-1; members.clear(); cname.clear();
    for(int i=0;i<512;i++){
        if(cls[i]>=0) continue;
        vector<int> orb;
        if(SYM==2) orb={i};
        else if(SYM==1) orb=d4orbit(i);
        else { for(int y=0;y<512;y++) if(((y>>4)&1)==((i>>4)&1) && __builtin_popcount(y&~16)==__builtin_popcount(i&~16)) orb.push_back(y); }
        int id=members.size(); members.push_back(orb); for(int y:orb) cls[y]=id;
        if(SYM==1) cname.push_back(hname[hclass[i]]);
        else if(SYM==0) cname.push_back(string((i>>4)&1?"S":"B")+char('0'+__builtin_popcount(i&~16)));
        else cname.push_back("i"+to_string(i));
    }
    NC=members.size();
}
static bool respectsSym(const Table&T){ for(int c=0;c<NC;c++) for(int y:members[c]) if(T.t[y]!=T.t[members[c][0]]) return false; return true; }

// ===================== patterns =====================
static Cells parseRLE(const string&rle){
    Cells c; int x=0,y=0; string num;
    for(char ch:rle){
        if(isspace((unsigned char)ch)) continue;
        if(isdigit((unsigned char)ch)){num+=ch;continue;}
        int k=num.empty()?1:stoi(num); num.clear();
        if(ch=='b'||ch=='.'){x+=k;}
        else if(ch=='$'){y+=k;x=0;}
        else if(ch=='!') break;
        else if(isalpha((unsigned char)ch)){ for(int i=0;i<k;i++) c.push_back({x+i,y}); x+=k; } // any non-b letter = live
    }
    return c;
}
static void bbox(const Cells&c,int&mx,int&my,int&Mx,int&My){ mx=my=INT_MAX;Mx=My=INT_MIN; for(auto&p:c){mx=min(mx,p.first);my=min(my,p.second);Mx=max(Mx,p.first);My=max(My,p.second);} }
static string keyOf(const Cells&c){
    if(c.empty()) return string(2,'\0');
    int mx,my,Mx,My; bbox(c,mx,my,Mx,My); int w=Mx-mx+1,h=My-my+1;
    string k; k.push_back((char)(w&255)); k.push_back((char)(w>>8)); k.push_back((char)(h&255)); k.push_back((char)(h>>8));
    vector<uint8_t> b((w*h+7)/8,0); for(auto&p:c){int i=(p.second-my)*w+(p.first-mx); b[i>>3]|=1<<(i&7);}
    k.append(b.begin(),b.end()); return k;
}
static Cells keyCells(const string&k){ int w=(uint8_t)k[0]|((uint8_t)k[1]<<8),h=(uint8_t)k[2]|((uint8_t)k[3]<<8); Cells c; for(int y=0;y<h;y++) for(int x=0;x<w;x++){int i=y*w+x; if((uint8_t)k[4+(i>>3)]>>(i&7)&1) c.push_back({x,y});} return c; }
static string toRLE(const Cells&c0,const string&rule){
    Cells c=keyCells(keyOf(c0)); int w=0,h=0; for(auto&p:c){w=max(w,p.first+1);h=max(h,p.second+1);}
    set<pair<int,int>> S(c.begin(),c.end()); string body; int pend=0;
    vector<string> toks; auto run=[&](int n,char ch){ if(n<=0) return; toks.push_back((n>1?to_string(n):string())+ch); };
    for(int y=0;y<h;y++){
        int x=0; string row; vector<pair<int,char>> runs;
        while(x<w){ char ch=S.count({x,y})?'o':'b'; int n=0; while(x<w&&(S.count({x,y})?'o':'b')==ch){n++;x++;} runs.push_back({n,ch}); }
        if(!runs.empty()&&runs.back().second=='b') runs.pop_back();
        if(runs.empty()){ pend++; continue; }
        if(y>0){ run(pend+1,'$'); } pend=0;
        for(auto&r:runs) run(r.first,r.second);
    }
    toks.push_back("!");
    string out="x = "+to_string(w)+", y = "+to_string(h)+", rule = "+rule+"\n", line;
    for(auto&t:toks){ if(line.size()+t.size()>70){ out+=line+"\n"; line.clear(); } line+=t; }
    return out+line+"\n";
}
static Cells xform(const Cells&c,int o){ Cells r; r.reserve(c.size()); for(auto&p:c){int x=p.first,y=p.second; if(o&4) swap(x,y); if(o&1) x=-x; if(o&2) y=-y; r.push_back({x,y});} return r; }
static vector<int> groupOrients(){ if(SYM==2) return {0}; return {0,1,2,3,4,5,6,7}; }

static Cells stepCells(const Cells&c,const uint8_t*T,int8_t*used=nullptr){ // used[512]: table entries looked up
    Cells out; if(c.empty()) return out;
    int mx,my,Mx,My; bbox(c,mx,my,Mx,My);
    int W=Mx-mx+5,H=My-my+5; vector<uint8_t> g(W*H,0);
    for(auto&p:c) g[(p.second-my+2)*W+(p.first-mx+2)]=1;
    for(int y=1;y<H-1;y++) for(int x=1;x<W-1;x++){
        int i=(g[(y-1)*W+x-1]<<8)|(g[(y-1)*W+x]<<7)|(g[(y-1)*W+x+1]<<6)|(g[y*W+x-1]<<5)|(g[y*W+x]<<4)|(g[y*W+x+1]<<3)|(g[(y+1)*W+x-1]<<2)|(g[(y+1)*W+x]<<1)|g[(y+1)*W+x+1];
        if(used) used[i]=T[i];
        if(T[i]) out.push_back({x+mx-2,y+my-2});
    }
    return out;
}
// evolve until normalized state repeats (translation allowed) or dies. returns gens to repeat t and earlier index s (s=-1: died at t)
struct Evo { vector<Cells> gens; int s=-2,t=0; int dx=0,dy=0; }; // s: -1 died at t, -2 no repeat, -3 grew past maxpop at t
static Evo evolve(const Cells&c0,const uint8_t*T,int maxg,int8_t*used,int maxpop=INT_MAX){
    Evo E; map<string,int> seen; Cells c=c0;
    for(int g=0;g<=maxg;g++){
        E.gens.push_back(c);
        if(c.empty()&&g>0){ E.s=-1;E.t=g; return E; }
        if((int)c.size()>maxpop){ E.s=-3; E.t=g; return E; }
        string k=keyOf(c); auto it=seen.find(k);
        if(it!=seen.end()){ E.s=it->second;E.t=g; int a,b,A,B,a2,b2; bbox(E.gens[E.s],a,b,A,B); bbox(c,a2,b2,A,B); E.dx=a2-a; E.dy=b2-b; return E; }
        seen[k]=g;
        if(g<maxg) c=stepCells(c,T,used);
    }
    return E; // s=-2: no repeat
}
static vector<Cells> splitComponents(const Cells&c){ // clusters at Chebyshev distance <=2 joined
    set<pair<int,int>> S(c.begin(),c.end()); vector<Cells> out;
    while(!S.empty()){ Cells comp; vector<pair<int,int>> st={*S.begin()}; S.erase(S.begin());
        while(!st.empty()){ auto p=st.back(); st.pop_back(); comp.push_back(p);
            for(int dy=-2;dy<=2;dy++) for(int dx=-2;dx<=2;dx++){ auto it=S.find({p.first+dx,p.second+dy}); if(it!=S.end()){st.push_back(*it);S.erase(it);} } }
        out.push_back(comp); }
    sort(out.begin(),out.end(),[](const Cells&a,const Cells&b){int ax,ay,AX,AY,bx,by,BX,BY; bbox(a,ax,ay,AX,AY); bbox(b,bx,by,BX,BY); return make_pair(ax,ay)<make_pair(bx,by);});
    return out;
}

// ===================== classification =====================
struct Info { int type=0; int period=0,dx=0,dy=0; Cells env; }; // 1 osc/still, 2 ship, 3 dies, 4 unknown
static Info classify(const Cells&c0,const uint8_t*T,int maxp){
    Info I; Cells c=c0; string k0=keyOf(c0); int mx0,my0,a,b; bbox(c0,mx0,my0,a,b);
    set<pair<int,int>> env(c0.begin(),c0.end());
    for(int g=1;g<=maxp;g++){
        c=stepCells(c,T);
        if(c.empty()){I.type=3;return I;}
        if((int)c.size()>1000){I.type=4;return I;}
        for(auto&p:c) env.insert(p);
        if(keyOf(c)==k0){ int mx,my; bbox(c,mx,my,a,b); I.period=g; I.dx=mx-mx0; I.dy=my-my0; I.type=(I.dx||I.dy)?2:1;
            for(auto&p:env) I.env.push_back({p.first-mx0,p.second-my0}); return I; }
    }
    I.type=4; return I;
}
// pseudo-object: some partition of its 8-connected components evolves independently for one period
static vector<Cells> comps8(const Cells&c){
    set<pair<int,int>> S(c.begin(),c.end()); vector<Cells> out;
    while(!S.empty()){ Cells comp; vector<pair<int,int>> st={*S.begin()}; S.erase(S.begin());
        while(!st.empty()){ auto p=st.back(); st.pop_back(); comp.push_back(p);
            for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){ auto it=S.find({p.first+dx,p.second+dy}); if(it!=S.end()){st.push_back(*it);S.erase(it);} } }
        out.push_back(comp); }
    return out;
}
static bool isPseudo(const Cells&c,const uint8_t*T,int P){
    auto cs=comps8(c); int k=cs.size(); if(k<2||k>7) return false;
    // enumerate set partitions via restricted growth strings
    vector<int> a(k,0),mx(k,0);
    while(true){
        int ng=*max_element(a.begin(),a.end())+1;
        if(ng>=2){
            vector<Cells> gs(ng); for(int i=0;i<k;i++) gs[a[i]].insert(gs[a[i]].end(),cs[i].begin(),cs[i].end());
            Cells whole=c; bool ok=true; vector<string> k0; for(auto&G:gs) k0.push_back(keyOf(G));
            for(int g=0;g<P&&ok;g++){
                whole=stepCells(whole,T); set<pair<int,int>> W(whole.begin(),whole.end()); size_t tot=0; set<pair<int,int>> U;
                for(auto&G:gs){ G=stepCells(G,T); tot+=G.size(); U.insert(G.begin(),G.end()); }
                if(U!=W||tot!=W.size()) ok=false;
            }
            // each group must itself be a periodic object (returns to its own shape), not a dying spark
            for(size_t q=0;q<gs.size()&&ok;q++) if(gs[q].empty()||keyOf(gs[q])!=k0[q]) ok=false;
            if(ok) return true;
        }
        int i=k-1; while(i>0 && a[i]==mx[i-1]+1) i--; // next RGS
        if(i==0) break;
        a[i]++; for(int j=i+1;j<k;j++) a[j]=0;
        for(int j=1;j<k;j++) mx[j]=max(mx[j-1],a[j]);
    }
    return false;
}

// ===================== config =====================
struct Pat { string name; Cells cells; Evo evo; };
struct Target { string name; int pat; double w=1; double minRate=0; set<string> keys; int period=0,dx=0,dy=0; };
struct Cfg {
    string refRule; Table ref; vector<Pat> pats; vector<Target> tg; vector<pair<string,int>> fixes;
    int soupW=16,soupH=16; double dens=0.5; int N=1000; int grid=256; int gmax=2000; int pmax=2500; int maxp=256; int cap=4;
    double maxBad=1.0, costW=0, costRef=0, maxCost=0, minOscPeriod=3; int lahc=50; int reseed=200; double soupCostLimit=-1; u64 seed=1;
    string logFile="discoveries.txt", archiveFile="archive.txt", bestFile; vector<string> starts; int maxgCons=4000;
    map<string,double> reward;
    string soupSym="C1"; int gunShips=10; bool growth=true;
    double varyMin=0.002; int varyProbe=1; int varyStrict=0; double varyDecay=1.0, varyCommonShip=0, varyCommonRate=0.1, varyConfined=1.0;
    int probe=0, probeGrid=512, probeShips=4, probeStills=16, probeGens=600, wickK=6, wickCells=8, wickMax=400, fossilCells=12, fossilMax=12; double fossilRate=0.002; int repSeeds=1000, repSoups=200, repMax=64, repGens=192, varyRep=0; vector<string> repAdd; // explore: novelty-credit multipliers by kind (ship, oblique, diagonal, orthogonal, osc)
} C;
static vector<int> forcedV; // per class: -1 free, 0/1 forced

// pattern RLE given inline, or as @FILE (an .rle file; '#' comment lines and the 'x = ...' header are skipped)
static string patternText(const string&r){
    if(r.empty()||r[0]!='@') return r;
    ifstream f(r.substr(1)); if(!f) die("cannot open pattern file "+r.substr(1)); string line,out;
    while(getline(f,line)){ size_t i=line.find_first_not_of(" \t\r"); if(i==string::npos||line[i]=='#'||line[i]=='x') continue; out+=line; }
    return out;
}
// overrides: extra "key value" lines from the command line (key=value, commas for spaces), applied after the file
static void loadCfg(const string&fn,const vector<string>&overrides){
    ifstream f(fn); if(!f) die("cannot open config "+fn); string line; int ln=0;
    vector<tuple<string,string,bool>> pats; vector<tuple<string,double>> tgs; map<string,double> minr;
    vector<string> lines; while(getline(f,line)) lines.push_back(line);
    int nfile=lines.size(); for(auto&o:overrides) lines.push_back(o);
    for(size_t li=0;li<lines.size();li++){ line=lines[li];
        ln++; auto h=line.find('#'); if(h!=string::npos) line=line.substr(0,h);
        istringstream is(line); string k; if(!(is>>k)) continue;
        auto need=[&](bool ok){ if(!ok) die((ln>nfile?string("command-line override: "):"config line "+to_string(ln)+": ")+line); };
        if(k=="rule"){ need((bool)(is>>C.refRule)); }
        else if(k=="symmetry"){ string s; is>>s; SYM= s=="totalistic"?0: s=="isotropic"?1: s=="none"?2: -1; need(SYM>=0); }
        else if(k=="pattern"){ string n,r,opt; need((bool)(is>>n>>r)); is>>opt; pats.push_back({n,patternText(r),opt=="split"}); }
        else if(k=="target"){ string n; double w=1; need((bool)(is>>n)); is>>w; tgs.push_back({n,w}); }
        else if(k=="min_rate"){ string n; double r; need((bool)(is>>n>>r)); minr[n]=r; }
        else if(k=="fix"){ string n; int v; need((bool)(is>>n>>v)); C.fixes.push_back({n,v}); }
        else if(k=="soup"){ need((bool)(is>>C.soupW>>C.soupH>>C.dens)); }
        else if(k=="soups") is>>C.N; else if(k=="grid") is>>C.grid; else if(k=="gmax") is>>C.gmax; else if(k=="pmax") is>>C.pmax;
        else if(k=="maxp") is>>C.maxp; else if(k=="cap") is>>C.cap; else if(k=="max_bad") is>>C.maxBad;
        else if(k=="cost_weight") is>>C.costW; else if(k=="cost_ref") is>>C.costRef; else if(k=="max_cost") is>>C.maxCost;
        else if(k=="min_osc_period") is>>C.minOscPeriod; else if(k=="lahc") is>>C.lahc; else if(k=="reseed") is>>C.reseed;
        else if(k=="seed") is>>C.seed;
        else if(k=="log"||k=="archive"){ string v; is>>v; if(v=="none") v.clear(); (k=="log"?C.logFile:C.archiveFile)=v; } // empty or none: off
        else if(k=="start"){ string r; is>>r; C.starts.push_back(r); }
        else if(k=="constraint_gens") is>>C.maxgCons;
        else if(k=="soup_cost_limit") is>>C.soupCostLimit;
        else if(k=="best") is>>C.bestFile;
        else if(k=="soup_symmetry"){ need((bool)(is>>C.soupSym)); need(C.soupSym=="C1"||C.soupSym=="C2"||C.soupSym=="C4"||C.soupSym=="D2"||C.soupSym=="D4"||C.soupSym=="D8"); }
        else if(k=="gun_ships") is>>C.gunShips;
        else if(k=="vary_min") is>>C.varyMin; else if(k=="vary_probe") is>>C.varyProbe;
        else if(k=="vary_strict") is>>C.varyStrict; else if(k=="vary_decay") is>>C.varyDecay;
        else if(k=="vary_confined") is>>C.varyConfined;
        else if(k=="vary_common_ship"){ need((bool)(is>>C.varyCommonShip)); is>>C.varyCommonRate; }
        else if(k=="probe") is>>C.probe; else if(k=="probe_grid") is>>C.probeGrid; else if(k=="probe_ships") is>>C.probeShips; else if(k=="probe_stills") is>>C.probeStills;
        else if(k=="replicator_add"){ string r; need((bool)(is>>r)); C.repAdd.push_back(r); }
        else if(k=="vary_replicators") is>>C.varyRep; else if(k=="replicator_seeds") is>>C.repSeeds; else if(k=="replicator_soups") is>>C.repSoups; else if(k=="replicator_max") is>>C.repMax; else if(k=="replicator_gens") is>>C.repGens;
        else if(k=="probe_gens") is>>C.probeGens; else if(k=="wick_period") is>>C.wickK; else if(k=="wick_cells") is>>C.wickCells; else if(k=="wick_max") is>>C.wickMax;
        else if(k=="fossil_cells") is>>C.fossilCells; else if(k=="fossil_rate") is>>C.fossilRate; else if(k=="fossil_max") is>>C.fossilMax;
        else if(k=="growth"){ int v; is>>v; C.growth=v; }
        else if(k=="reward"){ string n; double w; need((bool)(is>>n>>w)); need(n=="ship"||n=="oblique"||n=="diagonal"||n=="orthogonal"||n=="osc"||n=="growth"||n=="gun"||n=="catalyst"||n=="signal"||n=="active"); C.reward[n]=w; }
        else die("unknown config key '"+k+"'"+(ln>nfile?string(" (command line)"):" at line "+to_string(ln)));
    }
    if(C.refRule.empty()) die("config needs 'rule'");
    initClasses(); C.ref=parseRule(C.refRule);
    if(!respectsSym(C.ref)) die("reference rule is not symmetric under the chosen symmetry");
    if(C.ref.t[0]) die("B0 reference rules are not supported");
    for(auto&[n,r,sp]:pats){
        Cells c=parseRLE(r); if(c.empty()) die("empty pattern "+n);
        if(!sp){ C.pats.push_back({n,c,{}}); continue; }
        auto comps=splitComponents(c);
        for(size_t i=0;i<comps.size();i++) C.pats.push_back({n+"."+to_string(i),comps[i],{}});
        // check components don't interact over the longest constraint horizon (checked after evolution below)
    }
    // constraints
    forcedV.assign(NC,-1);
    auto force=[&](int idx,int v,const string&why){ int c=cls[idx]; if(forcedV[c]>=0&&forcedV[c]!=v) die("conflicting constraint on "+cname[c]+" ("+why+")"); forcedV[c]=v; };
    force(0,0,"B0 off");
    for(auto&P:C.pats){
        int8_t used[512]; memset(used,-1,512); P.evo=evolve(P.cells,C.ref.t,C.maxgCons,used,C.pmax);
        if(P.evo.s==-3) die("pattern "+P.name+" grows past "+to_string(C.pmax)+" cells (gen "+to_string(P.evo.t)+") in the reference rule");
        if(P.evo.s==-2) die("pattern "+P.name+" does not repeat within "+to_string(C.maxgCons)+" gens; raise constraint_gens or split it");
        for(int i=0;i<512;i++) if(used[i]>=0) force(i,used[i],P.name);
    }
    // split sanity: whole pattern vs union of components
    for(auto&[n,r,sp]:pats){ if(!sp) continue; Cells whole=parseRLE(r); int H=0; vector<Cells> comps; for(auto&P:C.pats) if(P.name.rfind(n+".",0)==0){ comps.push_back(P.cells); H=max(H,P.evo.t);}
        for(int g=0;g<H;g++){ whole=stepCells(whole,C.ref.t); set<pair<int,int>> U; for(auto&cc:comps){cc=stepCells(cc,C.ref.t); U.insert(cc.begin(),cc.end());} if(U!=set<pair<int,int>>(whole.begin(),whole.end())) die("components of "+n+" interact at gen "+to_string(g+1)+"; don't split it"); } }
    for(auto&[n,v]:C.fixes){
        bool done=false;
        if(SYM!=2){ for(int c=0;c<NC;c++) if(cname[c]==n){ force(members[c][0],v,"fix"); done=true; } }
        if(!done && SYM==1){ for(int i=0;i<512;i++) if(hname[hclass[i]]==n){ force(i,v,"fix"); done=true; } }
        if(!done && n[0]=='i'){ force(stoi(n.substr(1)),v,"fix"); done=true; }
        if(!done) die("fix: unknown transition "+n+" (isotropic: B2a, S3; totalistic: B3; none: i<index>)");
    }
    for(auto&[n,w]:tgs){
        int pi=-1; for(size_t i=0;i<C.pats.size();i++) if(C.pats[i].name==n) pi=i;
        if(pi<0) die("target "+n+": no such pattern (run 'check' to see split names)");
        Target T; T.name=n; T.pat=pi; T.w=w; T.minRate=minr.count(n)?minr[n]:0;
        auto&E=C.pats[pi].evo; if(E.s<0) die("target "+n+" dies");
        T.period=E.t-E.s; T.dx=E.dx; T.dy=E.dy;
        for(int o:groupOrients()) for(int g=E.s;g<E.t;g++) T.keys.insert(keyOf(xform(E.gens[g],o)));
        C.tg.push_back(T);
    }
    for(auto&kv:minr){ bool ok=false; for(auto&t:C.tg) ok|=t.name==kv.first; if(!ok) die("min_rate for non-target "+kv.first); }
    if(C.tg.size()>8) die("at most 8 targets");
    if(C.soupCostLimit<0) C.soupCostLimit = C.maxCost>0 ? min(20.0,10*C.maxCost) : 20.0;
}
static Table project(Table T){ for(int c=0;c<NC;c++){ int v= forcedV[c]>=0?forcedV[c]:T.t[members[c][0]]; for(int y:members[c]) T.t[y]=v; } return T; }
static bool satisfies(const Table&T){ for(int c=0;c<NC;c++) if(forcedV[c]>=0) for(int y:members[c]) if(T.t[y]!=forcedV[c]) return false; return true; }

// ===================== soup engine =====================
// The sweep builds neighbourhood indices column by column ("T2 order"); P2[j] is the standard (MAP-order) index.
static int P2[512];
static void initP2(){ for(int j=0;j<512;j++){ int m=0; for(int c=0;c<3;c++) for(int r=0;r<3;r++){ int col=(j>>(6-3*c))&7; if((col>>(2-r))&1) m|=1<<(8-(r*3+c)); } P2[j]=m; } }
// soup cells in [0,fw)x[0,fh): a soupW x soupH random block, completed by the configured symmetry
// (C2: 180-degree turn, 2W x H... see README). C1 reproduces the original soups bit for bit.
static Cells makeSoup(u64 seed,int&fw,int&fh){
    Cells c; u64 s=seed; int W=C.soupW,H=C.soupH;
    if(W<=64 && C.dens==0.5){ for(int y=0;y<H;y++){u64 r=splitmix(s); for(int x=0;x<W;x++) if(r>>x&1) c.push_back({x,y});} }
    else { u64 thr=(u64)(C.dens*18446744073709551615.0); for(int y=0;y<H;y++) for(int x=0;x<W;x++) if(splitmix(s)<thr) c.push_back({x,y}); }
    const string&S=C.soupSym; fw=W; fh=H; if(S=="C1") return c;
    Cells o; auto add=[&](int x,int y){ o.push_back({x,y}); };
    if(S=="D2"){ fw=2*W; for(auto&p:c){ add(p.first,p.second); add(2*W-1-p.first,p.second); } }
    else if(S=="C2"){ fh=2*H; for(auto&p:c){ add(p.first,p.second); add(W-1-p.first,2*H-1-p.second); } }
    else { int n=max(W,H); fw=fh=2*n; // C4, D4, D8 use a square quadrant
        for(auto&p:c){ int x=p.first,y=p.second;
            if(S=="C4"){ add(x,y); add(2*n-1-y,x); add(2*n-1-x,2*n-1-y); add(y,2*n-1-x); }
            else { int xs[2]={x,2*n-1-x}, ys[2]={y,2*n-1-y}; for(int i=0;i<2;i++) for(int j=0;j<2;j++){ add(xs[i],ys[j]); if(S=="D8") add(ys[j],xs[i]); } } } }
    sort(o.begin(),o.end()); o.erase(unique(o.begin(),o.end()),o.end()); return o;
}
#ifdef HAVE_NEON
// 16 cells at once: the 512-entry table is packed as 64 bytes of bits; the top 6 neighbourhood bits (NW..E) pick the
// byte (one tbl lookup), the bottom 3 (SW,S,SE) the bit. r0,r1,r2 = rows y-1,y,y+1; returns 0/1 bytes for x..x+15.
static inline uint8x16_t step16(const uint8_t*r0,const uint8_t*r1,const uint8_t*r2,int x,const uint8x16x4_t&tb){
    uint8x16_t hi=vorrq_u8(vorrq_u8(vshlq_n_u8(vld1q_u8(r0+x-1),5),vshlq_n_u8(vld1q_u8(r0+x),4)),
                  vorrq_u8(vorrq_u8(vshlq_n_u8(vld1q_u8(r0+x+1),3),vshlq_n_u8(vld1q_u8(r1+x-1),2)),vorrq_u8(vshlq_n_u8(vld1q_u8(r1+x),1),vld1q_u8(r1+x+1))));
    uint8x16_t lo=vorrq_u8(vorrq_u8(vshlq_n_u8(vld1q_u8(r2+x-1),2),vshlq_n_u8(vld1q_u8(r2+x),1)),vld1q_u8(r2+x+1));
    return vandq_u8(vshlq_u8(vqtbl4q_u8(tb,hi),vnegq_s8(vreinterpretq_s8_u8(lo))),vdupq_n_u8(1));
}
// one bit per nonzero lane, lane k at bit 4k
static inline u64 laneBits(uint8x16_t v){ return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(vtstq_u8(v,v)),4)),0)&0x1111111111111111ULL; }
#endif
struct Cand { long count=0; Cells ex; Info I; u64 seed=0; };
struct SoupRes { uint8_t cnt[8]; uint8_t bad; double cost; uint16_t growP; int16_t growD; uint16_t gunN; uint8_t gunKey; uint16_t activeN; }; // activeN: bounded region that never settled // growth: pop +growD every growP gens; gun: gunN copies of one ship
// everything one soup produces (kept per soup so merging is in soup order, independent of thread count)
struct SoupOut { SoupRes r; vector<pair<string,Cand>> cands; vector<uint32_t> usage; Cand gunShip; vector<pair<string,Cand>> ash; };
struct ShipEvt { Cells c; Info I; int t; }; // a ship deleted by the scan (probe runs)
struct Engine {
    const uint8_t*T=nullptr; uint8_t T2[512]; int W,H; vector<uint8_t> a,b; vector<u64> Z; vector<int> stamp; int st=0;
    unordered_map<string,Info> cache, cache2; unordered_map<string,Cand> cands; u64 curSeed=0; bool countUsage=false; uint32_t*usage=nullptr;
    static const int LG=256; vector<uint8_t> ca,cb; vector<int> hx,hy,hoff,hbx,hby; vector<u64> hh; // classification scratch
    uint8_t TB[64];
#ifdef HAVE_NEON
    uint8x16x4_t TBv;
#endif
    void setRule(const uint8_t*t){ T=t; cache.clear(); cache2.clear(); cands.clear(); for(int j=0;j<512;j++) T2[j]=t[P2[j]];
        for(int k=0;k<64;k++){ TB[k]=0; for(int b=0;b<8;b++) TB[k]|=t[k*8+b]<<b; }
#ifdef HAVE_NEON
        TBv=vld1q_u8_x4(TB);
#endif
    }
    Engine(int g):W(g),H(g),a(g*g,0),b(g*g,0),Z(g*g),stamp(g*g,0),ca(LG*LG,0),cb(LG*LG,0){ u64 s=12345; for(auto&z:Z) z=splitmix(s); }

    // Same Info as classify(c0,T,maxp), on a local byte grid. Additionally stops as soon as the evolution enters a
    // cycle that does not contain c0 (c0 can then never recur, so the answer is type 4).
    // ahead>=0 (only with maxp=64): may run 'ahead' generations past maxp, and returns in 'fut' the generations m
    // (multiples of 8) whose state's own classify(.,64) result is proven by this history to be 3 or 4 (type in fut).
    // A cluster that stays isolated is seen in exactly that state by the scan 8 generations later.
    long nskip=0;
    Info classifyFast(const Cells&c0,int maxp,vector<pair<int,int>>*fut=nullptr){
        Info I; int mx0,my0,Mx0,My0; bbox(c0,mx0,my0,Mx0,My0);
        int ox=LG/2-(Mx0-mx0)/2, oy=LG/2-(My0-my0)/2;
        hx.clear(); hy.clear(); hoff.clear(); hbx.clear(); hby.clear(); hh.clear();
        { vector<pair<int,int>> s; s.reserve(c0.size()); for(auto&p:c0) s.push_back({p.second-my0+oy,p.first-mx0+ox}); sort(s.begin(),s.end());
          hoff.push_back(0); for(auto&q:s){ hy.push_back(q.first); hx.push_back(q.second); ca[q.first*LG+q.second]=1; } }
        auto hashRange=[&](int k0,int k1,int x0,int y0){ u64 h=0x9E3779B97F4A7C15ULL^(u64)(k1-k0); for(int k=k0;k<k1;k++) h=(h^(u64)(((hy[k]-y0)<<10)|(hx[k]-x0)))*0x100000001B3ULL; return h; };
        auto endOf=[&](int g){ return g+1<(int)hoff.size()? hoff[g+1] : (int)hx.size(); };
        auto same=[&](int g1,int g2){ int a0=hoff[g1],a1=endOf(g1),b0=hoff[g2],b1=endOf(g2); if(a1-a0!=b1-b0) return false;
            for(int k=0;k<a1-a0;k++) if(hx[a0+k]-hbx[g1]!=hx[b0+k]-hbx[g2]||hy[a0+k]-hby[g1]!=hy[b0+k]-hby[g2]) return false; return true; };
        int bx0=ox,bx1=Mx0-mx0+ox,by0=oy,by1=My0-my0+oy, pbx0=bx0,pbx1=bx1,pby0=by0,pby1=by1;
        int ux0=bx0,ux1=bx1,uy0=by0,uy1=by1;
        hh.push_back(hashRange(0,hx.size(),bx0,by0)); hbx.push_back(bx0); hby.push_back(by0);
        uint8_t*A=ca.data(),*B=cb.data(); bool fallback=false; I.type=4;
        int glim = maxp;
        int E=0, why=0, cyc=-1; // last generation in history; why: 1 died, 2 too big, 3 cycle back to gen cyc, 4 limit
        for(int g=1;;g++){
            if(g>glim){ E=g-1; why=4; break; }
            int sx0=min(bx0,pbx0)-1, sx1=max(bx1,pbx1)+1, sy0=min(by0,pby0)-1, sy1=max(by1,pby1)+1;
            if(sx0<1||sy0<1||sx1>LG-2||sy1>LG-2){ if(g<=maxp) fallback=true; E=g-1; why=4; break; }
            ux0=min(ux0,sx0); ux1=max(ux1,sx1); uy0=min(uy0,sy0); uy1=max(uy1,sy1);
            pbx0=bx0;pbx1=bx1;pby0=by0;pby1=by1;
            int start=hx.size(); hoff.push_back(start); int nx0=INT_MAX,nx1=-1,ny0=INT_MAX,ny1=-1;
            for(int y=sy0;y<=sy1;y++){
                const uint8_t*r0=A+(y-1)*LG,*r1=A+y*LG,*r2=A+(y+1)*LG; uint8_t*o=B+y*LG;
                int x=sx0;
#ifdef HAVE_NEON
                // the last chunk may run past sx1 while inside the row: those lanes have empty neighbourhoods and get 0
                for(;x<=sx1&&(x+15<=sx1||x+16<=LG-1);x+=16){ uint8x16_t v=step16(r0,r1,r2,x,TBv); vst1q_u8(o+x,v); u64 m=laneBits(v);
                    while(m){ int xx=x+(__builtin_ctzll(m)>>2); m&=m-1; hx.push_back(xx); hy.push_back(y); if(xx<nx0)nx0=xx; if(xx>nx1)nx1=xx; if(y<ny0)ny0=y; ny1=y; } }
#endif
                int j= x<=sx1 ? (((r0[x-1]<<2)|(r1[x-1]<<1)|r2[x-1])<<3)|(r0[x]<<2)|(r1[x]<<1)|r2[x] : 0;
                for(;x<=sx1;x++){ j=((j<<3)&0x1ff)|(r0[x+1]<<2)|(r1[x+1]<<1)|r2[x+1];
                    uint8_t v=T2[j]; o[x]=v; if(v){ hx.push_back(x); hy.push_back(y); if(x<nx0)nx0=x; if(x>nx1)nx1=x; if(y<ny0)ny0=y; ny1=y; } }
            }
            swap(A,B);
            int n=hx.size()-start;
            if(n==0){ if(g<=maxp) I.type=3; E=g; why=1; break; }
            if(n>1000){ E=g; why=2; break; }
            bx0=nx0;bx1=nx1;by0=ny0;by1=ny1;
            u64 hg=hashRange(start,hx.size(),nx0,ny0); hh.push_back(hg); hbx.push_back(nx0); hby.push_back(ny0);
            if(hg==hh[0]&&same(0,g)){
                if(g<=maxp){ I.period=g; I.dx=nx0-hbx[0]; I.dy=ny0-hby[0]; I.type=(I.dx||I.dy)?2:1;
                    vector<pair<int,int>> env(hx.size()); for(size_t k=0;k<hx.size();k++) env[k]={hx[k]-ox,hy[k]-oy};
                    sort(env.begin(),env.end()); env.erase(unique(env.begin(),env.end()),env.end()); I.env=std::move(env); fut=nullptr; }
                E=g; why=3; cyc=0; break;
            }
            for(int k=1;k<g&&cyc<0;k++) if(hh[k]==hg&&same(k,g)) cyc=k;
            if(cyc>=0){ E=g; why=3; break; }
        }
        for(int y=uy0-1;y<=uy1+1;y++){ memset(&ca[y*LG+ux0-1],0,ux1-ux0+3); memset(&cb[y*LG+ux0-1],0,ux1-ux0+3); }
        if(fallback){ I=classify(c0,T,maxp); return I; }
        if(fut) for(int m=8;m<E;m+=8){
            // s_m never recurs unless it lies on the detected cycle; it dies iff why==1
            if(why==1){ fut->push_back({m,E-m<=64?3:4}); continue; }
            if(why==2){ fut->push_back({m,4}); continue; }
            if(why==3){ if(m<cyc || E-cyc>64) fut->push_back({m,4}); continue; }
            if(E-m>=64) fut->push_back({m,4});
        }
        return I;
    }
    Info& cls1(const Cells&c){ string k=keyOf(c); auto it=cache.find(k); if(it!=cache.end()) return it->second; if(cache.size()>200000) cache.clear(); return cache[k]=classifyFast(c,64); }
    Info& cls2(const Cells&c){ string k=keyOf(c); auto it=cache2.find(k); if(it!=cache2.end()) return it->second; if(cache2.size()>50000) cache2.clear(); return cache2[k]=classifyFast(c,C.maxp); }
    bool collectAsh=false; unordered_map<string,Cand> ash; vector<ShipEvt>*shipLog=nullptr; // ash: every isolated final still life / oscillator
    void note(const Cells&c,const Info&I){ string k=keyOf(c); auto&e=cands[k]; if(!e.count){e.ex=c;e.I=I;e.seed=curSeed;} e.count++; }
    void run(u64 seed, SoupOut&out){
        cands.clear();
        if(countUsage){ out.usage.assign(512,0); usage=out.usage.data(); } else { out.usage.clear(); usage=nullptr; }
        ash.clear(); run0(seed,out.r);
        out.cands.assign(cands.begin(),cands.end()); cands.clear(); out.gunShip=std::move(gunShip);
        out.ash.assign(ash.begin(),ash.end()); ash.clear();
    }
    vector<int> popH; map<string,int> shipN; Cand gunShip; int lastShipT=-1000000, curT=0; string lastShipK;
    void run0(u64 seed, SoupRes&R){
        memset(&R,0,sizeof R); curSeed=seed; shipN.clear(); gunShip=Cand(); lastShipT=-1000000;
        fill(a.begin(),a.end(),0); fill(b.begin(),b.end(),0);
        u64 h=0; int fw,fh; Cells soup=makeSoup(seed,fw,fh); int x0=W/2-fw/2,y0=H/2-fh/2;
        for(auto&p:soup){ a[(y0+p.second)*W+x0+p.first]=1; h^=Z[(y0+p.second)*W+x0+p.first]; }
        int bx0=x0,bx1=x0+fw-1,by0=y0,by1=y0+fh-1; int pbx0=bx0,pbx1=bx1,pby0=by0,pby1=by1;
        u64 hist[256]; int t=0; bool stable=false, edgeHit=false; long lastpop=1; double cost=0;
        popH.assign(C.gmax+1,0);
        for(t=0;t<C.gmax;t++){
            hist[t&255]=h;
            if(t>=20){ for(int p=1;p<=60&&2*p<=t;p++) if(hist[(t-p)&255]==h&&hist[(t-2*p)&255]==h){stable=true;
                    // the state recurred although a ship was deleted inside the window: something keeps emitting it
                    if(C.growth&&lastShipT>t-2*p&&shipN[lastShipK]>=2){ R.gunN=max(1,shipN[lastShipK]); gunShip.ex=keyCells(lastShipK); gunShip.I=cls1(gunShip.ex); gunShip.seed=seed; gunShip.count=R.gunN; }
                    break;} if(stable) break; }
            curT=t;
            if(t%8==0&&t>0&&lastpop<800) scan(R,false,bx0,bx1,by0,by1,h);
            int sx0=max(1,min(bx0,pbx0)-1), sx1=min(W-2,max(bx1,pbx1)+1), sy0=max(1,min(by0,pby0)-1), sy1=min(H-2,max(by1,pby1)+1);
            pbx0=bx0;pbx1=bx1;pby0=by0;pby1=by1;
            int nx0=INT_MAX,nx1=-1,ny0=INT_MAX,ny1=-1; long pop=0;
            cost+=(double)(sx1-sx0+1)*(sy1-sy0+1);
            if(cost>C.soupCostLimit*1e6){ R.bad=3; R.cost=cost; detectGrowth(R,t); return; }
            for(int y=sy0;y<=sy1;y++){
                const uint8_t*r0=&a[(y-1)*W],*r1=&a[y*W],*r2=&a[(y+1)*W]; uint8_t*o=&b[y*W];
                int x=sx0-1; int j=(r0[x]<<2)|(r1[x]<<1)|r2[x]; x++; j=(j<<3)|(r0[x]<<2)|(r1[x]<<1)|r2[x];
                if(usage){ for(x=sx0;x<=sx1;x++){ j=((j<<3)&0x1ff)|(r0[x+1]<<2)|(r1[x+1]<<1)|r2[x+1]; usage[j]++;
                        uint8_t v=T2[j]; o[x]=v; if(v!=r1[x]) h^=Z[y*W+x]; if(v){pop++; nx0=min(nx0,x); nx1=max(nx1,x); ny0=min(ny0,y); ny1=max(ny1,y);} } }
                else { x=sx0;
#ifdef HAVE_NEON
                    for(;x<=sx1&&(x+15<=sx1||x+16<=W-1);x+=16){ uint8x16_t v=step16(r0,r1,r2,x,TBv); vst1q_u8(o+x,v); // may overrun sx1, see classifyFast
                        u64 d=laneBits(veorq_u8(v,vld1q_u8(r1+x))); const u64*z=&Z[y*W+x]; while(d){ h^=z[__builtin_ctzll(d)>>2]; d&=d-1; }
                        u64 m=laneBits(v); if(m){ pop+=__builtin_popcountll(m); int f=x+(__builtin_ctzll(m)>>2), l=x+((63-__builtin_clzll(m))>>2);
                            if(f<nx0)nx0=f; if(l>nx1)nx1=l; if(y<ny0)ny0=y; if(y>ny1)ny1=y; } }
                    if(x<=sx1) j=(((r0[x-1]<<2)|(r1[x-1]<<1)|r2[x-1])<<3)|(r0[x]<<2)|(r1[x]<<1)|r2[x];
#endif
                    for(;x<=sx1;x++){ j=((j<<3)&0x1ff)|(r0[x+1]<<2)|(r1[x+1]<<1)|r2[x+1];
                        uint8_t v=T2[j]; o[x]=v; if(v!=r1[x]) h^=Z[y*W+x]; if(v){pop++; if(x<nx0)nx0=x; if(x>nx1)nx1=x; if(y<ny0)ny0=y; if(y>ny1)ny1=y;} } }
            }
            swap(a,b);
            if(pop==0){ R.cost=cost; return; }
            bx0=nx0;bx1=nx1;by0=ny0;by1=ny1; lastpop=pop; popH[t+1]=pop;
            // first contact with the grid border (cells there are never updated): check for linear growth now
            if(!R.growP && !edgeHit && (bx0<=2||by0<=2||bx1>=W-3||by1>=H-3)){ edgeHit=true; detectGrowth(R,t+1); }
            if(pop>C.pmax){ R.bad=1; R.cost=cost; detectGrowth(R,t+1); return; }
        }
        R.cost=cost;
        if(!stable){ R.bad=2; detectGrowth(R,t); }
        scan(R,true,bx0,bx1,by0,by1,h);
        if(!stable&&C.growth){ // never settled: a soup that keeps emitting copies of one ship looks like a gun
            const string*bk=nullptr; int bn=0; for(auto&kv:shipN) if(kv.second>bn){bn=kv.second;bk=&kv.first;}
            if(bk && bn>=C.gunShips){ R.gunN=min(bn,65535); gunShip.ex=keyCells(*bk); gunShip.I=cls1(gunShip.ex); gunShip.seed=seed; gunShip.count=bn; }
            // neither growing nor emitting nor at the border, yet still changing: a persistent bounded active region
            int e=t, h0=max(0,e-C.gmax/2); long mx=0,mx0=0; for(int g=h0;g<=e;g++) mx=max(mx,(long)popH[g]); for(int g=h0;g<=h0+50&&g<=e;g++) mx0=max(mx0,(long)popH[g]);
            if(!R.growP && !R.gunN && !edgeHit && e>=C.gmax && mx<=mx0*5/4+20) R.activeN=min<long>(popH[e],65535); }
    }
    // Probe runs: evolve an arbitrary pattern like a soup (escaping ships are deleted and logged in shipLog) until the
    // state repeats (period <= 60, twice) or gmax. The pattern is shifted by (offx,offy) onto the grid centre.
    // perGen(t) runs after every generation and may return false to stop. The grid is left clean afterwards, so
    // callers must read what they need inside perGen or via finalCells before the next run.
    vector<pair<int,int>> fin; // live cells (grid coords) at the end of the last runPattern
    bool runPattern(const Cells&init,int gmax,int&offx,int&offy,const function<bool(int)>&perGen){
        int mx,my,Mx,My; bbox(init,mx,my,Mx,My); offx=W/2-(mx+Mx)/2; offy=H/2-(my+My)/2;
        u64 h=0; SoupRes R; memset(&R,0,sizeof R); shipN.clear(); lastShipT=-1000000; curSeed=0; cands.clear();
        for(auto&p:init){ int i=(p.second+offy)*W+p.first+offx; a[i]=1; h^=Z[i]; }
        int bx0=mx+offx,bx1=Mx+offx,by0=my+offy,by1=My+offy, pbx0=bx0,pbx1=bx1,pby0=by0,pby1=by1;
        int ux0=bx0,ux1=bx1,uy0=by0,uy1=by1; u64 hist[256]; bool stable=false; long lastpop=init.size(); int t;
        for(t=0;t<gmax;t++){
            hist[t&255]=h;
            if(t>=20){ for(int p=1;p<=60&&2*p<=t;p++) if(hist[(t-p)&255]==h&&hist[(t-2*p)&255]==h){stable=true;break;} if(stable) break; }
            curT=t;
            if(t%8==0&&t>0&&lastpop<800) scan(R,false,bx0,bx1,by0,by1,h);
            int sx0=max(1,min(bx0,pbx0)-1), sx1=min(W-2,max(bx1,pbx1)+1), sy0=max(1,min(by0,pby0)-1), sy1=min(H-2,max(by1,pby1)+1);
            ux0=min(ux0,sx0); ux1=max(ux1,sx1); uy0=min(uy0,sy0); uy1=max(uy1,sy1);
            pbx0=bx0;pbx1=bx1;pby0=by0;pby1=by1; int nx0=INT_MAX,nx1=-1,ny0=INT_MAX,ny1=-1; long pop=0;
            for(int y=sy0;y<=sy1;y++){
                const uint8_t*r0=&a[(y-1)*W],*r1=&a[y*W],*r2=&a[(y+1)*W]; uint8_t*o=&b[y*W]; int x=sx0;
#ifdef HAVE_NEON
                for(;x+15<=sx1;x+=16){ uint8x16_t v=step16(r0,r1,r2,x,TBv); vst1q_u8(o+x,v);
                    u64 d=laneBits(veorq_u8(v,vld1q_u8(r1+x))); const u64*z=&Z[y*W+x]; while(d){ h^=z[__builtin_ctzll(d)>>2]; d&=d-1; }
                    u64 m=laneBits(v); if(m){ pop+=__builtin_popcountll(m); int f=x+(__builtin_ctzll(m)>>2), l=x+((63-__builtin_clzll(m))>>2);
                        if(f<nx0)nx0=f; if(l>nx1)nx1=l; if(y<ny0)ny0=y; if(y>ny1)ny1=y; } }
#endif
                int j= x<=sx1 ? (((r0[x-1]<<2)|(r1[x-1]<<1)|r2[x-1])<<3)|(r0[x]<<2)|(r1[x]<<1)|r2[x] : 0;
                for(;x<=sx1;x++){ j=((j<<3)&0x1ff)|(r0[x+1]<<2)|(r1[x+1]<<1)|r2[x+1];
                    uint8_t v=T2[j]; o[x]=v; if(v!=r1[x]) h^=Z[y*W+x]; if(v){pop++; if(x<nx0)nx0=x; if(x>nx1)nx1=x; if(y<ny0)ny0=y; if(y>ny1)ny1=y;} }
            }
            swap(a,b);
            if(pop==0){ bx0=1;bx1=0;by0=1;by1=0; t++; break; }
            bx0=nx0;bx1=nx1;by0=ny0;by1=ny1; lastpop=pop;
            if(pop>C.pmax||bx0<=2||by0<=2||bx1>=W-3||by1>=H-3){ t++; break; }
            if(perGen && !perGen(t+1)){ t++; break; }
        }
        fin.clear(); for(int y=by0;y<=by1;y++) for(int x=bx0;x<=bx1;x++) if(a[y*W+x]) fin.push_back({x,y});
        for(int y=max(0,uy0-1);y<=min(H-1,uy1+1);y++){ int x0=max(0,ux0-1),n=min(W-1,ux1+1)-x0+1; memset(&a[y*W+x0],0,n); memset(&b[y*W+x0],0,n); }
        return stable;
    }
    // population grows by exactly d>0 every P generations over the last 4 periods (puffer, rake, or blocked gun)
    void detectGrowth(SoupRes&R,int e){
        if(!C.growth||R.growP) return;
        for(int P=1;P<=300&&5*P<=e;P++){ int d=popH[e]-popH[e-P]; if(d<=0) continue; bool ok=true;
            for(int k=1;k<4&&ok;k++) ok = popH[e-k*P]-popH[e-(k+1)*P]==d;
            if(ok){ R.growP=P; R.growD=min(d,32767); return; } }
    }
    void scan(SoupRes&R,bool final,int bx0,int bx1,int by0,int by1,u64&h){
        if(st>1000000000){fill(stamp.begin(),stamp.end(),0);st=0;}
        st++; vector<int> live; for(int y=by0;y<=by1;y++) for(int x=bx0;x<=bx1;x++) if(a[y*W+x]) live.push_back(y*W+x);
        int sbase=st; vector<vector<int>> clusters; vector<int> q;
        for(int i:live){ if(stamp[i]>=sbase) continue; int id=clusters.size(); clusters.push_back({}); auto&cl=clusters.back();
            q.clear(); q.push_back(i); stamp[i]=sbase+id;
            for(size_t k=0;k<q.size();k++){ int p=q[k]; cl.push_back(p); int py=p/W,px=p%W;
                for(int dy=-2;dy<=2;dy++){int yy=py+dy; if(yy<0||yy>=H) continue; for(int dx=-2;dx<=2;dx++){int xx=px+dx; if(xx<0||xx>=W) continue; int n=yy*W+xx; if(a[n]&&stamp[n]<sbase){stamp[n]=sbase+id;q.push_back(n);} }} } }
        st=sbase+clusters.size();
        for(size_t id=0;id<clusters.size();id++){
            auto&cl=clusters[id]; if(cl.size()>400) continue;
            Cells c; for(int p:cl) c.push_back({p%W,p/W}); int mx,my,Mx,My; bbox(c,mx,my,Mx,My);
            if(Mx-mx>60||My-my>60) continue;
            // a ship counts only if nothing lies ahead in a cone along its direction of travel
            auto threatened=[&](int dx,int dy){
                int back=INT_MAX; for(auto&p:c) back=min(back,p.first*dx+p.second*dy);
                int cx=(mx+Mx)/2, cy=(my+My)/2;
                for(int p:live){ if(stamp[p]==(int)(sbase+id)) continue; int x=p%W,y=p/W; int along=x*dx+y*dy-back; if(along<-2) continue;
                    int lat=abs((x-cx)*dy-(y-cy)*dx)-(Mx-mx+My-my); if(lat<=along+6) return true; }
                return false; };
            // before the final scan only ships matter: skip classification when every direction is blocked
            if(!final && clusters.size()>1){ bool open=false; for(int d=0;d<9&&!open;d++) if(d!=4) open=!threatened(d%3-1,d/3-1); if(!open){ nskip++; continue; } }
            Info*Ip=&cls1(c); if(final&&Ip->type==4&&C.maxp>64) Ip=&cls2(c); Info&I=*Ip;
            string k;
            if(I.type==2){
                if(threatened((I.dx>0)-(I.dx<0),(I.dy>0)-(I.dy<0))) continue;
                k=keyOf(c); bool tgt=false; shipN[k]++; lastShipT=curT; lastShipK=k; if(shipLog) shipLog->push_back({c,I,curT});
                for(size_t t=0;t<C.tg.size();t++) if(C.tg[t].keys.count(k)){ if(R.cnt[t]<255) R.cnt[t]++; tgt=true; }
                if(!tgt) note(c,I);
                for(int p:cl){a[p]=0; h^=Z[p];}
            } else if(final && I.type==1){
                bool ok=true;
                for(auto&e:I.env){int ex=e.first+mx,ey=e.second+my;
                    for(int yy=ey-2;yy<=ey+2&&ok;yy++) for(int xx=ex-2;xx<=ex+2;xx++){ if(xx<0||yy<0||xx>=W||yy>=H) continue; int n=yy*W+xx; if(a[n]&&stamp[n]!=(int)(sbase+id)){ok=false;break;} } if(!ok)break; }
                if(!ok) continue;
                k=keyOf(c); bool tgt=false;
                if(collectAsh){ auto&e=ash[k]; if(!e.count){e.ex=c;e.I=I;e.seed=curSeed;} e.count++; }
                for(size_t t=0;t<C.tg.size();t++) if(C.tg[t].keys.count(k)){ if(R.cnt[t]<255) R.cnt[t]++; tgt=true; }
                if(!tgt && I.period>=C.minOscPeriod) note(c,I);
            }
        }
    }
};
static vector<Engine*> engines;

// ===================== evaluation =====================
struct Eval {
    Table T; int n=0; vector<SoupOut> res; vector<double> usage; // per class
    unordered_map<string,Cand> cands, ash; double f=-1e9; vector<double> rate; double bad=0, cost=0; u64 seedEpoch=0;
};
static double tSim=0,tProc=0; static double nowS(){ return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }
static u64 soupSeed(u64 epoch,int i){ u64 s=epoch*1000003ULL+i; return splitmix(s); }
static void evalRange(Eval&E,int i0,int i1){
    double t0=nowS(); if((int)E.res.size()<i1) E.res.resize(i1);
    #pragma omp parallel
    { Engine&G=*engines[omp_get_thread_num()];
      #pragma omp for schedule(dynamic,1)
      for(int i=i0;i<i1;i++) G.run(soupSeed(E.seedEpoch,i),E.res[i]); }
    E.n=i1; tSim+=nowS()-t0;
}
static void beginEval(Eval&E,const Table&T,u64 epoch){
    E.T=T; E.n=0; E.res.clear(); E.cands.clear(); E.ash.clear(); E.seedEpoch=epoch;
    for(auto*g:engines) g->setRule(E.T.t);
}
static void finishEval(Eval&E){
    // merge candidates in soup order
    E.usage.clear();
    for(int i=0;i<E.n;i++){ auto&R=E.res[i];
        for(auto&kv:R.cands){ auto&e=E.cands[kv.first]; if(!e.count){e.ex=kv.second.ex;e.I=kv.second.I;e.seed=kv.second.seed;} e.count+=kv.second.count; }
        for(auto&kv:R.ash){ auto&e=E.ash[kv.first]; if(!e.count){e.ex=kv.second.ex;e.I=kv.second.I;e.seed=kv.second.seed;} e.count+=kv.second.count; } }
}
// transition usage per class over the first 64 soups (mutation weights). Computed on demand, only for rules that
// are actually mutated from, by re-running those soups with counting on (the counting sweep is scalar).
static const vector<double>& usageOf(Eval&E){
    if(!E.usage.empty()) return E.usage;
    double t0=nowS(); int n=min(64,E.n); vector<vector<uint32_t>> per(n);
    for(auto*g:engines) g->setRule(E.T.t);
    #pragma omp parallel
    { Engine&G=*engines[omp_get_thread_num()]; G.countUsage=true;
      #pragma omp for schedule(dynamic,1)
      for(int i=0;i<n;i++){ SoupOut o; G.run(soupSeed(E.seedEpoch,i),o); per[i]=std::move(o.usage); }
      G.countUsage=false; }
    vector<u64> u(512,0); for(auto&p:per) for(int j=0;j<512;j++) u[j]+=p[j];
    E.usage.assign(NC,0); for(int j=0;j<512;j++) E.usage[cls[P2[j]]]+=u[j];
    tSim+=nowS()-t0; return E.usage;
}
static void score(Eval&E){
    int n=E.n, nt=C.tg.size(); E.rate.assign(nt,0); double bad=0,cost=0;
    for(int i=0;i<n;i++){ auto&r=E.res[i].r; cost+=r.cost; if(r.bad){bad++;continue;} for(int t=0;t<nt;t++) E.rate[t]+=min((int)r.cnt[t],C.cap); }
    for(auto&x:E.rate) x/=n; E.bad=bad/n; E.cost=cost/n/1e6;
    double f=0; for(int t=0;t<nt;t++) f+=C.tg[t].w*log(E.rate[t]+1e-3);
    if(C.costW>0 && C.costRef>0) f-=C.costW*max(0.0,log(E.cost/C.costRef));
    // feasibility penalties
    for(int t=0;t<nt;t++) if(C.tg[t].minRate>0) f-=5*max(0.0,(C.tg[t].minRate-E.rate[t])/C.tg[t].minRate);
    f-=5*max(0.0,E.bad-C.maxBad);
    if(C.maxCost>0) f-=5*max(0.0,log(E.cost/C.maxCost));
    E.f=f;
}
static bool feasible(const Eval&E){
    for(size_t t=0;t<C.tg.size();t++) if(E.rate[t]<C.tg[t].minRate) return false;
    if(E.bad>C.maxBad) return false; if(C.maxCost>0&&E.cost>C.maxCost) return false; return true;
}
static void fullEval(Eval&E,const Table&T,u64 epoch){ beginEval(E,T,epoch); evalRange(E,0,C.N); finishEval(E); score(E); }

// ===================== discoveries =====================
static set<string> archiveIds; static map<string,int> descCount; static FILE* logFp=nullptr;
static string descOf(const string&id){ return id.substr(0,id.find(':')); }
static string canonId(const Cells&c,const Info&I,const uint8_t*T){
    string best; int P=max(1,I.period);
    for(int o:groupOrients()){ Cells x=xform(c,o); for(int g=0;g<P;g++){ string k=keyOf(x); if(best.empty()||k<best) best=k; x=stepCells(x,T);} }
    int ax=abs(I.dx),ay=abs(I.dy); if(SYM!=2&&ax<ay) swap(ax,ay); int dx= SYM==2?I.dx:ax, dy= SYM==2?I.dy:ay;
    char buf[64]; if(I.type==2) snprintf(buf,64,"Q%d(%d,%d)",I.period,dx,dy); else snprintf(buf,64,"P%d",I.period);
    string hex; for(unsigned char ch:best){ char b2[3]; snprintf(b2,3,"%02x",ch); hex+=b2; }
    return string(buf)+":"+hex;
}
static string describe(const Info&I){
    char b[96]; if(I.type==1){ snprintf(b,96,"p%d oscillator",I.period); return b; }
    int ax=abs(I.dx),ay=abs(I.dy); if(ax<ay) swap(ax,ay); int g=std::gcd(ax,std::gcd(ay,I.period));
    if(ay==0) snprintf(b,96,"%s/%d orthogonal ship (p%d)", ax/g==1?"c":(to_string(ax/g)+"c").c_str(), I.period/g, I.period);
    else if(ax==ay) snprintf(b,96,"%s/%d diagonal ship (p%d)", ax/g==1?"c":(to_string(ax/g)+"c").c_str(), I.period/g, I.period);
    else snprintf(b,96,"(%d,%d)c/%d oblique ship",ax,ay,I.period);
    return b;
}
struct Find { string id; long count; Cells phase; Info I; u64 seed; };
// returns list of genuine (non-pseudo) interesting objects in this evaluation, aggregated by canonical id
static vector<Find> processCands(const Eval&E){
    double t0=nowS(); struct TT{double t0; ~TT(){ tProc+=nowS()-t0; }} tt{t0};
    map<string,Find> by;
    vector<const pair<const string,Cand>*> ord; for(auto&kv:E.cands) ord.push_back(&kv); sort(ord.begin(),ord.end(),[](auto*a,auto*b){return a->first<b->first;}); // deterministic representative
    for(auto*kvp:ord){
        const Cand&c=kvp->second; if(isPseudo(c.ex,E.T.t,c.I.period)) continue;
        string id=canonId(c.ex,c.I,E.T.t);
        auto it=by.find(id);
        if(it==by.end()){ // minimal-population phase
            Cells x=c.ex,m=c.ex; for(int g=0;g<c.I.period;g++){ x=stepCells(x,E.T.t); if(x.size()<m.size()) m=x; }
            by[id]={id,c.count,m,c.I,c.seed};
        } else it->second.count+=c.count;
    }
    vector<Find> v; for(auto&kv:by) v.push_back(kv.second);
    sort(v.begin(),v.end(),[](const Find&a,const Find&b){return a.count>b.count;});
    return v;
}
// novelty credit: new behaviour (type, period, velocity) scores high; a new shape of a known behaviour scores little
static double rewardMult(const Info&I){
    if(C.reward.empty()) return 1; double m=1; auto f=[&](const char*k){ auto it=C.reward.find(k); if(it!=C.reward.end()) m*=it->second; };
    if(I.type==1) f("osc"); else { f("ship"); int ax=abs(I.dx),ay=abs(I.dy); f(ax==0||ay==0?"orthogonal": ax==ay?"diagonal":"oblique"); }
    return m; }
static double interest(const Find&F,bool newDesc){ return rewardMult(F.I)*(!newDesc ? 0.1 : F.I.type==2 ? 5.0 : 1.0+log2((double)F.I.period)/2); }
static double recordFinds(const Eval&E,const vector<Find>&fs,const string&ctx){
    double nnew=0; string rule=ruleString(E.T);
    for(auto&F:fs){
        if(archiveIds.count(F.id)) continue;
        bool nd = descCount[descOf(F.id)]++==0;
        archiveIds.insert(F.id); nnew+=interest(F,nd);
        if(!C.archiveFile.empty()){ FILE*a=fopen(C.archiveFile.c_str(),"a"); if(a){ fprintf(a,"%s\n",F.id.c_str()); fclose(a);} }
        string msg=string(nd?"# NEW BEHAVIOUR ":"# new shape ")+describe(F.I)+", "+to_string(F.phase.size())+" cells, "+to_string(F.count)+" in "+to_string(E.n)+" soups ("+ctx+")\n# first seen in soup seed "+to_string(F.seed)+"\n"+toRLE(F.phase,rule);
        if(nd||F.I.type==2){ printf("%s",msg.c_str()); fflush(stdout); }
        if(logFp){ fprintf(logFp,"%s\n",msg.c_str()); fflush(logFp); }
    }
    return nnew;
}
// soups that never settled: linear growth (puffer/rake/gun) and ship streams (gun). Logged with the soup itself.
struct SoupEvent { string id, desc; long count=0; u64 seed=0; double credit=0; };
static vector<SoupEvent> soupEvents(const Eval&E){
    map<string,SoupEvent> by;
    for(int i=0;i<E.n;i++){ auto&R=E.res[i]; u64 seed=soupSeed(E.seedEpoch,i);
        auto add=[&](const string&id,const string&desc,double credit){ auto&e=by[id]; if(!e.count){e.id=id;e.desc=desc;e.seed=seed;e.credit=credit;} e.count++; };
        if(R.r.growP){ char b[160]; snprintf(b,160,"linear growth +%d cells per %d gens (puffer, rake or gun)",R.r.growD,R.r.growP);
            add("G"+to_string(R.r.growP)+":"+to_string(R.r.growD),b,3.0*(C.reward.count("growth")?C.reward.at("growth"):1)); }
        if(R.r.activeN){ int b=0; while((2<<b)<=R.r.activeN) b++; char t[160]; snprintf(t,160,"persistent bounded active region (~%d cells; never settles: aperiodic, or a period > 60)",1<<b);
            add("A"+to_string(1<<b),t,2.0*(C.reward.count("active")?C.reward.at("active"):1)); }
        if(R.r.gunN && R.gunShip.count){ string cid=canonId(R.gunShip.ex,R.gunShip.I,E.T.t);
            add("U"+descOf(cid), describe(R.gunShip.I)+" emitter (gun?): the soup keeps producing it ("+to_string(R.r.gunN)+" copies seen)",
                10.0*(C.reward.count("gun")?C.reward.at("gun"):1)); } }
    vector<SoupEvent> v; for(auto&kv:by) v.push_back(kv.second); return v;
}
// re-run a soup without ship deletion: a gun keeps adding ships, so the population grows between gen 300 and 600
static int emitterGrowth(const Cells&soup,const uint8_t*T){ Cells c=soup; size_t p300=0; for(int g=1;g<=600&&!c.empty();g++){ c=stepCells(c,T); if(g==300) p300=c.size(); if(c.size()>20000) break; } return (int)c.size()-(int)p300; }
static double recordSoupEvents(const Eval&E,const string&ctx){
    double credit=0; string rule=ruleString(E.T);
    for(auto e:soupEvents(E)){
        if(archiveIds.count(e.id)) continue; bool nd = descCount[descOf(e.id)]++==0; archiveIds.insert(e.id); credit+= nd? e.credit : 0.1*e.credit;
        if(!C.archiveFile.empty()){ FILE*a=fopen(C.archiveFile.c_str(),"a"); if(a){ fprintf(a,"%s\n",e.id.c_str()); fclose(a);} }
        int fw,fh; Cells soup=makeSoup(e.seed,fw,fh);
        if(e.id[0]=='U'){ int g=emitterGrowth(soup,E.T.t);
            e.desc+= g>=20 ? "; without ship deletion the population grows by "+to_string(g)+" from gen 300 to 600 (consistent with a gun)"
                         : "; probably NOT a gun: without ship deletion the population does not grow (likely a larger moving object whose front is this ship)"; }
        string msg=string(nd?"# NEW SOUP BEHAVIOUR ":"# soup behaviour ")+e.desc+", "+to_string(e.count)+" in "+to_string(E.n)+" soups ("+ctx+")\n# soup seed "+to_string(e.seed)+" (the soup is below; run it in Golly)\n"+toRLE(soup,rule);
        printf("%s",msg.c_str()); fflush(stdout); if(logFp){ fprintf(logFp,"%s\n",msg.c_str()); fflush(logFp); }
    }
    return credit;
}
// ===================== probe: catalysts, signals, fossils =====================
static void printEval(const char*tag,const Eval&E);
// Engines for probe runs use a larger grid (C.probeGrid) so products and long wires fit.
static vector<Engine*> probeEngines;
static void probeSetRule(const Table&T){ int nt=omp_get_max_threads(); while((int)probeEngines.size()<nt) probeEngines.push_back(new Engine(C.probeGrid)); for(auto*g:probeEngines) g->setRule(T.t); }
static Cells normCells(const Cells&c){ int mx,my,Mx,My; bbox(c,mx,my,Mx,My); Cells r; for(auto&p:c) r.push_back({p.first-mx,p.second-my}); sort(r.begin(),r.end()); return r; }
struct LibObj { Cells c; Info I; long count=0; string canon; };
struct ProbeHit { string id, text; long lanes=0; Cells setup; double credit=0, score=0; long sig=0; };
static string noSpace(string s){ for(auto&ch:s) if(ch==' ') ch='_'; return s; }

// objects for the probe: common ships (targets + found ships) and common still lifes from the final ash
static void probeLibrary(const Eval&E,vector<LibObj>&ships,vector<LibObj>&stills){
    map<string,LibObj> sh,sl;
    for(size_t t=0;t<C.tg.size();t++){ auto&P=C.pats[C.tg[t].pat]; if(!(C.tg[t].dx||C.tg[t].dy)) continue;
        LibObj o; o.c=normCells(P.evo.gens[P.evo.s]); o.I=classify(o.c,E.T.t,256); o.count=lround(E.rate[t]*E.n); if(o.I.type!=2) continue;
        o.canon=canonId(o.c,o.I,E.T.t); auto&e=sh[o.canon]; if(!e.count) e=o; else e.count+=o.count; }
    for(auto&kv:E.cands){ auto&c=kv.second; if(c.I.type!=2||isPseudo(c.ex,E.T.t,c.I.period)) continue; LibObj o; o.c=normCells(c.ex); o.I=c.I; o.count=c.count; o.canon=canonId(o.c,o.I,E.T.t);
        auto&e=sh[o.canon]; if(!e.count) e=o; else e.count+=o.count; }
    for(auto&kv:E.ash){ auto&c=kv.second; if(c.I.type!=1||c.I.period!=1||isPseudo(c.ex,E.T.t,1)) continue; LibObj o; o.c=normCells(c.ex); o.I=c.I; o.count=c.count; o.canon=canonId(o.c,o.I,E.T.t);
        auto&e=sl[o.canon]; if(!e.count) e=o; else e.count+=o.count; }
    auto pick=[](map<string,LibObj>&m,int k,long minc){ vector<LibObj> v; for(auto&kv:m) if(kv.second.count>=minc) v.push_back(kv.second);
        sort(v.begin(),v.end(),[](const LibObj&a,const LibObj&b){ return a.count!=b.count? a.count>b.count : a.canon<b.canon; }); if((int)v.size()>k) v.resize(k); return v; };
    ships=pick(sh,C.probeShips,2); stills=pick(sl,C.probeStills,2);
}
static string objName(const Info&I,size_t n){ return (I.type==1&&I.period==1? string("still life"): describe(I))+" ("+to_string(n)+" cells)"; }

// Fire every library ship at every library still life, in every orientation of the still life and every lane
// (one ship phase suffices: a still life is time-invariant). Keep reactions that restore the still life exactly.
static vector<ProbeHit> probeCatalysts(const Eval&E,const vector<LibObj>&ships,const vector<LibObj>&stills){
    struct Setup { int si,li; Cells sl, ship; int sx,sy; };
    vector<Setup> setups;
    for(size_t si=0;si<ships.size();si++){
        const LibObj&S=ships[si]; int vx=S.I.dx,vy=S.I.dy, P=S.I.period;
        int smx=INT_MAX,smy=INT_MAX,sMx=INT_MIN,sMy=INT_MIN; { Cells x=S.c; int ox=0,oy=0; // union box over one period, relative to phase 0
            for(int g=0;g<P;g++){ int a,b,A,B; bbox(x,a,b,A,B); smx=min(smx,a);smy=min(smy,b);sMx=max(sMx,A);sMy=max(sMy,B); x=stepCells(x,E.T.t); } (void)ox;(void)oy; }
        for(size_t li=0;li<stills.size();li++){ set<string> seenO;
            for(int o:groupOrients()){ Cells L=normCells(xform(stills[li].c,o)); if(!seenO.insert(keyOf(L)).second) continue;
                int lw,lh; { int a,b,A,B; bbox(L,a,b,A,B); lw=A; lh=B; }
                auto dist=[&](int px,int py){ int ax0=px+smx,ax1=px+sMx,ay0=py+smy,ay1=py+sMy; int dx=max(0,max(ax0-lw,-ax1)), dy=max(0,max(ay0-lh,-ay1)); return max(dx,dy); };
                set<pair<int,int>> starts; int R=3+max(sMx-smx,sMy-smy), K=4+(max(lw,lh)+2*R)/max(1,max(abs(vx),abs(vy)));
                for(int py=-R-sMy;py<=lh+R-smy;py++) for(int px=-R-sMx;px<=lw+R-smx;px++){
                    int kmin=INT_MAX; for(int k=-K;k<=K;k++) if(dist(px+k*vx,py+k*vy)<=2){ kmin=k; break; }
                    if(kmin==INT_MAX) continue;
                    int k0=kmin; while(dist(px+k0*vx,py+k0*vy)<5) k0--; // first position along the track at distance >= 5
                    starts.insert({px+k0*vx,py+k0*vy}); }
                for(auto&st:starts) setups.push_back({(int)si,(int)li,L,S.c,st.first,st.second});
            } }
    }
    vector<vector<ProbeHit>> res(setups.size());
    probeSetRule(E.T);
    #pragma omp parallel
    { Engine&G=*probeEngines[omp_get_thread_num()]; vector<ShipEvt> log; G.shipLog=&log;
      #pragma omp for schedule(dynamic,4)
      for(size_t q=0;q<setups.size();q++){ auto&U=setups[q]; log.clear();
          Cells init=U.sl; for(auto&p:U.ship) init.push_back({p.first+U.sx,p.second+U.sy});
          int ox,oy; int lw=0,lh=0; for(auto&p:U.sl){ lw=max(lw,p.first); lh=max(lh,p.second); }
          // interaction = the region around the still life differs from still life + ship evolving independently
          bool touched=false; Cells sc; for(auto&p:U.ship) sc.push_back({p.first+U.sx,p.second+U.sy});
          bool settled=G.runPattern(init,C.probeGens,ox,oy,[&](int){ if(touched) return true;
              sc=stepCells(sc,E.T.t); set<pair<int,int>> ideal(U.sl.begin(),U.sl.end()); ideal.insert(sc.begin(),sc.end());
              for(int y=-2;y<=lh+2&&!touched;y++) for(int x=-2;x<=lw+2;x++) if(G.a[(y+oy)*G.W+x+ox]!=(uint8_t)ideal.count({x,y})){ touched=true; break; }
              return true; });
          if(!settled||!touched) continue;
          // final state: the still life must be back exactly, everything else must be identifiable
          Cells fc; for(auto&p:G.fin) fc.push_back({p.first-ox,p.second-oy});
          bool restored=false; vector<string> prods; vector<Info> shipsOut; bool unknown=false;
          set<pair<int,int>> SL(U.sl.begin(),U.sl.end());
          for(auto&cl:splitComponents(fc)){
              if(set<pair<int,int>>(cl.begin(),cl.end())==SL){ restored=true; continue; }
              Info I=G.cls1(cl); if(I.type!=1){ unknown=true; break; } prods.push_back(objName(I,cl.size())); }
          if(!restored||unknown) continue;
          for(auto&ev:log){ shipsOut.push_back(ev.I); prods.push_back(objName(ev.I,ev.c.size())); }
          sort(prods.begin(),prods.end());
          const LibObj&S=ships[U.si]; string cat; double credit=1;
          bool onlyShips=shipsOut.size()==prods.size();
          if(prods.empty()){ cat="eater"; credit=1; }
          else if(onlyShips && shipsOut.size()==1){
              Info o=shipsOut[0]; Cells oc=log[0].c; string oid=canonId(normCells(oc),o,E.T.t);
              if(oid==S.canon){ if(o.dx==S.I.dx&&o.dy==S.I.dy){ cat="pass-through"; credit=4; }
                  else { double a=atan2((double)o.dy,(double)o.dx)-atan2((double)S.I.dy,(double)S.I.dx); int deg=(int)lround(fabs(remainder(a,2*M_PI))*180/M_PI);
                      cat="reflector "+to_string(deg)+"deg"; credit=8; } }
              else { cat="converter to "+objName(o,oc.size()); credit=5; } }
          else if(onlyShips){ cat="splitter to "+to_string(shipsOut.size())+" ships"; credit=5; }
          else { cat="catalysed reaction"; credit=1.5; }
          string ptxt; for(auto&p:prods) ptxt+=(ptxt.empty()?"":" + ")+p;
          ProbeHit hit; hit.id="K"+noSpace(cat)+"_"+descOf(S.canon)+":"+stills[U.li].canon.substr(stills[U.li].canon.find(':')+1)+":"+noSpace(ptxt);
          hit.text=cat+": "+objName(S.I,S.c.size())+" + "+objName(stills[U.li].I,U.sl.size())+" -> still life restored"+(ptxt.empty()?"":" + "+ptxt);
          hit.lanes=1; hit.setup=init; hit.credit=credit*(C.reward.count("catalyst")?C.reward.at("catalyst"):1); res[q].push_back(hit); }
      G.shipLog=nullptr; }
    map<string,ProbeHit> by; for(auto&v:res) for(auto&h:v){ auto it=by.find(h.id); if(it==by.end()) by[h.id]=h; else it->second.lanes++; }
    vector<ProbeHit> out; for(auto&kv:by) out.push_back(kv.second);
    stable_sort(out.begin(),out.end(),[](const ProbeHit&a,const ProbeHit&b){ return a.credit>b.credit; });
    return out;
}

// Wicks: stable periodic lines built from a small tile repeated along (k,0) or (k,k). Signals: a one-cell disturbance
// in the middle of a long wick whose damage (cells differing from the intact wick) stays localized at its front(s)
// and moves at a constant speed along it. 1D CA: the damage spreads along the wick at a constant speed without ever
// leaving the wick's own rows, so the wick acts as a one-dimensional automaton (holes multiplying inside it).
// Fuses: the damage front moves at constant speed while the damage grows and spills out of the wick.
static vector<ProbeHit> probeWicks(const Eval&E,int*nWicks){
    const uint8_t*T=E.T.t; struct Wick { Cells tile; int sx,sy,bw,bh; string sig; };
    vector<Wick> wicks; set<string> sigs;
    auto line=[](const Wick&w,int n){ Cells c; for(int i=0;i<n;i++) for(auto&p:w.tile) c.push_back({p.first+i*w.sx,p.second+i*w.sy}); return c; };
    auto proj=[](const Wick&w,int x,int y){ return w.sy? x+y : x; };
    for(int diag=0;diag<2;diag++) for(int k=1;k<=C.wickK;k++) for(int h=1;h<=(diag?k:3);h++){
        int bw=k,bh=diag?k:h; int nb=bw*bh; if(nb>16) continue;
        for(u64 m=1;m<(1ULL<<nb);m++){ if(__builtin_popcountll(m)>C.wickCells) continue;
            Wick w; w.sx=k; w.sy=diag?k:0; w.bw=bw; w.bh=bh;
            for(int b=0;b<nb;b++) if(m>>b&1) w.tile.push_back({b%bw,b/bw});
            if(!diag){ bool top=false,bot=false; for(auto&p:w.tile){ top|=p.second==0; bot|=p.second==h-1; } if(!top||!bot) continue; }
            else { bool l=false,t=false; for(auto&p:w.tile){ l|=p.first==0; t|=p.second==0; } if(!l||!t) continue; }
            int Q=diag?2*k:k; Cells L=line(w,16), L2=stepCells(L,T);
            auto win=[&](const Cells&c,int lo,int hi){ set<pair<int,int>> s; for(auto&p:c){ int q=proj(w,p.first,p.second); if(q>=lo&&q<hi) s.insert(p);} return s; };
            if(win(L,5*Q,11*Q)!=win(L2,5*Q,11*Q)) continue;                 // interior must be still
            bool trivial=false; for(auto&cc:comps8(L)){ int lo=INT_MAX,hi=INT_MIN; for(auto&p:cc){ int q=proj(w,p.first,p.second); lo=min(lo,q); hi=max(hi,q);} 
                if(lo>=3*Q&&hi<13*Q){ Cells s=stepCells(cc,T); if(set<pair<int,int>>(s.begin(),s.end())==set<pair<int,int>>(cc.begin(),cc.end())){ trivial=true; break; } } }
            if(trivial) continue;                                           // just a row of independent still lifes
            string best; Cells L3=line(w,24);
            for(int o:groupOrients()) for(int j=0;j<Q;j++){ Cells x; for(auto&p:L3){ int q=proj(w,p.first,p.second); if(q>=6*Q+j&&q<6*Q+j+24) x.push_back(p);} string kk=keyOf(xform(x,o)); if(best.empty()||kk<best) best=kk; }
            if(!sigs.insert(best).second) continue;
            w.sig=best; wicks.push_back(w);
        } }
    if(nWicks) *nWicks=wicks.size();
    struct Run { int wi; pair<int,int> cell; };
    vector<Run> runs; for(size_t wi=0;wi<wicks.size()&&wi<(size_t)C.wickMax;wi++){ auto&w=wicks[wi]; int n=(C.probeGrid*3/5)/(w.sy?w.sx:w.sx); int mid=n/2;
        for(int y=-1;y<=w.bh;y++) for(int x=-1;x<=w.bw;x++) runs.push_back({(int)wi,{x+mid*w.sx,y+mid*w.sy}}); }
    vector<vector<ProbeHit>> res(runs.size()); probeSetRule(E.T);
    const int G=64; // speed-c fronts would reach the light cone of the wick's ends soon after
    #pragma omp parallel
    { Engine&Gn=*probeEngines[omp_get_thread_num()];
      #pragma omp for schedule(dynamic,2)
      for(size_t q=0;q<runs.size();q++){ auto&w=wicks[runs[q].wi]; int n=(C.probeGrid*3/5)/w.sx;
          Cells L=line(w,n); set<pair<int,int>> Ls(L.begin(),L.end());
          Cells init=L; auto pc=runs[q].cell; if(Ls.count(pc)) init.erase(find(init.begin(),init.end(),pc)); else init.push_back(pc);
          int qlo=INT_MAX,qhi=INT_MIN; for(auto&p:L){ int v=proj(w,p.first,p.second); qlo=min(qlo,v); qhi=max(qhi,v); }
          int q0=proj(w,pc.first,pc.second), sp=w.sy?2:1; int ox,oy; vector<int> lo(G+1,INT_MAX),hi(G+1,INT_MIN),cnt(G+1,0),wire(G+1,0);
          // last TL gens: damage positions along the wick, and damage outside the wick's band of rows
          const int TL=8; vector<vector<pair<int,int>>> dv(TL+1); vector<int> outBand(TL+1,0);
          auto perp=[&](int x,int y){ return w.sy? x-y : y; }; int pmin=INT_MAX,pmax=INT_MIN; for(auto&p:L){ int u=perp(p.first,p.second); pmin=min(pmin,u); pmax=max(pmax,u); }
          // ideal: the intact wick (stable away from its ends); compared inside the light cone of the ends
          int mx,my,Mx,My; bbox(init,mx,my,Mx,My); int IW=Mx-mx+40, IH=My-my+40; vector<uint8_t> ideal(IW*IH,0); for(auto&p:L) ideal[(p.second-my+20)*IW+p.first-mx+20]=1;
          Gn.runPattern(init,G,ox,oy,[&](int t){
              int a=qlo+sp*t+4, b=qhi-sp*t-4; if(a>=b) return false;
              for(int y=my-20;y<My+20;y++) for(int x=mx-20;x<Mx+20;x++){ int v=proj(w,x,y); if(v<a||v>b) continue;
                  int gx=x+ox, gy=y+oy; if(gx<0||gy<0||gx>=Gn.W||gy>=Gn.H) continue; uint8_t is=Gn.a[gy*Gn.W+gx], id=ideal[(y-my+20)*IW+x-mx+20];
                  if(is!=id){ cnt[t]++; if(id) wire[t]++; lo[t]=min(lo[t],v); hi[t]=max(hi[t],v);
                      if(t>=G-TL){ dv[t-(G-TL)].push_back({x,y}); int u=perp(x,y); if(u<pmin||u>pmax) outBand[t-(G-TL)]++; } } }
              return true; });
          auto ok=[&](int t){ return t<=G && cnt[t]>0; };
          if(!ok(G)||!ok(G/2)||!ok(3*G/4)) continue;
          int R1=hi[3*G/4]-hi[G/2], R2=hi[G]-hi[3*G/4], L1=lo[G/2]-lo[3*G/4], L2=lo[3*G/4]-lo[G];
          bool right=R1>0&&abs(R1-R2)<=1, left=L1>0&&abs(L1-L2)<=1; if(!right&&!left) continue;
          int c1=cnt[G/2],c2=cnt[3*G/4],c3=cnt[G]; bool bounded=c3<=c1+max(4,c1/5)&&c2<=c1+max(4,c1/5);
          bool growing=c2-c1>0&&c3-c2>0&&abs((c3-c2)-(c2-c1))<=max(2,(c2-c1)/5);
          if(wire[G]==0) continue;                                          // must actually alter the wick
          // localized: in each of the last TL gens all damage lies within Wd of the leading/trailing edge (a 1D CA
          // fills the space between them); inBand: no damage outside the wick's rows in those gens
          int Wd=12*sp+2*(w.sy?2*w.sx:w.sx); bool localized=true, inBand=true;
          for(int i=0;i<=TL;i++){ int t=G-TL+i; if(outBand[i]) inBand=false; if(cnt[t]==0) continue;
              for(auto&c:dv[i]){ int v=proj(w,c.first,c.second); if(v>lo[t]+Wd&&v<hi[t]-Wd){ localized=false; break; } } }
          // a 1D CA keeps changing behind its front; a fuse leaves static or short-period remains there
          bool periodicBehind=false; { int a=lo[G-TL]+Wd, b=hi[G-TL]-Wd;
              auto behind=[&](int i){ set<pair<int,int>> r; for(auto&c:dv[i]){ int v=proj(w,c.first,c.second); if(v>a&&v<b) r.insert(c);} return r; };
              auto now=behind(TL); for(int p=1;p<=TL&&!periodicBehind;p++) periodicBehind= behind(TL-p)==now; }
          string kind= bounded&&localized? "signal" : inBand&&!localized&&!periodicBehind? "1dca" : growing? "fuse" : ""; if(kind.empty()) continue;
          // exact speed: smallest period P with constant displacement of the leading edge over three periods
          int dq=0,P=0; for(int p=1;p<=40&&3*p<=G/2;p++){ int e=right? hi[G]:-lo[G], e1=right? hi[G-p]:-lo[G-p], e2=right? hi[G-2*p]:-lo[G-2*p], e3=right? hi[G-3*p]:-lo[G-3*p];
              if(e-e1>0&&e-e1==e1-e2&&e1-e2==e2-e3&&(kind!="signal"||cnt[G]==cnt[G-p])){ dq=e-e1; P=p; break; } }
          if(!P) continue;                                                   // no exact periodic front: just drifting damage
          if(kind=="signal"){ // confirm: follow the front further (until the ends' light cone closes in) at the same speed;
              // slow damage often creeps for a while and then stops
              const int G2=2*G; vector<int> lo2(G2+1,INT_MAX),hi2(G2+1,INT_MIN); int tE=0;
              Gn.runPattern(init,G2,ox,oy,[&](int t){
                  int a=qlo+sp*t+4, b=qhi-sp*t-4; if(a>=b) return false; if(t<=G) return true;
                  for(int y=my-20;y<My+20;y++) for(int x=mx-20;x<Mx+20;x++){ int v=proj(w,x,y); if(v<a||v>b) continue;
                      int gx=x+ox, gy=y+oy; if(gx<0||gy<0||gx>=Gn.W||gy>=Gn.H) continue;
                      if(Gn.a[gy*Gn.W+gx]!=ideal[(y-my+20)*IW+x-mx+20]){ lo2[t]=min(lo2[t],v); hi2[t]=max(hi2[t],v); } }
                  if(lo2[t]<=hi2[t]){ int e=right? hi2[t]:-lo2[t]; if(e+Wd>=(right? b:-a)) return false; tE=t; } // front near the window edge
                  return true; });
              if(tE>=G+P){ int e0=right? hi[G]:-lo[G], e1=right? hi2[tE]:-lo2[tE]; if(abs((e1-e0)*P-dq*(tE-G))>2*P) continue; } }
          // speed: on a diagonal wick proj = x+y moves 2 per cell step
          char sbuf[64]; if(P){ int d=dq, PP=w.sy? 2*P:P; int g=std::gcd(d,PP); snprintf(sbuf,64,"%s/%d",d/g==1?"c":(to_string(d/g)+"c").c_str(),PP/g); } else snprintf(sbuf,64,"~%.3fc",(double)(right?R1+R2:L1+L2)/(G/2)/(w.sy?2:1));
          string dir= left&&right? "both ways" : "one way";
          ProbeHit hit; hit.id="S"+kind+"_"+string(sbuf)+(w.sy?"_diag":"_orth")+":"; for(unsigned char ch:w.sig){ char b2[3]; snprintf(b2,3,"%02x",ch); hit.id+=b2; }
          hit.text=(kind=="1dca"?string("1D CA"):kind)+" at "+sbuf+" along a "+(w.sy?"diagonal":"orthogonal")+" wick ("+dir+")"+(kind=="signal"&&cnt[G]>0&&lo[G]<=q0&&hi[G]>=q0&&!(left&&right)?"":"");
          Cells seg; int m=n/2; for(auto&p:init){ int u=w.sy? (p.first)/w.sx : p.first/w.sx; if(abs(u-m)<=20) seg.push_back(p); } hit.setup=seg;
          hit.credit=(kind=="signal"?6:kind=="1dca"?5:2)*(C.reward.count("signal")?C.reward.at("signal"):1); hit.lanes=1; res[q].push_back(hit); } }
    map<string,ProbeHit> by; for(auto&v:res) for(auto&h:v){ auto it=by.find(h.id); if(it==by.end()) by[h.id]=h; else it->second.lanes++; }
    vector<ProbeHit> out; for(auto&kv:by) out.push_back(kv.second); return out;
}

// fossils: large isolated objects in the final ash that occur in few soups (markers of rare reactions)
static vector<ProbeHit> probeFossils(const Eval&E){
    map<string,ProbeHit> by;
    for(auto&kv:E.ash){ auto&c=kv.second; if((int)c.ex.size()<C.fossilCells) continue; if(isPseudo(c.ex,E.T.t,max(1,c.I.period))) continue; string cid=canonId(normCells(c.ex),c.I,E.T.t);
        auto it=by.find(cid); if(it!=by.end()){ it->second.lanes+=c.count; continue; }
        ProbeHit h; h.id=string(c.I.period==1?"Fstill":"Fosc")+":"+cid.substr(cid.find(':')+1); // behaviour = kind only; each size/shape is a new shape
        h.lanes=c.count; h.setup=normCells(c.ex);
        h.text=objName(c.I,c.ex.size())+" in the ash, first in soup seed "+to_string(c.seed); h.credit=0.3; by[cid]=h; }
    vector<ProbeHit> out; for(auto&kv:by) if(kv.second.lanes<=max(1L,(long)(C.fossilRate*E.n))) out.push_back(kv.second);
    sort(out.begin(),out.end(),[](const ProbeHit&a,const ProbeHit&b){ return a.setup.size()!=b.setup.size()? a.setup.size()>b.setup.size() : a.id<b.id; });
    if((int)out.size()>C.fossilMax) out.resize(C.fossilMax); return out;
}

// report: new ids go to the archive/log; returns novelty credit (behaviour = id up to ':')
static double recordProbe(const Eval&E,const vector<ProbeHit>&hits,const string&what,const string&ctx,bool printAll){
    double credit=0; string rule=ruleString(E.T);
    for(auto&h:hits){
        bool known=archiveIds.count(h.id); if(known&&!printAll) continue;
        bool nd=false; if(!known){ nd=descCount[descOf(h.id)]++==0; archiveIds.insert(h.id); credit+= nd? h.credit : 0.1*h.credit;
            if(!C.archiveFile.empty()){ FILE*a=fopen(C.archiveFile.c_str(),"a"); if(a){ fprintf(a,"%s\n",h.id.c_str()); fclose(a);} } }
        string msg=string(nd?"# NEW ":"# ")+what+": "+h.text+", "+to_string(h.lanes)+(what=="fossil"?" in "+to_string(E.n)+" soups":" setups")+" ("+ctx+")\n"+toRLE(h.setup,rule);
        printf("%s",msg.c_str()); fflush(stdout); if(logFp&&!known){ fprintf(logFp,"%s\n",msg.c_str()); fflush(logFp); }
    }
    return credit;
}
// replicators: a pattern that, on its own, makes exact copies of itself (in any orientation of the symmetry group).
// Stage 1 runs random small seeds and some census soups and keeps cluster shapes that occur twice at once and, alone,
// copy themselves within repGens/4 gens; stage 2 runs each such
// shape alone: it must not come back alone as itself (oscillator/ship), must reach >= 3 exact copies, and when the
// second copy appears at least half of all live cells must be exact copies.
// Controllability = the best share of live cells that are exact copies around gens 4R..8R (R = first time there are
// two copies): 1 for clean doubling like HighLife's replicator, near 0 when the copies dissolve into chaos.
static vector<ProbeHit> probeReplicators(const Eval&E){
    vector<ProbeHit> out; if(C.repSeeds<=0&&C.repAdd.empty()) return out; probeSetRule(E.T);
    auto canon=[](const Cells&c){ string best; for(int o:groupOrients()){ string k=keyOf(normCells(xform(c,o))); if(best.empty()||k<best) best=k; } return best; };
    // 8-connected clusters of the live cells within radius r of the grid centre, by flood fill on the grid itself;
    // stops (returns false) past cap live cells. n = number of live cells seen.
    auto clusters=[](Engine&G,int r,size_t cap,vector<Cells>&out,size_t&n){ out.clear(); n=0; static thread_local vector<uint8_t> vis; static thread_local vector<int> touched;
        if(vis.size()!=(size_t)G.W*G.H) vis.assign((size_t)G.W*G.H,0); int x0=max(1,G.W/2-r),x1=min(G.W-2,G.W/2+r),y0=max(1,G.H/2-r),y1=min(G.H-2,G.H/2+r); bool ok=true; vector<int> st;
        for(int y=y0;y<=y1&&ok;y++) for(int x=x0;x<=x1;x++){ int i=y*G.W+x; if(!G.a[i]||vis[i]) continue; Cells cc; st.assign(1,i); vis[i]=1; touched.push_back(i);
            while(!st.empty()){ int j=st.back(); st.pop_back(); cc.push_back({j%G.W,j/G.W});
                for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){ int k=j+dy*G.W+dx; if(k>=0&&k<(int)vis.size()&&G.a[k]&&!vis[k]){ vis[k]=1; touched.push_back(k); st.push_back(k); } } }
            n+=cc.size(); out.push_back(std::move(cc)); if(n>cap){ ok=false; break; } }
        for(int i:touched) vis[i]=0; touched.clear(); return ok; };
    map<string,pair<Cells,int>> cand;
    #pragma omp parallel
    { Engine&Gn=*probeEngines[omp_get_thread_num()]; map<string,bool> qt; // quick-test results, shared by this thread's runs
      #pragma omp for schedule(dynamic,8)
      for(int s=0;s<C.repSeeds+C.repSoups;s++){ Cells seed; int gens=C.repGens/2, rad0=8;
          if(s<C.repSeeds){ u64 r=0x5EED0000ULL+s; int b=3+s%4; for(int y=0;y<b;y++){ u64 z=splitmix(r); for(int x=0;x<b;x++) if(z>>x&1) seed.push_back({x,y}); } }
          else { int fw,fh; seed=makeSoup(soupSeed(E.seedEpoch,s-C.repSeeds),fw,fh); gens=C.repGens; rad0=max(fw,fh); }
          if(seed.empty()) continue; int ox,oy; map<string,Cells> found;
          // quick test: on its own the shape makes two copies of itself within repGens/4 gens
          auto copiesItself=[&](const Cells&cc,const string&k){ Cells x=cc; for(int g=1;g<=C.repGens/4;g++){ x=stepCells(x,E.T.t); if(x.empty()||x.size()>500) return false;
              if(x.size()<2*cc.size()) continue; int n=0; for(auto&q:comps8(x)) if(q.size()==cc.size()&&canon(q)==k&&++n==2) return true; } return false; };
          int every= s<C.repSeeds? 4 : 8; size_t cap=1500;
          vector<Cells> cs; size_t npop;
          Gn.runPattern(seed,gens,ox,oy,[&](int t){ if(t<8||t%every) return true; if(!clusters(Gn,t+rad0,cap,cs,npop)) return false;
              map<long,int> sg; vector<long> sgs(cs.size(),-1); // cheap signature first: size and bounding box
              for(size_t i=0;i<cs.size();i++){ auto&cc=cs[i]; if(cc.size()<4||cc.size()>300) continue; int mx,my,Mx,My; bbox(cc,mx,my,Mx,My); int w=Mx-mx+1,h=My-my+1;
                  sgs[i]=(long)cc.size()*1000000+min(w,h)*1000+max(w,h); sg[sgs[i]]++; }
              map<string,int> cnt; for(size_t i=0;i<cs.size();i++){ auto&cc=cs[i]; if(sgs[i]<0||sg[sgs[i]]<2) continue; string k=canon(cc);
                  if(++cnt[k]==2&&!found.count(k)){ auto it=qt.find(k); bool ok= it!=qt.end()? it->second : (qt[k]=copiesItself(cc,k)); if(!ok) continue; found[k]=cc; if(found.size()>=8) return false; } }
              return true; });
          if(!found.empty()){
              #pragma omp critical(repcand)
              for(auto&kv:found){ auto&e=cand[kv.first]; if(!e.second) e.first=normCells(kv.second); e.second++; } } } }
    for(auto&r:C.repAdd){ Cells c=normCells(parseRLE(r)); if(!c.empty()){ auto&e=cand[canon(c)]; if(!e.second) e.first=c; e.second+=1000000; } } // replicator_add RLE: always tested
    vector<pair<string,pair<Cells,int>>> cv(cand.begin(),cand.end());
    sort(cv.begin(),cv.end(),[](auto&a,auto&b){ return a.second.second!=b.second.second? a.second.second>b.second.second : a.first<b.first; });
    if((int)cv.size()>C.repMax) cv.resize(C.repMax);
    vector<vector<ProbeHit>> res(cv.size());
    #pragma omp parallel
    { Engine&Gn=*probeEngines[omp_get_thread_num()];
      #pragma omp for schedule(dynamic,1)
      for(size_t q=0;q<cv.size();q++){ const string&key=cv[q].first; const Cells&K=cv[q].second.first; int mx,my,Mx,My; bbox(K,mx,my,Mx,My); int rad=max(Mx-mx,My-my)+4;
          int R=0,maxC=0,last=0; bool periodic=false; vector<double> fr(C.repGens+1,-1); int ox,oy;
          vector<Cells> cs; size_t npop;
          Gn.runPattern(K,C.repGens,ox,oy,[&](int t){ last=t; if(!clusters(Gn,t+rad,20000,cs,npop)||npop==0) return false;
              int copies=0; for(auto&cc:cs) if(cc.size()==K.size()&&canon(cc)==key) copies++;
              if(copies==1&&npop==K.size()){ periodic=true; return false; }
              if(copies>=2){ if(!R) R=t; fr[t]=(double)copies*K.size()/npop; } maxC=max(maxC,copies); return true; });
          if(periodic||!R||maxC<3||fr[R]<0.5) continue; // first copies must be clean: chance repeats of a shape in chaos are not
          int a=min(4*R,last), e=min(8*R,last); if(a>e-R) a=max(2*R,e-R); double f=0; for(int t=a;t<=e;t++) f=max(f,fr[t]);
          char buf[256]; snprintf(buf,256,"%d cells, 2 copies at gen %d, up to %d copies by gen %d, controllability %.2f",(int)K.size(),R,maxC,last,f);
          ProbeHit h; h.id="Rrep:"; for(unsigned char ch:key){ char b2[3]; snprintf(b2,3,"%02x",ch); h.id+=b2; }
          h.text=buf; h.setup=K; h.lanes=cv[q].second.second%1000000; h.credit=4+6*f; h.score=f; h.sig=R*1000+maxC; res[q].push_back(h); } }
    for(auto&v:res) for(auto&h:v) out.push_back(h);
    sort(out.begin(),out.end(),[](const ProbeHit&a,const ProbeHit&b){ return a.score!=b.score? a.score>b.score : a.id<b.id; });
    return out;
}
static double probeRule(const Eval&E,const string&ctx,bool printAll){
    double t0=nowS(); vector<LibObj> ships,stills; probeLibrary(E,ships,stills);
    if(printAll){ printf("probe library: %zu ships, %zu still lifes\n",ships.size(),stills.size());
        for(auto&s:ships) printf("  ship  %6ld  %s\n",s.count,objName(s.I,s.c.size()).c_str());
        for(auto&s:stills) printf("  still %6ld  %s\n",s.count,objName(s.I,s.c.size()).c_str()); }
    double cr=0; auto cat=probeCatalysts(E,ships,stills); cr+=recordProbe(E,cat,"catalyst",ctx,printAll);
    int nw=0; auto sig=probeWicks(E,&nw); cr+=recordProbe(E,sig,"wick",ctx,printAll);
    auto fos=probeFossils(E); cr+=recordProbe(E,fos,"fossil",ctx,printAll);
    auto rep=probeReplicators(E); cr+=recordProbe(E,rep,"replicator",ctx,printAll);
    if(printAll) printf("probe: %zu catalytic reactions, %d wicks tested -> %zu signals/fuses, %zu fossils, %zu replicators; %.1fs\n",cat.size(),nw,sig.size(),fos.size(),rep.size(),nowS()-t0);
    return cr;
}
static void modeProbe(const string&rule,int N,u64 seed){
    Table T=parseRule(rule); if(!satisfies(T)) printf("warning: rule violates constraints\n");
    C.N=N; for(auto*g:engines) g->collectAsh=true; Eval E; fullEval(E,T,seed); printEval("census",E);
    recordFinds(E,processCands(E),"census"); recordSoupEvents(E,"census"); // ships/oscillators/puffers too (all shapes go to the log)
    probeRule(E,"probe",true);
}
static void loadArchive(){ if(C.archiveFile.empty()) return; ifstream f(C.archiveFile); string s; while(f>>s){ if(archiveIds.insert(s).second) descCount[descOf(s)]++; } }

// ===================== mutation =====================
static vector<int> freeClasses(){ vector<int> v; for(int c=0;c<NC;c++) if(forcedV[c]<0) v.push_back(c); return v; }
static Table mutate(const Table&T,const vector<double>&usage,u64&rng,string&desc){
    static vector<int> fc=freeClasses();
    vector<double> w(fc.size()); double tot=0, tu=0; for(size_t i=0;i<fc.size();i++) tu+=sqrt(usage.empty()?0:usage[fc[i]]);
    double eps=max(1e-9,0.02*tu/fc.size());
    for(size_t i=0;i<fc.size();i++){ w[i]=sqrt(usage.empty()?0:usage[fc[i]])+eps; tot+=w[i]; }
    double r=(splitmix(rng)>>11)*(1.0/9007199254740992.0); int k= r<0.5?1: r<0.8?2:3;
    Table M=T; desc.clear(); set<int> chosen;
    auto pick=[&](bool usedOnly)->int{ double tt=0; for(size_t i=0;i<fc.size();i++) if(!usedOnly||usage[fc[i]]>0) tt+=w[i];
        double u=(splitmix(rng)>>11)*(1.0/9007199254740992.0)*tt;
        for(size_t i=0;i<fc.size();i++){ if(usedOnly&&!(usage[fc[i]]>0)) continue; if(u<w[i]) return fc[i]; u-=w[i]; }
        return fc.back(); };
    bool anyUsed=false; for(int c:fc) if(!usage.empty()&&usage[c]>0) anyUsed=true;
    int guard=0; while((int)chosen.size()<min(k,(int)fc.size()) && guard++<1000) chosen.insert(pick(false));
    bool has=false; for(int c:chosen) if(!usage.empty()&&usage[c]>0) has=true;
    if(anyUsed&&!has){ chosen.erase(chosen.begin()); int g2=0; int c; do c=pick(true); while(chosen.count(c)&&g2++<1000); chosen.insert(c); }
    for(int c:chosen){ int v=!M.t[members[c][0]]; for(int y:members[c]) M.t[y]=v; desc+= string(v?"+":"-")+cname[c]+" "; }
    return M;
}

// paired race: evaluate mutant in batches, reject early when clearly worse than parent (linearised objective)
static bool race(Eval&M,const Eval&P,int batch){
    int nt=C.tg.size(); vector<double> d;
    for(int i0=0;i0<C.N;){
        int i1=min(C.N, i0==0 ? min(32,batch) : i0+batch); evalRange(M,i0,i1);
        for(int i=i0;i<i1;i++){
            auto g=[&](const SoupRes&r,bool par){ double s=0; if(!r.bad) for(int t=0;t<nt;t++) s+=C.tg[t].w*min((int)r.cnt[t],C.cap)/(P.rate[t]+1e-3);
                s-=5*(r.bad?1:0)*(P.bad>C.maxBad?1:0.2); if(C.costW>0&&C.costRef>0) s-=C.costW*r.cost/1e6/max(P.cost,1e-9); return s; };
            d.push_back(g(M.res[i].r,false)-g(P.res[i].r,true));
        }
        i0=i1; if(i1>=C.N) break;
        { double cs=0,bd=0; for(int i=0;i<i1;i++){ cs+=M.res[i].r.cost; bd+=M.res[i].r.bad?1:0; } cs/=i1*1e6; bd/=i1;
          if(C.maxCost>0 && cs>1.5*C.maxCost && cs>P.cost) return false;
          if(bd>C.maxBad+0.1 && bd>P.bad+0.05) return false; }
        if(d.size()>=(size_t)batch*2){ double m=0,v=0; for(double x:d) m+=x; m/=d.size(); for(double x:d) v+=(x-m)*(x-m); v/=(d.size()-1);
            if(m+2.5*sqrt(v/d.size())<0) return false; }
    }
    return true;
}

static volatile sig_atomic_t stopFlag=0;
static void onSigint(int){ if(stopFlag) _exit(130); stopFlag=1; const char m[]="\n[interrupt: finishing this iteration; Ctrl-C again to quit immediately]\n"; if(write(2,m,sizeof m-1)<0){} }
static string evalLine(const char*tag,const Eval&E);
static void writeBest(const vector<const Eval*>&es,const char*what){
    if(C.bestFile.empty()) return; FILE*f=fopen(C.bestFile.c_str(),"w"); if(!f) return;
    fprintf(f,"# %s\n",what); for(auto*e:es) fprintf(f,"%s",evalLine("#",*e).c_str()), fprintf(f,"%s\n",ruleString(e->T).c_str()); fclose(f); }
static void printEval(const char*tag,const Eval&E){ printf("%s",evalLine(tag,E).c_str()); fflush(stdout); }
static string evalLine(const char*tag,const Eval&E){
    string o; char b[256]; snprintf(b,256,"%s f=%.4f",tag,E.f); o+=b; for(size_t t=0;t<C.tg.size();t++){ snprintf(b,256,E.rate[t]>=0.01||E.rate[t]==0?" %s=%.3f":" %s=%.2e",C.tg[t].name.c_str(),E.rate[t]); o+=b; }
    snprintf(b,256," bad=%.3f cost=%.3fM %s ",E.bad,E.cost,feasible(E)?"feasible":"INFEASIBLE"); o+=b; return o+ruleString(E.T)+"\n"; }

// ===================== modes =====================
static void modeCheck(const vector<string>&rules){
    printf("symmetry %s: %d classes, %d forced, %d free\n",SYM==0?"totalistic":SYM==1?"isotropic":"none",NC,(int)count_if(forcedV.begin(),forcedV.end(),[](int v){return v>=0;}),(int)freeClasses().size());
    for(auto&P:C.pats){ auto&E=P.evo; printf("pattern %-10s %3zu cells: ",P.name.c_str(),P.cells.size());
        if(E.s==-1) printf("dies at gen %d\n",E.t); else printf("gen %d repeats gen %d: period %d, shift (%d,%d)\n",E.t,E.s,E.t-E.s,E.dx,E.dy); }
    for(auto&T:C.tg) printf("target %s: %s period %d shift (%d,%d), %zu phase keys, weight %.2f, min_rate %.3f\n",T.name.c_str(),(T.dx||T.dy)?"ship":"oscillator",T.period,T.dx,T.dy,T.keys.size(),T.w,T.minRate);
    if(SYM!=2){ printf("forced:"); for(int c=0;c<NC;c++) if(forcedV[c]>=0) printf(" %s=%d",cname[c].c_str(),forcedV[c]); printf("\n"); }
    else printf("forced: %d of 512 indices\n",(int)count_if(forcedV.begin(),forcedV.end(),[](int v){return v>=0;}));
    for(auto&r:rules){ Table T=parseRule(r); printf("%s: %s, symmetric=%s, canonical %s\n",r.c_str(),satisfies(T)?"satisfies constraints":"VIOLATES constraints",respectsSym(T)?"yes":"no",ruleString(T).c_str()); }
}
static void modeCensus(const string&rule,int N,u64 seed){
    Table T=parseRule(rule); if(!satisfies(T)) printf("warning: rule violates constraints\n");
    C.N=N; Eval E; auto t0=chrono::steady_clock::now(); fullEval(E,T,seed);
    printEval("census",E); printf("time %.2fs\n",chrono::duration<double>(chrono::steady_clock::now()-t0).count());
    if(getenv("BADSTAT")){ int nb[4]={0,0,0,0}; for(int i=0;i<E.n;i++) nb[E.res[i].r.bad]++; printf("bad types: exploded %d, unsettled %d, too slow %d\n",nb[1],nb[2],nb[3]); }
    for(auto&e:soupEvents(E)) printf("%6ld  %s  first in soup seed %llu (dump RULE SEED shows it)\n",e.count,e.desc.c_str(),(unsigned long long)e.seed);
    auto fs=processCands(E);
    for(auto&F:fs){ printf("%6ld  %-34s %3zu cells  %s",F.count,describe(F.I).c_str(),F.phase.size(),archiveIds.count(F.id)?"":"(not in archive) "); string r=toRLE(F.phase,ruleString(T)); printf("\n%s",r.c_str()); }
}
static void modeDump(const string&rule,u64 seed){
    Table T=parseRule(rule); int fw,fh; Cells c=makeSoup(seed,fw,fh);
    printf("%s",toRLE(c,ruleString(T)).c_str());
}
static Table startRule(int i){ return project(parseRule(C.starts.empty()? C.refRule : C.starts[i%C.starts.size()])); }

static void modeClimb(int iters){
    u64 rng=C.seed*7919+17, epoch=C.seed*1000+1; Eval cur; fullEval(cur,startRule(0),epoch);
    if(C.costW>0&&C.costRef<=0){ C.costRef=cur.cost; score(cur); printf("cost_ref set to start cost %.3fM\n",C.costRef); }
    printEval("start",cur); recordFinds(cur,processCands(cur),"start"); recordSoupEvents(cur,"start");
    Eval best=cur; vector<double> hist(max(1,C.lahc),cur.f); long evals=0, early=0;
    double tStart=nowS();
    for(int it=1;it<=iters&&!stopFlag;it++){
        if(it%50==0){ printf("[it %d: %ld evals, %ld rejected early, current f %.4f, best f %.4f; time %.0fs: sim %.0fs, objects %.0fs]\n",it,evals,early,cur.f,best.f,nowS()-tStart,tSim,tProc); fflush(stdout); }
        if(C.reseed>0&&it%C.reseed==0){ epoch++; Table T=cur.T; fullEval(cur,T,epoch); fill(hist.begin(),hist.end(),cur.f); printEval(("reseed@"+to_string(it)).c_str(),cur); }
        string desc; Table MT=mutate(cur.T,usageOf(cur),rng,desc);
        Eval M; beginEval(M,MT,epoch); evals++;
        if(!race(M,cur,max(32,C.N/8))){ early++; continue; }
        finishEval(M); score(M);
        double h=hist[it%hist.size()];
        if(M.f>=cur.f||M.f>=h){
            bool improved=M.f>cur.f; cur=std::move(M); hist[it%hist.size()]=cur.f;
            if(improved){ printEval(("it "+to_string(it)+" ["+desc+"]").c_str(),cur); }
            if(feasible(cur)) recordFinds(cur,processCands(cur),"climb it "+to_string(it));
            recordSoupEvents(cur,"climb it "+to_string(it));
            // best = highest f among feasible rules (any rule only while none is feasible)
            if((feasible(cur)&&(!feasible(best)||cur.f>best.f)) || (!feasible(best)&&cur.f>best.f)){ best=cur; writeBest({&best},"climb: best rule so far"); }
        } else hist[it%hist.size()]=cur.f;
    }
    printEval("best",best); writeBest({&best},"climb: best rule");
}

struct Member { Eval E; double credit=0; int visits=0; vector<string> finds; };
// ===================== variety: how many different interesting things one rule has =====================
// Each distinct behaviour scores weight * presence, presence = min(1, count / (vary_min * soups)), so reliable
// behaviours count fully and one-off flukes little. Behaviours: ship velocities (oblique worth more), oscillator
// periods >= min_osc_period, extra shapes of them, still-life diversity, puffer/gun/active-region soups, and (with
// the probe) catalyst kinds, signals and fuses. Returns the score; 'why' gets a breakdown.
struct Variety { double score=0, ships=0, oscs=0, shapes=0, stills=0, soup=0, cat=0, sig=0, ca=0, rep=0; int nShip=0,nOsc=0; bool confined=false; };
static Variety varietyOf(const Eval&E,bool probe,const string&ctx,bool log){
    Variety V; double need=max(1.0,C.varyMin*E.n); auto pres=[&](double c){ return min(1.0,c/need); };
    map<string,pair<double,int>> beh; // behaviour -> (count, shapes)
    for(auto&F:processCands(E)) if(F.I.type==2||F.I.period>=C.minOscPeriod){ auto&b=beh[descOf(F.id)]; b.first+=F.count; b.second++; }
    for(size_t t=0;t<C.tg.size();t++){ auto&T=C.tg[t]; Info I; I.type=(T.dx||T.dy)?2:1; I.period=T.period; I.dx=T.dx; I.dy=T.dy;
        Cells c=C.pats[T.pat].evo.gens[C.pats[T.pat].evo.s]; auto&b=beh[descOf(canonId(c,I,E.T.t))]; b.first+=E.rate[t]*E.n; b.second++; }
    // vary_strict N: an oscillator period seen in fewer than N soups scores nothing (one-off periods are mostly random
    // geometry: signals in walled boxes, chaotic agar fragments). vary_decay d: the k-th best period, with its shape
    // bonus, is worth d^k.
    vector<pair<double,double>> oscTerms; // (period term, shape bonus), only with vary_decay
    for(auto&kv:beh){ const string&d=kv.first; double p=pres(kv.second.first); int P=atoi(d.c_str()+1); double shp;
        if(d[0]=='Q'){ int ax=0,ay=0; sscanf(d.c_str(),"Q%*d(%d,%d)",&ax,&ay); double w= (ax&&ay&&ax!=ay)? 6 : 3;
            if(C.varyCommonShip>0) w*=1+C.varyCommonShip*min(1.0,kv.second.first/(C.varyCommonRate*E.n)); // vary_common_ship m r
            V.ships+=w*p; if(p>=1) V.nShip++; }
        else { if(C.varyStrict>0&&kv.second.first<C.varyStrict) p=0; if(p>=1) V.nOsc++;
            double t=(1+log2((double)max(1,P))/2)*p;
            if(C.varyDecay!=1.0){ oscTerms.push_back({t,0.2*min(5,kv.second.second-1)*p}); continue; }
            V.oscs+=t; }
        shp=0.2*min(5,kv.second.second-1)*p; V.shapes+=shp; }
    sort(oscTerms.rbegin(),oscTerms.rend()); { double w=1; for(auto&t:oscTerms){ V.oscs+=w*t.first; V.shapes+=w*t.second; w*=C.varyDecay; } }
    { set<string> sl; for(auto&kv:E.ash) if(kv.second.I.type==1&&kv.second.I.period==1&&kv.second.count>=need) sl.insert(canonId(normCells(kv.second.ex),kv.second.I,E.T.t)); V.stills=0.5*log2(1.0+sl.size()); }
    { map<string,double> sb; for(auto&e:soupEvents(E)){ string k=e.id.substr(0,1)+(e.id[0]=='G'? e.id.substr(1,e.id.find(':')-1):""); sb[k]=max(sb[k],e.credit/ (e.id[0]=='U'?10.0:e.id[0]=='G'?3.0:2.0)); }
      int ng=0; for(auto&kv:sb){ char c=kv.first[0]; if(c=='G'){ if(ng++<5) V.soup+=2; } else if(c=='U') V.soup+=4; else if(c=='A') V.soup+=2; }
      bool moves=false; for(auto&kv:beh) if(kv.first[0]=='Q') moves=true; for(auto&kv:sb) if(kv.first[0]=='G'||kv.first[0]=='U') moves=true; V.confined=!moves; }
    if(probe){ vector<LibObj> ships,stills; probeLibrary(E,ships,stills);
        auto cat=probeCatalysts(E,ships,stills); set<string> kinds;
        for(auto&h:cat){ string k=h.id.substr(1,h.id.find('_')-1); string ship=h.id.substr(h.id.find('_')+1); ship=ship.substr(0,ship.find(':'));
            if(k.rfind("reflector",0)==0) k="reflector"; if(k.rfind("converter",0)==0) k="converter"; if(k.rfind("splitter",0)==0) k="splitter"; kinds.insert(k+"|"+ship); }
        for(auto&k:kinds){ string c=k.substr(0,k.find('|')); V.cat+= c=="reflector"?5: c=="converter"||c=="splitter"?2: c=="catalysed"?0.5:1; }
        auto sig=probeWicks(E,nullptr); set<string> sg,fu,ca; for(auto&h:sig){ string d=descOf(h.id); (d[1]=='s'?sg:d[1]=='1'?ca:fu).insert(d); }
        V.sig=4*sg.size()+min<size_t>(3,fu.size()); V.ca= ca.empty()? 0 : 4+min<size_t>(2,ca.size()-1); // one behaviour: a 1D CA wick usually shows several speeds
        vector<ProbeHit> rp; if(C.varyRep) rp=probeReplicators(E); if(!rp.empty()){ set<long> kinds; for(auto&h:rp) kinds.insert(h.sig); // phases of one replicator share R and copy count
            V.rep=4+6*rp[0].score+min<size_t>(2,kinds.size()-1); } // best one counts by controllability
        if(log){ recordProbe(E,cat,"catalyst",ctx,false); recordProbe(E,sig,"wick",ctx,false); recordProbe(E,rp,"replicator",ctx,false); } }
    V.score=V.ships+V.oscs+V.shapes+V.stills+V.soup+V.cat+V.sig+V.ca+V.rep;
    if(V.rep>0) V.confined=false; // a replicator spreads
    if(V.confined&&C.varyConfined!=1.0) V.score*=C.varyConfined; // vary_confined f: no ship and no growth at all, everything stays in its box
    return V;
}
static string varietyLine(const Variety&V){ char b[300]; snprintf(b,256,"variety %.2f [ships %.1f (%d), oscillators %.1f (%d), shapes %.1f, still lifes %.1f, soup %.1f, catalysts %.1f, signals %.1f, 1D CA %.1f, replicators %.1f]",
    V.score,V.ships,V.nShip,V.oscs,V.nOsc,V.shapes,V.stills,V.soup,V.cat,V.sig,V.ca,V.rep); string r=b; if(V.confined&&C.varyConfined!=1.0) r+=" (confined: x"+to_string(C.varyConfined).substr(0,4)+")"; return r; }
// vary: late-acceptance hill climb on the variety score (constraints as usual; targets/min_rate/max_bad/max_cost
// only as feasibility). No racing: a mutant is dropped after the first 1/8 of the soups if it is clearly too bad.
static void modeVary(int iters){
    for(auto*g:engines) g->collectAsh=true;
    bool probe=C.varyProbe; u64 rng=C.seed*15485863+5, epoch=C.seed*1000+9; double tStart=nowS();
    auto value=[&](Eval&E,Variety&V,const string&ctx){ V=varietyOf(E,probe,ctx,true); recordFinds(E,processCands(E),ctx); recordSoupEvents(E,ctx);
        return feasible(E)? V.score : V.score-10-20*max(0.0,E.bad-C.maxBad); };
    Eval cur; fullEval(cur,startRule(0),epoch); Variety cv; double cf=value(cur,cv,"start");
    printf("start %s\n      ",varietyLine(cv).c_str()); printEval("",cur);
    Eval best=cur; Variety bv=cv; double bf=cf; vector<double> hist(max(1,C.lahc),cf); long evals=0,early=0;
    for(int it=1;it<=iters&&!stopFlag;it++){
        if(it%25==0){ printf("[it %d: %ld evals, %ld rejected early, current %.2f, best %.2f, %.0fs]\n",it,evals,early,cf,bf,nowS()-tStart); fflush(stdout); }
        if(C.reseed>0&&it%C.reseed==0){ epoch++; Table T=cur.T; fullEval(cur,T,epoch); cf=value(cur,cv,"reseed"); fill(hist.begin(),hist.end(),cf); printf("reseed@%d %s\n",it,varietyLine(cv).c_str()); }
        string desc; Table MT=mutate(cur.T,usageOf(cur),rng,desc);
        Eval M; beginEval(M,MT,epoch); evals++;
        int b=max(32,C.N/8); evalRange(M,0,b); score(M);
        if(M.bad>2*C.maxBad+0.05 || (C.maxCost>0&&M.cost>2*C.maxCost)){ early++; continue; }
        evalRange(M,b,C.N); finishEval(M); score(M);
        Variety mv; double mf=value(M,mv,"vary it "+to_string(it)+" ["+desc+"]");
        double h=hist[it%hist.size()];
        if(mf>=cf||mf>=h){ bool imp=mf>cf; cur=std::move(M); cf=mf; cv=mv; hist[it%hist.size()]=cf;
            if(imp){ printf("it %d [%s] %s\n      ",it,desc.c_str(),varietyLine(cv).c_str()); printEval("",cur); }
            if(cf>bf&&feasible(cur)){ best=cur; bv=cv; bf=cf; writeBest({&best},("vary: most varied rule so far, "+varietyLine(bv)).c_str()); } }
        else hist[it%hist.size()]=cf;
    }
    printf("best %s\n     ",varietyLine(bv).c_str()); printEval("",best); writeBest({&best},("vary: most varied rule, "+varietyLine(bv)).c_str());
}
static void modeExplore(int iters,int cap){
    if(C.probe) for(auto*g:engines) g->collectAsh=true;
    u64 rng=C.seed*104729+3, epoch=C.seed*1000+7;
    vector<Member> pop;
    for(size_t i=0;i<max<size_t>(1,C.starts.size());i++){ Member m; fullEval(m.E,startRule(i),epoch);
        if(C.costW>0&&C.costRef<=0){ C.costRef=m.E.cost; score(m.E); printf("cost_ref set to start cost %.3fM\n",C.costRef); }
        printEval("start",m.E); auto fs=processCands(m.E); m.credit=recordFinds(m.E,fs,"start")+recordSoupEvents(m.E,"start"); if(C.probe) m.credit+=probeRule(m.E,"start",false); pop.push_back(std::move(m)); }
    long evals=0,early=0,infeas=0,added=0; double tStart=nowS();
    for(int it=1;it<=iters&&!stopFlag;it++){
        // parent: weight (1+credit)/sqrt(1+visits) * feasibility
        vector<double> w(pop.size()); double tot=0; for(size_t i=0;i<pop.size();i++){ w[i]=(1+pop[i].credit)/sqrt(1.0+pop[i].visits)*(feasible(pop[i].E)?1:0.2); tot+=w[i]; }
        double u=(splitmix(rng)>>11)*(1.0/9007199254740992.0)*tot; size_t pi=0; while(pi+1<pop.size()&&u>w[pi]){u-=w[pi];pi++;}
        Member&P=pop[pi]; P.visits++;
        string desc; Table MT=mutate(P.E.T,usageOf(P.E),rng,desc);
        Member m; beginEval(m.E,MT,epoch); evals++;
        // quick feasibility check on first batch
        int b=max(32,C.N/8); evalRange(m.E,0,32); score(m.E);
        if(m.E.bad>0.25 || (C.maxCost>0&&m.E.cost>3*C.maxCost)){ early++; continue; }
        evalRange(m.E,32,b); score(m.E);
        bool hopeless=m.E.bad>2*C.maxBad+0.05; for(size_t t=0;t<C.tg.size();t++) hopeless|= m.E.rate[t]<0.6*C.tg[t].minRate; if(C.maxCost>0) hopeless|=m.E.cost>1.5*C.maxCost;
        if(hopeless){ early++; continue; }
        evalRange(m.E,b,C.N); finishEval(m.E); score(m.E);
        double sc=recordSoupEvents(m.E,"explore it "+to_string(it)+", parent ["+desc+"]");
        if(!feasible(m.E)){ infeas++; continue; }
        auto fs=processCands(m.E);
        if(C.probe) sc+=probeRule(m.E,"explore it "+to_string(it)+", parent ["+desc+"]",false);
        m.credit=sc+recordFinds(m.E,fs,"explore it "+to_string(it)+", parent ["+desc+"]");
        for(auto&F:fs) m.finds.push_back(F.id);
        if(m.credit>0 || m.E.f>=P.E.f-0.1){
            P.credit*=0.9; // parent shares a little decay so the search keeps moving
            pop.push_back(std::move(m)); added++;
            if((int)pop.size()>cap){ size_t worst=0; double ws=1e18; for(size_t i=0;i<pop.size()-1;i++){ double s=pop[i].E.f+pop[i].credit-0.05*pop[i].visits; if(s<ws){ws=s;worst=i;} } pop.erase(pop.begin()+worst); }
        }
        if(it%25==0){ double bf=-1e9; for(auto&q:pop) bf=max(bf,q.E.f); printf("[it %d: pop %zu, %ld evals, %ld early-rejected, %ld infeasible, %ld added, archive %zu, best f %.3f, %.0fs]\n",it,pop.size(),evals,early,infeas,added,archiveIds.size(),bf,nowS()-tStart); fflush(stdout);
            vector<const Eval*> es; for(auto&q:pop) es.push_back(&q.E); sort(es.begin(),es.end(),[](auto*a,auto*b){return a->f>b->f;}); writeBest(es,"explore: current population, best f first"); }
    }
    printf("final population:\n"); sort(pop.begin(),pop.end(),[](const Member&a,const Member&b){return a.E.f>b.E.f;});
    for(auto&q:pop) printEval(("  credit "+to_string(q.credit).substr(0,5)).c_str(),q.E);
    vector<const Eval*> es; for(auto&q:pop) es.push_back(&q.E); writeBest(es,"explore: final population, best f first");
}

static void selfTest(){
    // Life via Hensel and totalistic parse agree; glider moves (1,1) per 4 gens; MAP round trip
    Table L=parseRule("B3/S23"); Cells g=parseRLE("bo$2bo$3o!"); auto E=evolve(g,L.t,10,nullptr);
    if(!(E.s==0&&E.t==4&&abs(E.dx)==1&&abs(E.dy)==1)) die("selftest: glider");
    if(!(parseRule(mapString(L))==L)) die("selftest: MAP round trip");
    if(ruleString(L)!="B3/S23") die("selftest: life string "+ruleString(L));
    { u64 s=7; for(int k=0;k<300;k++){ Cells c; int w=1+splitmix(s)%20,h=1+splitmix(s)%20; for(int y=0;y<h;y++) for(int x=0;x<w;x++) if(splitmix(s)%5==0) c.push_back({x,y}); if(c.empty()) continue;
        string r=toRLE(c,"B3/S23"); string body=r.substr(r.find('\n')+1); Cells d=parseRLE(body);
        if(keyOf(d)!=keyOf(c)) die("selftest: RLE round trip\n"+r); } }
    // Golly's Life MAP string (Golly Help: MAP rules) must parse to Life
    u64 s=99; for(int k=0;k<200;k++){ Table T; for(int i=0;i<102;i++){ int v=splitmix(s)&1; for(int x=0;x<512;x++) if(hclass[x]==i) T.t[x]=v; } T.t[0]=0;
        Table T2=parseRule(ruleString(T)); if(!(T2==T)) die("selftest: hensel round trip "+ruleString(T)); if(!(parseRule(mapString(T))==T)) die("selftest: map rt"); }
}

int main(int argc,char**argv){
    initHensel(); initP2(); selfTest();
    // key=value arguments anywhere after CONFIG override config lines (commas stand for spaces: soup=16,16,0.4)
    vector<string> args, ov; for(int i=1;i<argc;i++){ string a=argv[i]; size_t e=a.find('='); if(i>1&&e!=string::npos&&e>0&&isalpha((unsigned char)a[0])&&a.compare(0,3,"MAP")!=0){ string v=a.substr(e+1); replace(v.begin(),v.end(),',',' '); ov.push_back(a.substr(0,e)+" "+v); } else args.push_back(a); }
    argc=args.size()+1; vector<char*> av={argv[0]}; for(auto&a:args) av.push_back((char*)a.c_str()); argv=av.data();
    if(argc<3){ fprintf(stderr,"usage: rulesearch CONFIG MODE [ARGS] [key=value ...]\n  modes: check [RULE...] | census RULE [N] [SEED] | climb ITERS | explore ITERS [POPCAP] | vary ITERS | probe RULE [N] [SEED] | dump RULE SEED\n  key=value overrides a config line, e.g. seed=3 soups=400 soup=16,16,0.4\n"); return 1; }
    loadCfg(argv[1],ov); string mode=argv[2];
    signal(SIGINT,onSigint);
    int nt=omp_get_max_threads(); for(int i=0;i<nt;i++) engines.push_back(new Engine(C.grid));
    loadArchive(); if(!C.logFile.empty()) logFp=fopen(C.logFile.c_str(),"a");
    if(mode=="check"){ vector<string> r; for(int i=3;i<argc;i++) r.push_back(argv[i]); modeCheck(r); }
    else if(mode=="census"){ if(argc<4) die("census RULE [N] [SEED]"); modeCensus(argv[3],argc>4?atoi(argv[4]):C.N,argc>5?strtoull(argv[5],0,10):1); }
    else if(mode=="vary"){ modeVary(argc>3?atoi(argv[3]):300); }
    else if(mode=="probe"){ if(argc<4) die("probe RULE [N] [SEED]"); modeProbe(argv[3],argc>4?atoi(argv[4]):C.N,argc>5?strtoull(argv[5],0,10):1); }
    else if(mode=="dump"){ if(argc<5) die("dump RULE SEED"); modeDump(argv[3],strtoull(argv[4],0,10)); }
    else if(mode=="climb"){ modeClimb(argc>3?atoi(argv[3]):1000); }
    else if(mode=="explore"){ modeExplore(argc>3?atoi(argv[3]):1000,argc>4?atoi(argv[4]):64); }
    else die("unknown mode "+mode);
    if(logFp) fclose(logFp);
}
