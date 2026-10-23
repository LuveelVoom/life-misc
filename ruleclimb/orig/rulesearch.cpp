// rulesearch: constrained search over range-1 Moore rules (totalistic / isotropic / MAP).
// Build: g++ -O2 -march=native -fopenmp -std=c++17 -o rulesearch rulesearch.cpp
//
// Conventions (checked against Golly gollybase/liferules.cpp):
//   512-entry table index: row-major from NW, NW = bit 8 (256) ... SE = bit 0, centre = bit 4 (16).
//   MAP: 86 base64 chars, 6 bits per char MSB first, bit k of the stream = table[k].
//   Hensel letters/neighbourhoods: Golly's entry tables; counts 5..8 are complements of 3..0.
//   Patterns: x to the right, y down (RLE order).
#include <bits/stdc++.h>
#include <omp.h>
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

static Cells stepCells(const Cells&c,const uint8_t*T,vector<pair<int,int>>*used=nullptr){
    Cells out; if(c.empty()) return out;
    int mx,my,Mx,My; bbox(c,mx,my,Mx,My);
    int W=Mx-mx+5,H=My-my+5; vector<uint8_t> g(W*H,0);
    for(auto&p:c) g[(p.second-my+2)*W+(p.first-mx+2)]=1;
    for(int y=1;y<H-1;y++) for(int x=1;x<W-1;x++){
        int i=(g[(y-1)*W+x-1]<<8)|(g[(y-1)*W+x]<<7)|(g[(y-1)*W+x+1]<<6)|(g[y*W+x-1]<<5)|(g[y*W+x]<<4)|(g[y*W+x+1]<<3)|(g[(y+1)*W+x-1]<<2)|(g[(y+1)*W+x]<<1)|g[(y+1)*W+x+1];
        if(used) used->push_back({i,T[i]});
        if(T[i]) out.push_back({x+mx-2,y+my-2});
    }
    return out;
}
// evolve until normalized state repeats (translation allowed) or dies. returns gens to repeat t and earlier index s (s=-1: died at t)
struct Evo { vector<Cells> gens; int s=-2,t=0; int dx=0,dy=0; };
static Evo evolve(const Cells&c0,const uint8_t*T,int maxg,vector<pair<int,int>>*used){
    Evo E; map<string,int> seen; Cells c=c0;
    for(int g=0;g<=maxg;g++){
        E.gens.push_back(c);
        if(c.empty()&&g>0){ E.s=-1;E.t=g; return E; }
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
    string logFile="discoveries.txt", archiveFile="archive.txt"; vector<string> starts; int maxgCons=4000;
} C;
static vector<int> forcedV; // per class: -1 free, 0/1 forced

static void loadCfg(const string&fn){
    ifstream f(fn); if(!f) die("cannot open config "+fn); string line; int ln=0;
    vector<tuple<string,string,bool>> pats; vector<tuple<string,double>> tgs; map<string,double> minr;
    while(getline(f,line)){
        ln++; auto h=line.find('#'); if(h!=string::npos) line=line.substr(0,h);
        istringstream is(line); string k; if(!(is>>k)) continue;
        auto need=[&](bool ok){ if(!ok) die("config line "+to_string(ln)+": "+line); };
        if(k=="rule"){ need((bool)(is>>C.refRule)); }
        else if(k=="symmetry"){ string s; is>>s; SYM= s=="totalistic"?0: s=="isotropic"?1: s=="none"?2: -1; need(SYM>=0); }
        else if(k=="pattern"){ string n,r,opt; need((bool)(is>>n>>r)); is>>opt; pats.push_back({n,r,opt=="split"}); }
        else if(k=="target"){ string n; double w=1; need((bool)(is>>n)); is>>w; tgs.push_back({n,w}); }
        else if(k=="min_rate"){ string n; double r; need((bool)(is>>n>>r)); minr[n]=r; }
        else if(k=="fix"){ string n; int v; need((bool)(is>>n>>v)); C.fixes.push_back({n,v}); }
        else if(k=="soup"){ need((bool)(is>>C.soupW>>C.soupH>>C.dens)); }
        else if(k=="soups") is>>C.N; else if(k=="grid") is>>C.grid; else if(k=="gmax") is>>C.gmax; else if(k=="pmax") is>>C.pmax;
        else if(k=="maxp") is>>C.maxp; else if(k=="cap") is>>C.cap; else if(k=="max_bad") is>>C.maxBad;
        else if(k=="cost_weight") is>>C.costW; else if(k=="cost_ref") is>>C.costRef; else if(k=="max_cost") is>>C.maxCost;
        else if(k=="min_osc_period") is>>C.minOscPeriod; else if(k=="lahc") is>>C.lahc; else if(k=="reseed") is>>C.reseed;
        else if(k=="seed") is>>C.seed; else if(k=="log") is>>C.logFile; else if(k=="archive") is>>C.archiveFile;
        else if(k=="start"){ string r; is>>r; C.starts.push_back(r); }
        else if(k=="constraint_gens") is>>C.maxgCons;
        else if(k=="soup_cost_limit") is>>C.soupCostLimit;
        else die("unknown config key '"+k+"' at line "+to_string(ln));
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
        vector<pair<int,int>> used; P.evo=evolve(P.cells,C.ref.t,C.maxgCons,&used);
        if(P.evo.s==-2) die("pattern "+P.name+" does not repeat within "+to_string(C.maxgCons)+" gens; raise constraint_gens or split it");
        for(auto&u:used) force(u.first,u.second,P.name);
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
struct Cand { long count=0; Cells ex; Info I; u64 seed=0; };
struct SoupRes { uint8_t cnt[8]; uint8_t bad; double cost; };
struct Engine {
    const uint8_t*T=nullptr; uint8_t T2[512]; int W,H; vector<uint8_t> a,b; vector<u64> Z; vector<int> stamp; int st=0;
    unordered_map<string,Info> cache, cache2; unordered_map<string,Cand> cands; u64 curSeed=0; bool countUsage=false; vector<u64> usage; // usage per 512 index
    void setRule(const uint8_t*t){ T=t; cache.clear(); cache2.clear(); cands.clear();
        for(int j=0;j<512;j++){ int m=0; for(int c=0;c<3;c++) for(int r=0;r<3;r++){ int col=(j>>(6-3*c))&7; if((col>>(2-r))&1) m|=1<<(8-(r*3+c)); } T2[j]=t[m]; } }
    Engine(int g):W(g),H(g),a(g*g),b(g*g),Z(g*g),stamp(g*g,0),usage(512,0){ u64 s=12345; for(auto&z:Z) z=splitmix(s); }
    Info& cls1(const Cells&c){ string k=keyOf(c); auto it=cache.find(k); if(it!=cache.end()) return it->second; if(cache.size()>200000) cache.clear(); return cache[k]=classify(c,T,64); }
    Info& cls2(const Cells&c){ string k=keyOf(c); auto it=cache2.find(k); if(it!=cache2.end()) return it->second; if(cache2.size()>50000) cache2.clear(); return cache2[k]=classify(c,T,C.maxp); }
    void note(const Cells&c,const Info&I){ string k=keyOf(c); auto&e=cands[k]; if(!e.count){e.ex=c;e.I=I;e.seed=curSeed;} e.count++; }
    void run(u64 seed, SoupRes&R){
        memset(&R,0,sizeof R); curSeed=seed;
        fill(a.begin(),a.end(),0); fill(b.begin(),b.end(),0);
        u64 s=seed; u64 h=0; int x0=W/2-C.soupW/2,y0=H/2-C.soupH/2;
        auto set1=[&](int x,int y){ a[(y0+y)*W+x0+x]=1; h^=Z[(y0+y)*W+x0+x]; };
        if(C.soupW<=64 && C.dens==0.5){ for(int y=0;y<C.soupH;y++){u64 r=splitmix(s); for(int x=0;x<C.soupW;x++) if(r>>x&1) set1(x,y);} }
        else { u64 thr=(u64)(C.dens*18446744073709551615.0); for(int y=0;y<C.soupH;y++) for(int x=0;x<C.soupW;x++) if(splitmix(s)<thr) set1(x,y); }
        int bx0=x0,bx1=x0+C.soupW-1,by0=y0,by1=y0+C.soupH-1; int pbx0=bx0,pbx1=bx1,pby0=by0,pby1=by1;
        u64 hist[256]; int t=0; bool stable=false; long lastpop=1; double cost=0;
        for(t=0;t<C.gmax;t++){
            hist[t&255]=h;
            if(t>=20){ for(int p=1;p<=60&&2*p<=t;p++) if(hist[(t-p)&255]==h&&hist[(t-2*p)&255]==h){stable=true;break;} if(stable) break; }
            if(t%8==0&&t>0&&lastpop<800) scan(R,false,bx0,bx1,by0,by1,h);
            int sx0=max(1,min(bx0,pbx0)-1), sx1=min(W-2,max(bx1,pbx1)+1), sy0=max(1,min(by0,pby0)-1), sy1=min(H-2,max(by1,pby1)+1);
            pbx0=bx0;pbx1=bx1;pby0=by0;pby1=by1;
            int nx0=INT_MAX,nx1=-1,ny0=INT_MAX,ny1=-1; long pop=0;
            cost+=(double)(sx1-sx0+1)*(sy1-sy0+1);
            if(cost>C.soupCostLimit*1e6){ R.bad=3; R.cost=cost; return; }
            for(int y=sy0;y<=sy1;y++){
                const uint8_t*r0=&a[(y-1)*W],*r1=&a[y*W],*r2=&a[(y+1)*W]; uint8_t*o=&b[y*W];
                int x=sx0-1; int j=(r0[x]<<2)|(r1[x]<<1)|r2[x]; x++; j=(j<<3)|(r0[x]<<2)|(r1[x]<<1)|r2[x];
                if(countUsage){ for(x=sx0;x<=sx1;x++){ j=((j<<3)&0x1ff)|(r0[x+1]<<2)|(r1[x+1]<<1)|r2[x+1]; usage[j]++;
                        uint8_t v=T2[j]; o[x]=v; if(v!=r1[x]) h^=Z[y*W+x]; if(v){pop++; nx0=min(nx0,x); nx1=max(nx1,x); ny0=min(ny0,y); ny1=max(ny1,y);} } }
                else { for(x=sx0;x<=sx1;x++){ j=((j<<3)&0x1ff)|(r0[x+1]<<2)|(r1[x+1]<<1)|r2[x+1];
                        uint8_t v=T2[j]; o[x]=v; if(v!=r1[x]) h^=Z[y*W+x]; if(v){pop++; if(x<nx0)nx0=x; if(x>nx1)nx1=x; if(y<ny0)ny0=y; if(y>ny1)ny1=y;} } }
            }
            swap(a,b);
            if(pop==0){ R.cost=cost; return; }
            bx0=nx0;bx1=nx1;by0=ny0;by1=ny1; lastpop=pop;
            if(pop>C.pmax){ R.bad=1; R.cost=cost; return; }
        }
        R.cost=cost;
        if(!stable) R.bad=2;
        scan(R,true,bx0,bx1,by0,by1,h);
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
            Info*Ip=&cls1(c); if(final&&Ip->type==4&&C.maxp>64) Ip=&cls2(c); Info&I=*Ip;
            string k;
            if(I.type==2){
                int dx=(I.dx>0)-(I.dx<0), dy=(I.dy>0)-(I.dy<0);
                int back=INT_MAX; for(auto&p:c) back=min(back,p.first*dx+p.second*dy);
                int cx=(mx+Mx)/2, cy=(my+My)/2; bool threat=false;
                for(int p:live){ if(stamp[p]==(int)(sbase+id)) continue; int x=p%W,y=p/W; int along=x*dx+y*dy-back; if(along<-2) continue;
                    int lat=abs((x-cx)*dy-(y-cy)*dx)-(Mx-mx+My-my); if(lat<=along+6){threat=true;break;} }
                if(threat) continue;
                k=keyOf(c); bool tgt=false;
                for(size_t t=0;t<C.tg.size();t++) if(C.tg[t].keys.count(k)){ if(R.cnt[t]<255) R.cnt[t]++; tgt=true; }
                if(!tgt) note(c,I);
                for(int p:cl){a[p]=0; h^=Z[p];}
            } else if(final && I.type==1){
                bool ok=true;
                for(auto&e:I.env){int ex=e.first+mx,ey=e.second+my;
                    for(int yy=ey-2;yy<=ey+2&&ok;yy++) for(int xx=ex-2;xx<=ex+2;xx++){ if(xx<0||yy<0||xx>=W||yy>=H) continue; int n=yy*W+xx; if(a[n]&&stamp[n]!=(int)(sbase+id)){ok=false;break;} } if(!ok)break; }
                if(!ok) continue;
                k=keyOf(c); bool tgt=false;
                for(size_t t=0;t<C.tg.size();t++) if(C.tg[t].keys.count(k)){ if(R.cnt[t]<255) R.cnt[t]++; tgt=true; }
                if(!tgt && I.period>=C.minOscPeriod) note(c,I);
            }
        }
    }
};
static vector<Engine*> engines;

// ===================== evaluation =====================
struct Eval {
    Table T; int n=0; vector<SoupRes> res; vector<double> usage; // per class
    unordered_map<string,Cand> cands; double f=-1e9; vector<double> rate; double bad=0, cost=0; u64 seedEpoch=0;
};
static double tSim=0,tProc=0,tMisc=0; static double nowS(){ return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }
static u64 soupSeed(u64 epoch,int i){ u64 s=epoch*1000003ULL+i; return splitmix(s); }
static void evalRange(Eval&E,int i0,int i1){
    double t0=nowS(); if((int)E.res.size()<i1) E.res.resize(i1);
    #pragma omp parallel
    { Engine&G=*engines[omp_get_thread_num()];
      #pragma omp for schedule(dynamic,4)
      for(int i=i0;i<i1;i++){ G.countUsage = i<64; G.run(soupSeed(E.seedEpoch,i),E.res[i]); } }
    E.n=i1; tSim+=nowS()-t0;
}
static void beginEval(Eval&E,const Table&T,u64 epoch){
    E.T=T; E.n=0; E.res.clear(); E.cands.clear(); E.seedEpoch=epoch;
    for(auto*g:engines){ g->setRule(E.T.t); fill(g->usage.begin(),g->usage.end(),0); }
}
static void finishEval(Eval&E){
    // merge candidates and usage
    E.usage.assign(NC,0);
    for(auto*g:engines){
        for(int j=0;j<512;j++){ int m=0; for(int c=0;c<3;c++) for(int r=0;r<3;r++){ int col=(j>>(6-3*c))&7; if((col>>(2-r))&1) m|=1<<(8-(r*3+c)); } E.usage[cls[m]]+=g->usage[j]; }
        for(auto&kv:g->cands){ auto&e=E.cands[kv.first]; if(!e.count){e.ex=kv.second.ex;e.I=kv.second.I;e.seed=kv.second.seed;} e.count+=kv.second.count; }
        g->cands.clear();
    }
}
static void score(Eval&E){
    int n=E.n, nt=C.tg.size(); E.rate.assign(nt,0); double bad=0,cost=0;
    for(int i=0;i<n;i++){ auto&r=E.res[i]; cost+=r.cost; if(r.bad){bad++;continue;} for(int t=0;t<nt;t++) E.rate[t]+=min((int)r.cnt[t],C.cap); }
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
    int ax=abs(I.dx),ay=abs(I.dy); if(ax<ay) swap(ax,ay); int g=__gcd(ax,__gcd(ay,I.period));
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
    for(auto&kv:E.cands){
        const Cand&c=kv.second; if(isPseudo(c.ex,E.T.t,c.I.period)) continue;
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
static double interest(const Find&F,bool newDesc){ if(!newDesc) return 0.1; return F.I.type==2 ? 5.0 : 1.0+log2((double)F.I.period)/2; }
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
            d.push_back(g(M.res[i],false)-g(P.res[i],true));
        }
        i0=i1; if(i1>=C.N) break;
        { double cs=0,bd=0; for(int i=0;i<i1;i++){ cs+=M.res[i].cost; bd+=M.res[i].bad?1:0; } cs/=i1*1e6; bd/=i1;
          if(C.maxCost>0 && cs>1.5*C.maxCost && cs>P.cost) return false;
          if(bd>C.maxBad+0.1 && bd>P.bad+0.05) return false; }
        if(d.size()>=(size_t)batch*2){ double m=0,v=0; for(double x:d) m+=x; m/=d.size(); for(double x:d) v+=(x-m)*(x-m); v/=(d.size()-1);
            if(m+2.5*sqrt(v/d.size())<0) return false; }
    }
    return true;
}

static void printEval(const char*tag,const Eval&E){
    printf("%s f=%.4f",tag,E.f); for(size_t t=0;t<C.tg.size();t++) printf(" %s=%.3f",C.tg[t].name.c_str(),E.rate[t]);
    printf(" bad=%.3f cost=%.3fM %s %s\n",E.bad,E.cost,feasible(E)?"feasible":"INFEASIBLE",ruleString(E.T).c_str()); fflush(stdout);
}

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
    auto fs=processCands(E);
    for(auto&F:fs){ printf("%6ld  %-34s %3zu cells  %s",F.count,describe(F.I).c_str(),F.phase.size(),archiveIds.count(F.id)?"":"(not in archive) "); string r=toRLE(F.phase,ruleString(T)); printf("\n%s",r.c_str()); }
}
static void modeDump(const string&rule,u64 seed){
    Table T=parseRule(rule); Cells c; u64 s=seed;
    if(C.soupW<=64&&C.dens==0.5){ for(int y=0;y<C.soupH;y++){u64 r=splitmix(s); for(int x=0;x<C.soupW;x++) if(r>>x&1) c.push_back({x,y});} }
    else { u64 thr=(u64)(C.dens*18446744073709551615.0); for(int y=0;y<C.soupH;y++) for(int x=0;x<C.soupW;x++) if(splitmix(s)<thr) c.push_back({x,y}); }
    printf("%s",toRLE(c,ruleString(T)).c_str());
}
static Table startRule(int i){ return project(parseRule(C.starts.empty()? C.refRule : C.starts[i%C.starts.size()])); }

static void modeClimb(int iters){
    u64 rng=C.seed*7919+17, epoch=C.seed*1000+1; Eval cur; fullEval(cur,startRule(0),epoch);
    if(C.costW>0&&C.costRef<=0){ C.costRef=cur.cost; score(cur); printf("cost_ref set to start cost %.3fM\n",C.costRef); }
    printEval("start",cur); recordFinds(cur,processCands(cur),"start");
    Eval best=cur; vector<double> hist(max(1,C.lahc),cur.f); long evals=0, early=0;
    for(int it=1;it<=iters;it++){
        if(it%50==0){ printf("[it %d: %ld evals, %ld rejected early, current f %.4f, best f %.4f; time sim %.0fs, objects %.0fs]\n",it,evals,early,cur.f,best.f,tSim,tProc); fflush(stdout); }
        if(C.reseed>0&&it%C.reseed==0){ epoch++; Table T=cur.T; fullEval(cur,T,epoch); fill(hist.begin(),hist.end(),cur.f); printEval(("reseed@"+to_string(it)).c_str(),cur); }
        string desc; Table MT=mutate(cur.T,cur.usage,rng,desc);
        Eval M; beginEval(M,MT,epoch); evals++;
        if(!race(M,cur,max(32,C.N/8))){ early++; continue; }
        finishEval(M); score(M);
        double h=hist[it%hist.size()];
        if(M.f>=cur.f||M.f>=h){
            bool improved=M.f>cur.f; cur=std::move(M); hist[it%hist.size()]=cur.f;
            if(improved){ printEval(("it "+to_string(it)+" ["+desc+"]").c_str(),cur); }
            if(feasible(cur)) recordFinds(cur,processCands(cur),"climb it "+to_string(it));
            if(cur.f>best.f){ best=cur; }
        } else hist[it%hist.size()]=cur.f;
    }
    printEval("best",best);
}

struct Member { Eval E; double credit=0; int visits=0; vector<string> finds; };
static void modeExplore(int iters,int cap){
    u64 rng=C.seed*104729+3, epoch=C.seed*1000+7;
    vector<Member> pop;
    for(size_t i=0;i<max<size_t>(1,C.starts.size());i++){ Member m; fullEval(m.E,startRule(i),epoch);
        if(C.costW>0&&C.costRef<=0){ C.costRef=m.E.cost; score(m.E); printf("cost_ref set to start cost %.3fM\n",C.costRef); }
        printEval("start",m.E); auto fs=processCands(m.E); m.credit=recordFinds(m.E,fs,"start"); pop.push_back(std::move(m)); }
    long evals=0,early=0,infeas=0,added=0;
    for(int it=1;it<=iters;it++){
        // parent: weight (1+credit)/sqrt(1+visits) * feasibility
        vector<double> w(pop.size()); double tot=0; for(size_t i=0;i<pop.size();i++){ w[i]=(1+pop[i].credit)/sqrt(1.0+pop[i].visits)*(feasible(pop[i].E)?1:0.2); tot+=w[i]; }
        double u=(splitmix(rng)>>11)*(1.0/9007199254740992.0)*tot; size_t pi=0; while(pi+1<pop.size()&&u>w[pi]){u-=w[pi];pi++;}
        Member&P=pop[pi]; P.visits++;
        string desc; Table MT=mutate(P.E.T,P.E.usage,rng,desc);
        Member m; beginEval(m.E,MT,epoch); evals++;
        // quick feasibility check on first batch
        int b=max(32,C.N/8); evalRange(m.E,0,32); score(m.E);
        if(m.E.bad>0.25 || (C.maxCost>0&&m.E.cost>3*C.maxCost)){ early++; continue; }
        evalRange(m.E,32,b); score(m.E);
        bool hopeless=m.E.bad>2*C.maxBad+0.05; for(size_t t=0;t<C.tg.size();t++) hopeless|= m.E.rate[t]<0.6*C.tg[t].minRate; if(C.maxCost>0) hopeless|=m.E.cost>1.5*C.maxCost;
        if(hopeless){ early++; continue; }
        evalRange(m.E,b,C.N); finishEval(m.E); score(m.E);
        if(!feasible(m.E)){ infeas++; continue; }
        auto fs=processCands(m.E);
        m.credit=recordFinds(m.E,fs,"explore it "+to_string(it)+", parent ["+desc+"]");
        for(auto&F:fs) m.finds.push_back(F.id);
        if(m.credit>0 || m.E.f>=P.E.f-0.1){
            P.credit*=0.9; // parent shares a little decay so the search keeps moving
            pop.push_back(std::move(m)); added++;
            if((int)pop.size()>cap){ size_t worst=0; double ws=1e18; for(size_t i=0;i<pop.size()-1;i++){ double s=pop[i].E.f+pop[i].credit-0.05*pop[i].visits; if(s<ws){ws=s;worst=i;} } pop.erase(pop.begin()+worst); }
        }
        if(it%25==0){ double bf=-1e9; for(auto&q:pop) bf=max(bf,q.E.f); printf("[it %d: pop %zu, %ld evals, %ld early-rejected, %ld infeasible, %ld added, archive %zu, best f %.3f]\n",it,pop.size(),evals,early,infeas,added,archiveIds.size(),bf); fflush(stdout); }
    }
    printf("final population:\n"); sort(pop.begin(),pop.end(),[](const Member&a,const Member&b){return a.E.f>b.E.f;});
    for(auto&q:pop) printEval(("  credit "+to_string(q.credit).substr(0,5)).c_str(),q.E);
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
    initHensel(); selfTest();
    if(argc<3){ fprintf(stderr,"usage: rulesearch CONFIG check [RULE...] | census RULE [N] [SEED] | climb ITERS | explore ITERS [POPCAP] | dump RULE SEED\n"); return 1; }
    loadCfg(argv[1]); string mode=argv[2];
    int nt=omp_get_max_threads(); for(int i=0;i<nt;i++) engines.push_back(new Engine(C.grid));
    loadArchive(); if(!C.logFile.empty()) logFp=fopen(C.logFile.c_str(),"a");
    if(mode=="check"){ vector<string> r; for(int i=3;i<argc;i++) r.push_back(argv[i]); modeCheck(r); }
    else if(mode=="census"){ if(argc<4) die("census RULE [N] [SEED]"); modeCensus(argv[3],argc>4?atoi(argv[4]):C.N,argc>5?strtoull(argv[5],0,10):1); }
    else if(mode=="dump"){ if(argc<5) die("dump RULE SEED"); modeDump(argv[3],strtoull(argv[4],0,10)); }
    else if(mode=="climb"){ modeClimb(argc>3?atoi(argv[3]):1000); }
    else if(mode=="explore"){ modeExplore(argc>3?atoi(argv[3]):1000,argc>4?atoi(argv[4]):64); }
    else die("unknown mode "+mode);
    if(logFp) fclose(logFp);
}
