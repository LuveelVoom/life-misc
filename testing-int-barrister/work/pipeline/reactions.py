# usage: reactions.py [--exact] [--min-leftover N] rule base.rle out_prefix barrister.out [...]
# --min-leftover N: list only reactions that leave at least N cells besides the catalyst at gen 300
# (N=1 drops plain eaters). Each reaction is labelled with what its leftover is: spaceships with
# direction (y down), still lifes, oscillators, other.
# Default grouping ("same main reaction"): same generation of the catalyst's first hit, identical
# leftover at gen 300, and the first 8 generations after the hit (pattern XOR catalyst) differing by
# at most 6 cells per generation (merged transitively). Variants that differ only in secondary still
# lifes, or in small details of the catalyst, merge; a different outcome always separates them.
# --exact uses the full signature below (every generation until recovery).
# Groups results by the actual reaction. Signature: for every generation from the catalyst's first hit
# until it has recovered (10 quiet gens), the pattern XOR the catalyst (= active cells plus catalyst
# cells that changed), plus the leftover at gen 300. Catalyst parts that are never touched cancel, so
# decorated variants of one mechanism collapse. Signatures are taken up to the symmetries of the base
# object (e.g. mirror images for the pi). Keeps catalysts stable in the rule that recover for >= 10
# gens (later damage is allowed and shows in the leftover). Writes out_prefix.rle (smallest catalyst per reaction; compact grid, split per 400) and .txt.
import sys, os, math, hashlib, collections, multiprocessing as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from align import offset
if "--exact" in sys.argv:
    sys.argv.remove("--exact"); os.environ["REACTIONS_EXACT"]="1"   # worker processes re-read argv, so use env
EXACT = os.environ.get("REACTIONS_EXACT")=="1"
MINLEFT = 0
if "--min-leftover" in sys.argv:
    k=sys.argv.index("--min-leftover"); MINLEFT=int(sys.argv[k+1]); del sys.argv[k:k+2]
R=sys.argv[1]; T=300
def S(q): return set(map(tuple,q.coords() if callable(q.coords) else q.coords))
def zoi(s): return {(x+a,y+b) for x,y in s for a in (-1,0,1) for b in (-1,0,1)}
base=S(C.parse_rle(open(sys.argv[2]).read()))
def orient(cells,t):
    out=set()
    for x,y in cells:
        if t&4: x,y=y,x
        if t&1: x=-x
        if t&2: y=-y
        out.add((x,y))
    return out
# symmetries of the base: transforms t with t(base) == base + shift
SYMS=[]
for t in range(8):
    tb=orient(base,t); dx=min(x for x,y in base)-min(x for x,y in tb); dy=min(y for x,y in base)-min(y for x,y in tb)
    if {(x+dx,y+dy) for x,y in tb}==base: SYMS.append((t,dx,dy))
def signature(cat):
    q=C.Pattern.from_coords(list(base|cat),R); Z=zoi(cat); seq=[]; hit=None; streak=0; recovered=False
    for g in range(T):
        s=S(q)
        if not recovered:
            if any(((c in s)!=(c in cat)) for c in Z):
                if hit is None: hit=g
                streak=0
            elif hit is not None:
                streak+=1
            if hit is not None: seq.append(frozenset(s ^ cat))
            if streak>=10: recovered=True   # later damage only shows in the leftover
        q.step(1)
    if not recovered: return None
    return (hit, tuple(seq), frozenset(S(q) ^ cat))
def islands(cat):
    left=set(cat); out=[]
    while left:
        st=[left.pop()]; comp=set(st)
        while st:
            x,y=st.pop()
            for c in [c for c in left if max(abs(c[0]-x),abs(c[1]-y))<=2]:
                left.discard(c); comp.add(c); st.append(c)
        out.append(frozenset(comp))
    return out
W_MAIN=8
TOL=6   # cells per generation
def coarse(cat):
    # main reaction: the first W_MAIN gens after the first hit (pattern XOR catalyst), plus the leftover
    q=C.Pattern.from_coords(list(base|cat),R); Z=zoi(cat); seq=[]; hit=None; streak=0; recovered=False
    for g in range(T):
        s=S(q)
        if not recovered:
            if any(((c in s)!=(c in cat)) for c in Z):
                if hit is None: hit=g
                streak=0
            elif hit is not None:
                streak+=1
                if streak>=10: recovered=True
        if hit is not None and g < hit+W_MAIN: seq.append(frozenset(s ^ cat))
        q.step(1)
    if not recovered: return None
    return (hit, tuple(seq), frozenset(S(q) ^ cat))
def digest(sig):
    hit,seq,left=sig
    h=hashlib.sha1(str(hit).encode())
    for d in tuple(seq)+(left,): h.update(repr(sorted(d)).encode()); h.update(b"|")
    return h.hexdigest()
def work(cat):
    q=C.Pattern.from_coords(list(cat),R); q.step(1)
    if S(q)!=set(cat): return ("unstable",None,None,None)
    best=None; first=None
    for t,dx,dy in SYMS:
        c={(x+dx,y+dy) for x,y in orient(cat,t)}
        sig=signature(c) if EXACT else coarse(c)
        if sig is None: return ("norecover",None,None,None)
        if first is None: first=sig
        if best is None or digest(sig)<digest(best): best=sig
    return ("ok",best,len(first[2]),first[0],label(first[2],set(cat)))
DIRS={(0,-1):"N",(0,1):"S",(1,0):"E",(-1,0):"W",(1,-1):"NE",(-1,-1):"NW",(1,1):"SE",(-1,1):"SW"}
def label(left,cat):
    """What the leftover at gen 300 is (catalyst cells missing at gen 300 are reported separately)."""
    alive=set(left)-cat; missing=set(left)&cat
    if not alive and not missing: return "eater (nothing left)"
    parts=collections.Counter()
    for comp in islands(alive):
        try: pr=C.Pattern.from_coords(list(comp),R).period(400)
        except Exception: pr=None
        if not pr: parts[f"{len(comp)}-cell active/growing"]+=1; continue
        k=pr["kind"]
        if k=="spaceship":
            d=(int(math.copysign(1,pr["dx"])) if pr["dx"] else 0, int(math.copysign(1,pr["dy"])) if pr["dy"] else 0)
            parts[f"ship {len(comp)} cells p{pr['period']} ({pr['dx']},{pr['dy']}) going {DIRS[d]}"]+=1
        elif k=="still life": parts[f"{len(comp)}-cell still life"]+=1
        elif k=="oscillator": parts[f"p{pr['period']} oscillator ({len(comp)} cells)"]+=1
        elif k=="dies": parts["dying debris"]+=1
        else: parts[k]+=1
    out=", ".join((f"{n}x " if n>1 else "")+k for k,n in sorted(parts.items(), key=lambda kv:("ship" not in kv[0], kv[0])))
    if missing: out+=(", " if out else "")+f"catalyst lost {len(missing)} cells later"
    return out
if __name__=="__main__":
    cats=set(); tot=0
    for f in sys.argv[4:]:
        L=open(f).read().split("\n")
        for i,l in enumerate(L):
            if not (l.startswith("Completed:") and L[i+1].strip()): continue
            tot+=1
            full=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+L[i+1]))
            off=offset(L, i, base)
            if off is None: continue
            dx,dy=off
            cats.add(frozenset((x-dx,y-dy) for x,y in full)-frozenset(base))
    cats=list(cats)
    with mp.Pool(10) as p: res=p.map(work,[set(c) for c in cats],chunksize=20)
    groups=collections.defaultdict(list); bad=collections.Counter(); sigof={}
    for c,r in zip(cats,res):
        st=r[0]
        if st!="ok": bad[st]+=1; continue
        _,sig,left,hit,lab=r
        if left<MINLEFT: bad["leftover<%d"%MINLEFT]+=1; continue
        key=digest(sig); sigof[key]=sig
        groups[key].append((len(c),c,left,hit,lab))
    if not EXACT:
        # merge main reactions with the same first-hit gen and leftover whose first 8 gens after the
        # hit never differ by more than TOL cells (transitively)
        keys=list(groups); parent={k:k for k in keys}
        def find(k):
            while parent[k]!=k: parent[k]=parent[parent[k]]; k=parent[k]
            return k
        bucket=collections.defaultdict(list)
        for k in keys: bucket[(sigof[k][0], sigof[k][2])].append(k)
        for b in bucket.values():
            for i in range(len(b)):
                for j in range(i+1,len(b)):
                    A=sigof[b[i]][1]; B=sigof[b[j]][1]
                    if len(A)==len(B) and all(len(x^y)<=TOL for x,y in zip(A,B)): parent[find(b[i])]=find(b[j])
        merged=collections.defaultdict(list)
        for k in keys: merged[find(k)]+=groups[k]
        groups=merged
    order=sorted(groups.values(), key=lambda v:(min(e[2] for e in v), min(e[0] for e in v)))
    ents=[]; lines=[]
    for i,v in enumerate(order):
        n,c,left,hit,lab=min(v,key=lambda e:e[0])
        ents.append(base|c); lines.append(f"{i+1}: catalyst {n} cells, first hit gen {hit}, leftover {left} cells: {lab}; {len(v)} variant(s)")
    summary=f"{R}: {tot} completions, {len(cats)} distinct catalysts, {sum(len(v) for v in order)} pass, {len(order)} {'exact' if EXACT else 'main'} reactions; rejected {dict(bad)}; symmetries used {len(SYMS)}"
    print(summary)
    open(sys.argv[3]+".txt","w").write(summary+"\n"+"\n".join(lines)+"\n")
    norm=[]
    for e in ents:
        mx=min(x for x,y in e); my=min(y for x,y in e); norm.append({(x-mx,y-my) for x,y in e})
    if norm:
        W=max(max(x for x,y in e) for e in norm)+15; H=max(max(y for x,y in e) for e in norm)+15
        PER=400; parts=[norm[k:k+PER] for k in range(0,len(norm),PER)]
        for pi,part in enumerate(parts):
            cols=max(1,math.ceil(math.sqrt(len(part)*H/W))); out=set()
            for j,e in enumerate(part): out|={(x+(j%cols)*W,y+(j//cols)*H) for x,y in e}
            name=sys.argv[3]+(".rle" if len(parts)==1 else f"_part{pi+1}.rle")
            hdr=f"#C {summary}\n#C one entry per {'exact reaction' if EXACT else 'main reaction (same hit gen and leftover, first 8 gens within 6 cells)'} (smallest catalyst), sorted by leftover then catalyst size; {cols} per row, spacing {W}x{H}, row-major; entries {pi*PER+1}-{pi*PER+len(part)} (details in {os.path.basename(sys.argv[3])}.txt)\n"
            open(name,"w").write(hdr+C.Pattern.from_coords(list(out),R).to_rle())
    for l in lines[:8]: print(" ",l)
