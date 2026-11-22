# usage: headcats.py [--window N] out_prefix barrister.out [...]
# Catalysts for the head 2bo$b3o$2o2bo (B2i3/S23-c) that stay intact >= 30 gens after the collision.
# For each Completed result: drop catalyst parts (islands 2+ cells apart) not disturbed within
# N (default 40) + 30 gens of the first contact (Barrister's filler), then re-simulate with casim: the catalyst
# must be a still life and be back intact (its cells and 1-cell ring) for >= 30 consecutive gens
# within 41 gens of the first contact. Survivors are grouped and labelled by ../pipeline/reactions.py
# (same first-hit gen and leftover at gen 300, first 8 gens after the hit within 6 cells).
# Each line gives the smallest catalyst's RLE (head + catalyst); leftovers identical to what the head
# leaves alone (a 12-cell p2 oscillator) are marked "same as the head alone".
import sys, os, math, collections, multiprocessing as mp
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim")); import casim as C
R="B2i3/S23-c"; WIN=40; STAY=30
if "--window" in sys.argv:   # the TOML's active-window maximum
    k=sys.argv.index("--window"); WIN=int(sys.argv[k+1]); del sys.argv[k:k+2]
prefix=sys.argv[1]; files=sys.argv[2:]
sys.argv=["reactions.py", R, os.path.join(HERE,"head.base.rle"), prefix]
sys.path.insert(0, os.path.join(HERE,"..","pipeline")); import reactions as RX
HEAD=set(C.parse_rle(f"x=0,y=0,rule={R}\n2bo$b3o$2o2bo!").coords())
def ring(s): return {(x+a,y+b) for x,y in s for a in(-1,0,1) for b in(-1,0,1)}
def islands(s):
    left=set(s); out=[]
    while left:
        st=[left.pop()]; c=set(st)
        while st:
            x,y=st.pop()
            for nb in [(x+i,y+j) for i in(-2,-1,0,1,2) for j in(-2,-1,0,1,2)]:
                if nb in left: left.discard(nb); c.add(nb); st.append(nb)
        out.append(frozenset(c))
    return out
def reduce_check(cat):
    """(core catalyst relative to the head, first contact, gen it is back) or (None, reason)."""
    cat=set(cat)
    q=C.Pattern.from_coords(list(cat),R); q.step(1)
    if set(q.coords())!=cat: return None,"unstable"
    parts=islands(cat); touched=set(); hit=None
    q=C.Pattern.from_coords(list(HEAD|cat),R)
    for g in range(400):
        s=set(q.coords())
        for i,c in enumerate(parts):
            if i not in touched and any((p in s)!=(p in c) for p in ring(c)): touched.add(i)
        if hit is None and touched: hit=g
        if hit is not None and g>hit+WIN+STAY: break
        q.step(1)
    if hit is None: return None,"never touched"
    core=frozenset().union(*[parts[i] for i in touched])
    Z=ring(core); q=C.Pattern.from_coords(list(HEAD|core),R); hit2=None; streak=0
    for g in range(hit+WIN+STAY+2):
        s=set(q.coords())
        ok=all((p in s)==(p in core) for p in Z)
        if hit2 is None and not ok: hit2=g
        if hit2 is not None:
            streak = streak+1 if ok else 0
            if streak>=STAY: return (core,hit2,g-STAY+1),None
        q.step(1)
    return None,"not intact for %d gens"%STAY
if __name__=="__main__":
    cats=set(); tot=0
    for fn in files:
        L=open(fn).read().split("\n")
        for i,l in enumerate(L):
            if not (l.startswith("Completed:") and L[i+1].strip()): continue
            tot+=1
            full=set(C.parse_rle(f"x=0,y=0,rule={R}\n"+L[i+1]).coords())
            a=min(HEAD)
            d=next(((fx-a[0],fy-a[1]) for fx,fy in full if all((x+fx-a[0],y+fy-a[1]) in full for x,y in HEAD)),None)
            if d is None: continue
            cats.add(frozenset((x-d[0],y-d[1]) for x,y in full)-frozenset(HEAD))
    cats=list(cats)
    with mp.Pool(10) as p: chk=p.map(reduce_check,cats,chunksize=10)
    bad=collections.Counter(); cores={}
    for c,(ok,why) in zip(cats,chk):
        if ok is None: bad[why]+=1; continue
        cores.setdefault(ok[0],ok)
    cl=list(cores)
    with mp.Pool(10) as p: res=p.map(RX.work,[set(c) for c in cl],chunksize=10)
    groups=collections.defaultdict(list)
    for c,r in zip(cl,res):
        if r[0]!="ok": bad["reactions.py: "+r[0]]+=1; continue
        _,sig,left,hit,lab=r
        groups[RX.digest(sig)].append((len(c),c,left,hit,lab,cores[c][2]))
    # merge like reactions.py (same hit gen + leftover, first 8 gens within 6 cells)
    sigs={}
    for c,r in zip(cl,res):
        if r[0]=="ok": sigs[RX.digest(r[1])]=r[1]
    keys=list(groups); parent={k:k for k in keys}
    def find(k):
        while parent[k]!=k: parent[k]=parent[parent[k]]; k=parent[k]
        return k
    bucket=collections.defaultdict(list)
    for k in keys: bucket[(sigs[k][0],sigs[k][2])].append(k)
    for b in bucket.values():
        for i in range(len(b)):
            for j in range(i+1,len(b)):
                A=sigs[b[i]][1]; B=sigs[b[j]][1]
                if len(A)==len(B) and all(len(x^y)<=RX.TOL for x,y in zip(A,B)): parent[find(b[i])]=find(b[j])
    merged=collections.defaultdict(list)
    for k in keys: merged[find(k)]+=groups[k]
    order=sorted(merged.values(), key=lambda v:(min(e[2] for e in v), min(e[0] for e in v)))
    q=C.Pattern.from_coords(list(HEAD),R); q.step(300); natural=frozenset(q.coords())
    lines=[]; ents=[]
    for i,v in enumerate(order):
        n,c,left,hit,lab,back=min(v,key=lambda e:e[0])
        q=C.Pattern.from_coords(list(HEAD|c),R); q.step(300); rest=frozenset(set(q.coords())-c)
        same=" (same as the head alone)" if rest==natural else ""
        ents.append(HEAD|c)
        lines.append(f"{i+1}: catalyst {n} cells, first contact gen {hit}, back from gen {back}, leftover at 300 {left} cells: {lab}{same}; {len(v)} variant(s)\n    "
                     + C.Pattern.from_coords(list(HEAD|c),R).to_rle().split("\n",1)[1].replace("\n",""))
    summary=(f"{R} head catalysts (intact >= {STAY} gens): {tot} completions, {len(cats)} distinct, "
             f"{len(cores)} pass after removing filler, {len(order)} reactions; rejected {dict(bad)}")
    print(summary)
    open(prefix+".txt","w").write(summary+"\n"+"\n".join(lines)+"\n")
    if ents:
        norm=[]
        for e in ents:
            mx=min(x for x,y in e); my=min(y for x,y in e); norm.append({(x-mx,y-my) for x,y in e})
        W=max(max(x for x,y in e) for e in norm)+12; H=max(max(y for x,y in e) for e in norm)+12
        cols=max(1,math.ceil(math.sqrt(len(norm)*H/W))); out=set()
        for j,e in enumerate(norm): out|={(x+(j%cols)*W,y+(j//cols)*H) for x,y in e}
        open(prefix+".rle","w").write(f"#C {summary}\n#C entries in .txt order, {cols} per row, spacing {W}x{H}\n"
                                       +C.Pattern.from_coords(list(out),R).to_rle())
    for l in lines[:10]: print(" ",l)
