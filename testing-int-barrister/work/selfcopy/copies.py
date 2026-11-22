# usage: copies.py [--window N] [--input RLE] barrister.out [...]
# --input: the pattern the search started from (default: the target). --target: the pattern whose
# copies are looked for (default 2bo$b3o$2o2bo2$4b2o$4b2o$4b2o).
# Strict: the copy must appear within N gens (default 40, the TOML's active window) of the
# catalyst's first contact; catalyst parts (islands 2+ cells apart) not disturbed before the copy are
# removed, the copy must still form without them, and that catalyst must be back intact for 10 gens
# at some point within N + 14 gens of the first contact (it may still be recovering at the copy).
# (copies_loose.py is the old version, which also counted copies made later by debris.)
# For every Completed result: re-simulate with casim (B2i3/S23-c) for 400 gens and look, after the
# catalyst's first contact, for an isolated copy of the input pattern (any of 8 orientations, every
# cell of its 1-cell ring dead). Prints the catalysts that give one, with gen, position, orientation,
# whether the catalyst is intact then and at gen 400.
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B2i3/S23-c"; BODY="2bo$b3o$2o2bo2$4b2o$4b2o$4b2o!"
if "--target" in sys.argv:
    k=sys.argv.index("--target"); BODY=sys.argv[k+1]; del sys.argv[k:k+2]
P0=set(C.parse_rle(f"x=0,y=0,rule={R}\n{BODY}").coords())   # target
IN=BODY
if "--input" in sys.argv:
    k=sys.argv.index("--input"); IN=sys.argv[k+1]; del sys.argv[k:k+2]
I0=set(C.parse_rle(f"x=0,y=0,rule={R}\n{IN}").coords())     # input
def norm(s):
    mx=min(x for x,y in s); my=min(y for x,y in s); return frozenset((x-mx,y-my) for x,y in s)
ORIENTS=[]
for t in range(8):
    o=set()
    for x,y in P0:
        if t&4: x,y=y,x
        if t&1: x=-x
        if t&2: y=-y
        o.add((x,y))
    o=norm(o)
    if o not in [q for q,_ in ORIENTS]: ORIENTS.append((o,t))
def ring(s): return {(x+a,y+b) for x,y in s for a in(-1,0,1) for b in(-1,0,1)}-set(s)
RINGS=[(o,t,ring(o),min(o)) for o,t in ORIENTS]
def find_copies(s):
    out=[]
    for o,t,rg,a in RINGS:
        for (x,y) in s:
            dx,dy=x-a[0],y-a[1]
            if all((u+dx,v+dy) in s for u,v in o) and not any((u+dx,v+dy) in s for u,v in rg):
                out.append((dx,dy,t))
    return out
WIN=40
if "--window" in sys.argv:
    k=sys.argv.index("--window"); WIN=int(sys.argv[k+1]); del sys.argv[k:k+2]
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
def first_copy(pat, cat, gens):
    """(gen, copy, catalyst intact, first contact, parts disturbed by then) for the first copy."""
    Z=ring(cat)|cat; parts=islands(cat); touched=set()
    q=C.Pattern.from_coords(list(pat|cat),R); hit=None
    for g in range(gens+1):
        s=set(q.coords())
        for i,c in enumerate(parts):
            if i not in touched and any((p in s)!=(p in c) for p in ring(c)|c): touched.add(i)
        if hit is None and touched: hit=g
        if hit is not None and g>hit:
            if g>hit+WIN: return None
            cs=find_copies(s)
            if cs:
                intact=all((c2 in s)==(c2 in cat) for c2 in Z)
                return g, cs[0], intact, hit, [parts[i] for i in touched]
        q.step(1)
    return None
seen=set(); n=0; hits=0
for fn in sys.argv[1:]:
    L=open(fn).read().split("\n")
    for i,l in enumerate(L):
        if not l.startswith("Completed:"): continue
        full=set(C.parse_rle(f"x=0,y=0,rule={R}\n"+L[i+1]).coords())
        ir=ring(I0)-I0   # the input starts isolated: match it with an empty ring
        a=min(I0); d=next(((fx-a[0],fy-a[1]) for fx,fy in full
                         if all((x+fx-a[0],y+fy-a[1]) in full for x,y in I0)
                         and not any((x+fx-a[0],y+fy-a[1]) in full for x,y in ir)), None)
        if d is None: continue
        pat={(x+d[0],y+d[1]) for x,y in I0}; cat=full-pat
        key=norm(cat) if cat else None
        k2=frozenset((x-d[0],y-d[1]) for x,y in cat)
        if k2 in seen: continue
        seen.add(k2); n+=1
        q=C.Pattern.from_coords(list(cat),R); q.step(1); still=set(q.coords())==cat
        if not still: continue
        r=first_copy(pat, cat, 400)
        if not r: continue
        g,(dx,dy,t),intact,hit,used=r
        core=frozenset().union(*used)
        r2=first_copy(pat, set(core), g)            # only the parts that took part
        if not r2 or r2[0]!=g: continue
        Zc=ring(core)|core; q=C.Pattern.from_coords(list(pat|core),R); streak=0; back=None
        for gg in range(hit+WIN+15):
            s_=set(q.coords())
            streak = streak+1 if gg>hit and all((c2 in s_)==(c2 in core) for c2 in Zc) else 0
            if streak>=10: back=gg-9; break
            q.step(1)
        if back is None: continue
        intact=r2[2]
        hits+=1
        o=[q for q,tt in ORIENTS if tt==t][0]
        reuse=len({(u+dx,v+dy) for u,v in o} & cat)
        print(f"{fn}: catalyst {len(core)} cells (of {len(cat)}), first contact gen {hit}; copy at gen {g}, "
              f"offset ({dx-d[0]},{dy-d[1]}) orient {t}, catalyst {'intact then' if intact else 'still recovering then'}, back from gen {back}"
              + (f"; {reuse} of the copy's cells were catalyst cells (pre-placed)" if reuse else ""))
        print("   ", C.Pattern.from_coords(list(pat|core),R).to_rle().split("\n",1)[1].replace("\n",""))
print(f"# {n} distinct completed catalysts checked, {hits} produce a copy", file=sys.stderr)
