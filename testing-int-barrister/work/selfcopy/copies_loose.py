# usage: copies.py barrister.out [...]
# For every Completed result: re-simulate with casim (B2i3/S23-c) for 400 gens and look, after the
# catalyst's first contact, for an isolated copy of the input pattern (any of 8 orientations, every
# cell of its 1-cell ring dead). Prints the catalysts that give one, with gen, position, orientation,
# whether the catalyst is intact then and at gen 400.
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B2i3/S23-c"; BODY="2bo$b3o$2o2bo2$4b2o$4b2o$4b2o!"
P0=set(C.parse_rle(f"x=0,y=0,rule={R}\n{BODY}").coords())
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
seen=set(); n=0; hits=0
for fn in sys.argv[1:]:
    L=open(fn).read().split("\n")
    for i,l in enumerate(L):
        if not l.startswith("Completed:"): continue
        full=set(C.parse_rle(f"x=0,y=0,rule={R}\n"+L[i+1]).coords())
        a=min(P0); d=next(((fx-a[0],fy-a[1]) for fx,fy in full
                         if all((x+fx-a[0],y+fy-a[1]) in full for x,y in P0)), None)
        if d is None: continue
        pat={(x+d[0],y+d[1]) for x,y in P0}; cat=full-pat
        key=norm(cat) if cat else None
        k2=frozenset((x-d[0],y-d[1]) for x,y in cat)
        if k2 in seen: continue
        seen.add(k2); n+=1
        q=C.Pattern.from_coords(list(cat),R); q.step(1); still=set(q.coords())==cat
        Z=ring(cat)|cat
        q=C.Pattern.from_coords(list(full),R); hit=None; found=[]
        for g in range(401):
            s=set(q.coords())
            if hit is None and any((c in s)!=(c in cat) for c in Z): hit=g
            if hit is not None and g>hit:
                for c in find_copies(s):
                    found.append((g,c,all((c2 in s)==(c2 in cat) for c2 in Z)))
                if found: break
            q.step(1)
        if found:
            hits+=1
            q.step(max(0,400-found[0][0]))
            end=set(q.coords()); catend=all((c in end)==(c in cat) for c in Z)
            g,(dx,dy,t),catok=found[0]
            print(f"{fn}: catalyst {len(cat)} cells (still={still}), first contact gen {hit}; copy at gen {g}, "
                  f"offset ({dx-d[0]},{dy-d[1]}) orient {t}, catalyst intact then={catok}, at 400={catend}")
            print("   ", C.Pattern.from_coords(list(full),R).to_rle().split("\n",1)[1].replace("\n",""))
print(f"# {n} distinct completed catalysts checked, {hits} produce a copy", file=sys.stderr)
