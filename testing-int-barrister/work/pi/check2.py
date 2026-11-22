import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B36/S23"
here=os.path.dirname(os.path.abspath(__file__))
p=C.parse_rle(open(os.path.join(here,sys.argv[2] if len(sys.argv)>2 else "rx.rle")).read())
co=set(map(tuple,p.coords() if callable(p.coords) else p.coords))
gl={(x,y) for x,y in co if y<=2}; st=co-gl
def S(s): return set(map(tuple,s.coords() if callable(s.coords) else s.coords))
def run(s,g):
    q=C.Pattern.from_coords(list(s),R); q.step(g); return S(q)
def zoi(s): return {(x+a,y+b) for x,y in s for a in (-1,0,1) for b in (-1,0,1)}
lines=open(sys.argv[1]).read().split("\n")
seen=set()
for i,l in enumerate(lines):
    if not l.startswith("Completed:"): continue
    full=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+lines[i+1]))
    # align original pattern into the output by translation
    for fx,fy in full:
        dx,dy=fx-min(co)[0],fy-min(co)[1]
        if {(x+dx,y+dy) for x,y in co}<=full: break
    cat={(x-dx,y-dy) for x,y in full}-co
    key=tuple(sorted(cat))
    if key in seen: continue
    seen.add(key)
    alls=st|cat
    still = run(cat,1)==cat
    lanegap=min(x-y for x,y in cat)
    pat=co|cat
    hits=[g for g in range(0,400) if not (run(pat,g)&zoi(cat))==(cat if True else 0) and (run(pat,g)&zoi(cat))!=cat]
    end=run(pat,400)
    stills_ok = st<=end and all(((c in end)==(c in st)) for c in zoi(st))
    cat_ok = all(((c in end)==(c in cat)) for c in zoi(cat))
    import itertools
    rng=" ".join(f"{g[0][1]}-{g[-1][1]}" for k,g in ((k,list(g)) for k,g in itertools.groupby(enumerate(hits),lambda t:t[1]-t[0])))
    q=C.Pattern.from_coords(list(cat),R)
    lane=min(abs(x-y+0.5) for x,y in cat); print(f"west={max(x-y for x,y in cat)<=-7} east={min(x-y for x,y in cat)>=6} cat pop {len(cat)} stable={still} min(x-y)={lanegap} catalyst-active-gens={rng[:60]} @400: stills intact={stills_ok} catalyst intact={cat_ok} total pop={len(end)}")
    print("  cat RLE (rel):", q.to_rle(pos=True).replace("\n"," "))
