import sys, re; import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B2i3/S23-c"
act=set([(0, 7),(1, 6),(1, 8),(2, 7),(4, 10),(4, 11),(5, 10),(5, 11),(9, 6),(9, 7),(10, 7),(10, 8),(11, 8),(11, 9)])
def cells_of(rle_body):
    p=C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+rle_body)
    return set(map(tuple,p.coords() if callable(p.coords) else p.coords))
def run(s,g):
    q=C.Pattern.from_coords(list(s),R); q.step(g); return set(map(tuple,q.coords() if callable(q.coords) else q.coords))
out=open(sys.argv[1]).read().split("\n")
ref = None
for i,l in enumerate(out):
    if l.startswith("Completed:"):
        full=cells_of(out[i+1])
        # active = the given A cells: take from the first line of the toml
        yield_=None
        acts=act
        # align: barrister centres the pattern; find translation via the active cells
        best=None
        for (ax,ay) in [min(full)]:
            pass
        # find shift d so that act+d subset full and maximizes
        cands=[(fx-min(act)[0], fy-min(act)[1]) for fx,fy in full]
        for dx,dy in cands:
            sh={(x+dx,y+dy) for x,y in act}
            if sh<=full: best=(dx,dy); break
        A={(x+best[0],y+best[1]) for x,y in act}
        cat=full-A
        still = run(cat,1)==cat
        zoi={(x+a,y+b) for x,y in cat for a in (-1,0,1) for b in (-1,0,1)}
        hits=[g for g in range(0,80) if any(((c in run(full,g))!=(c in cat)) for c in zoi)]
        final=run(full,100); ok100 = all(((c in final)==(c in cat)) for c in zoi)
        print(f"cat pop {len(cat)} still={still} active-gens={hits[:3]}..{hits[-3:] if hits else ''} restored@100={ok100}")
        print("  ", C.Pattern.from_coords(list(cat),R).to_rle().split("\n",1)[1].replace("\n",""))
