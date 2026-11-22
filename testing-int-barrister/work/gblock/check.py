# usage: check.py result.out [...]   re-simulates each Completed result with casim (B3/S23-q4i)
# Reports: catalyst still, first gen the catalyst's ZOI changes, block + empty ring present at gen 300,
# catalyst intact at 300, leftover (everything except catalyst and block) at 300.
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B3/S23-q4i"
GL=[(1,0),(2,1),(0,2),(1,2),(2,2)]; BL=[(6,4),(7,4),(6,5),(7,5)]
def S(q): return set(q.coords())
def run(s,g): q=C.Pattern.from_coords(list(s),R); q.step(g); return S(q)
seen=set(); n=0
for fn in sys.argv[1:]:
    L=open(fn).read().split("\n")
    for i,l in enumerate(L):
        if not l.startswith("Completed:"): continue
        full=S(C.parse_rle(f"x = 0, y = 0, rule = {R}\n"+L[i+1]))
        d=next((fx-1,fy) for fx,fy in full if all((x+fx-1,y+fy) in full for x,y in GL+BL))
        g={(x+d[0],y+d[1]) for x,y in GL}; b={(x+d[0],y+d[1]) for x,y in BL}
        cat=full-g-b
        key=C.Pattern.from_coords([(x-d[0],y-d[1]) for x,y in cat],R).to_rle(pos=True)
        if key in seen: continue
        seen.add(key); n+=1
        still=run(cat,1)==cat
        Z={(x+a,y+c) for x,y in cat for a in (-1,0,1) for c in (-1,0,1)}
        q=C.Pattern.from_coords(list(full),R); hit=None
        for gen in range(300):
            s=S(q)
            if hit is None and any((c in s)!=(c in cat) for c in Z): hit=gen
            q.step(1)
        s=S(q)
        ring={(x+a,y+c) for x,y in b for a in (-1,0,1) for c in (-1,0,1)}
        blockok=all((c in s)==(c in b) for c in ring)
        catok=all((c in s)==(c in cat) for c in Z)
        left=s-cat-b
        print(f"{fn}#{n} catpop={len(cat)} still={still} firsthit={hit} block@300={blockok} cat@300={catok} leftover={len(left)}")
        print("   ", C.Pattern.from_coords(list(full),R).to_rle().split("\n",1)[1].replace("\n",""))
