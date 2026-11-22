# usage: lonecheck.py rule base.rle summary.rle out1 [out2 ...]
# keeps results whose catalyst is stable in the rule and survives >= 10 gens after the first hit; groups by
# (interacting catalyst cells, leftover at gen 300); writes the smallest catalyst of each group.
import sys, os, collections; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from align import offset
R=sys.argv[1]
def S(q): return set(map(tuple,q.coords() if callable(q.coords) else q.coords))
def zoi(s,r=1): return {(x+a,y+b) for x,y in s for a in range(-r,r+1) for b in range(-r,r+1)}
base=S(C.parse_rle(open(sys.argv[2]).read()))
groups=collections.defaultdict(list); seen=set(); bad=collections.Counter(); tot=0
for f in sys.argv[4:]:
    L=open(f).read().split("\n")
    for i,l in enumerate(L):
        if not (l.startswith("Completed:") and L[i+1].strip()): continue
        tot+=1
        full=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+L[i+1]))
        off=offset(L, i, base)
        if off is None: continue
        dx,dy=off
        cat=frozenset((x-dx,y-dy) for x,y in full)-frozenset(base)
        if cat in seen: continue
        seen.add(cat)
        q=C.Pattern.from_coords(list(cat),R); q.step(1)
        if S(q)!=set(cat): bad["catalyst unstable in rule"]+=1; continue
        q=C.Pattern.from_coords(list(base|cat),R); Z=zoi(cat); E=set(); hit=False; streak=0; ok=False
        for g in range(300):
            s=S(q); d={c for c in Z if (c in s)!=(c in cat)}
            if d: hit=True; streak=0; E|=d
            elif hit:
                streak+=1
                if streak>=10: ok=True
            q.step(1)
        if not ok: bad["no 10-gen recovery"]+=1; continue
        s=S(q); left=frozenset(s-cat) if all(((c in s)==(c in cat)) for c in Z) else None
        key=(frozenset(cat & zoi(E,2)), left)
        groups[key].append(cat)
print(f"{R}: completions {tot}, distinct catalysts {len(seen)}, passing {sum(len(v) for v in groups.values())}, groups {len(groups)}, rejected {dict(bad)}")
cells=set(); notes=[]
for i,((inter,left),v) in enumerate(sorted(groups.items(), key=lambda kv:(kv[0][1] is None, len(kv[0][1] or ()), len(min(kv[1],key=len))))):
    c=min(v,key=len); p=base|c; mx=min(x for x,y in p); my=min(y for x,y in p)
    cells|={(x-mx+(i%8)*50, y-my+(i//8)*50) for x,y in p}
    notes.append(f"{i+1}:cat{len(c)}/" + ("destroyed@300" if left is None else f"left{len(left)}") + f"/x{len(v)}")
hdr=f"#C {R}: one entry per distinct reaction (smallest catalyst), sorted by leftover at gen 300. 8 per row, 50-cell spacing.\n"
for k in range(0,len(notes),8): hdr+="#C "+"  ".join(notes[k:k+8])+"\n"
if cells: open(sys.argv[3],"w").write(hdr+C.Pattern.from_coords(list(cells),R).to_rle())
for n in notes[:15]: print(" ",n)
