# usage: conduitcheck.py rule object_rle out1 [...]: find results where an isolated copy of the object
# (any orientation, nothing else within 2 cells) appears after the catalyst's first hit.
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from align import offset
R=sys.argv[1]
def S(q): return set(map(tuple,q.coords() if callable(q.coords) else q.coords))
def zoi(s,r=1): return {(x+a,y+b) for x,y in s for a in range(-r,r+1) for b in range(-r,r+1)}
obj=S(C.parse_rle(f"x = 0, y = 0, rule = {R}\n{sys.argv[2]}"))
def norm(c):
    mx=min(x for x,y in c); my=min(y for x,y in c); return frozenset((x-mx,y-my) for x,y in c)
# any phase of the object's own evolution before it settles (gens 0-5), any orientation
PHASES=int(os.environ.get("PHASES","6"))
forms=set(); q0=C.Pattern.from_coords(list(obj),R)
for ph in range(PHASES):
    o=S(q0)
    for t in range(8):
        c=set()
        for x,y in o:
            if t&4: x,y=y,x
            if t&1: x=-x
            if t&2: y=-y
            c.add((x,y))
        forms.add(norm(c))
    q0.step(1)
def clusters(s):
    left=set(s); out=[]
    while left:
        st=[left.pop()]; comp=set(st)
        while st:
            x,y=st.pop()
            for c in [c for c in left if max(abs(c[0]-x),abs(c[1]-y))<=2]: left.discard(c); comp.add(c); st.append(c)
        out.append(comp)
    return out
hits=[]; tot=0
for f in sys.argv[3:]:
    L=open(f).read().split("\n")
    for i,l in enumerate(L):
        if not (l.startswith("Completed:") and L[i+1].strip()): continue
        tot+=1
        full=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+L[i+1]))
        off=offset(L, i, obj)
        if off is None: continue
        dx,dy=off
        cat={(x-dx,y-dy) for x,y in full}-obj
        q=C.Pattern.from_coords(list(obj|cat),R); Z=zoi(cat); first=None; found=None
        for g in range(200):
            s=S(q)
            if first is None and any(((c in s)!=(c in cat)) for c in Z): first=g
            if first is not None and g>first:
                for comp in clusters(s-cat if all(((c in s)==(c in cat)) for c in cat) else s):
                    if norm(comp) in forms and not (zoi(comp,2)&(s-comp)) :
                        found=(g,min(comp)); break
            if found: break
            q.step(1)
        if found:
            q2=C.Pattern.from_coords(list(obj|cat),R); q2.step(300); s3=S(q2)
            hits.append((f,len(cat),found,cat<=s3, C.Pattern.from_coords(list(obj|cat),R).to_rle().split("\n",1)[1].replace("\n","")))
print(f"checked {tot} results, {len(hits)} produce an isolated copy of the object")
for h in hits[:20]: print(f"  {os.path.basename(h[0])}: catalyst {h[1]} cells, object at gen {h[2][0]} at {h[2][1]}, catalyst intact at 300: {h[3]}  {h[4]}")
