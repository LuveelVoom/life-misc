import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B3/S23-q4i"
here=os.path.dirname(os.path.abspath(__file__))
def S(q): return set(map(tuple,q.coords() if callable(q.coords) else q.coords))
def zoi(s,r=1): return {(x+a,y+b) for x,y in s for a in range(-r,r+1) for b in range(-r,r+1)}
on=S(C.parse_rle(open(os.path.join(here,"rx.rle")).read()))
gl={(2,13),(3,14),(1,15),(2,15),(3,15)}
assert gl<=on
st=on-gl
# lane: glider path from 40 cells upstream until it first diverges from the free glider
up={(x-40,y-40) for x,y in gl}
qa=C.Pattern.from_coords(list(up),R); qb=C.Pattern.from_coords(list(on),R)
lane=set(); hit=None
for g in range(160+60):
    lane|=S(qa)
    if g>=160:
        if S(qb)-st != S(qa) or not st<=S(qb): hit=g-160; break
        qb.step(1)
    qa.step(1)
def allowed(x,y):
    if (x,y) in on or (x,y) in lane: return False
    if any(max(abs(x-a),abs(y-b))<=2 for a,b in st): return False
    return True
def judge(cat,T=250):
    q=C.Pattern.from_coords(list(cat),R); q.step(1)
    if S(q)!=set(cat): return None
    q=C.Pattern.from_coords(list(on|cat),R); Z=zoi(cat); hits=[]; streak=0; rec=None
    for g in range(T):
        s=S(q)
        if any(((c in s)!=(c in cat)) for c in Z): hits.append(g); streak=0; rec=None
        elif hits:
            streak+=1
            if streak>=10 and rec is None: rec=g
        q.step(1)
    s=S(q)
    if not hits or hits[0]<hit or rec is None: return None
    if not all(((c in s)==(c in st)) for c in zoi(st)): return None
    cat_ok=all(((c in s)==(c in cat)) for c in Z)
    return (hits[0],hits[-1],len(s-st-cat),cat_ok)
def leftover(cat,T=250):
    q=C.Pattern.from_coords(list(on|cat),R); q.step(T); return S(q)-st-set(cat)
