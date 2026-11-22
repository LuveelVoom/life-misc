import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B3/S23-q4i"
here=os.path.dirname(os.path.abspath(__file__))
on=set(eval(open(os.path.join(here,"on.txt")).read())); E=set(eval(open(os.path.join(here,"E.txt")).read()))
gl={(x,y) for x,y in on if y<=2 and x<=2}; st=on-gl
def S(q): return set(map(tuple,q.coords() if callable(q.coords) else q.coords))
def zoi(s,r=1): return {(x+a,y+b) for x,y in s for a in range(-r,r+1) for b in range(-r,r+1)}
q=C.Pattern.from_coords(list(on),R); hitE=None
for g in range(60):
    s=S(q)
    if any(((c in s)!=(c in E)) for c in zoi(E)): hitE=g; break
    q.step(1)
def allowed(x,y):
    if (x,y) in on: return False
    if any(max(abs(x-a),abs(y-b))<=2 for a,b in on): return False
    if y<=4 and -6<=x-y<=5: return False      # incoming glider lane (moving SE) with margin
    return True
def judge(cat,T=250):
    q=C.Pattern.from_coords(list(cat),R); q.step(1)
    if S(q)!=set(cat): return None
    q=C.Pattern.from_coords(list(on|cat),R); Z=zoi(cat); hits=[]; streak=0; rec=None
    for g in range(T):
        s=S(q)
        if any(((c in s)!=(c in cat)) for c in Z):
            hits.append(g); streak=0; rec=None
        elif hits:
            streak+=1
            if streak>=10 and rec is None: rec=g
        q.step(1)
    s=S(q)
    if not hits or hits[0]<hitE or rec is None: return None
    if not all(((c in s)==(c in st)) for c in zoi(st)): return None   # E and other stills intact
    cat_ok=all(((c in s)==(c in cat)) for c in Z)
    rest=s-st-cat
    return (hits[0],hits[-1],len(rest),cat_ok)
