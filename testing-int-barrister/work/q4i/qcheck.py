import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B3/S23-q4i"
here=os.path.dirname(os.path.abspath(__file__))
on=set(eval(open(os.path.join(here,"on.txt")).read()))
gl={(x,y) for x,y in on if y<=2}; st=on-gl
D={(11,9),(11,10),(12,9),(12,10)}; E={(7,10),(7,11),(8,10),(8,11)}
def S(q): return set(map(tuple,q.coords() if callable(q.coords) else q.coords))
def zoi(s): return {(x+a,y+b) for x,y in s for a in (-1,0,1) for b in (-1,0,1)}
def judge(cat, T=200):
    """returns None if fails, else dict"""
    q=C.Pattern.from_coords(list(cat),R); q.step(1)
    if S(q)!=set(cat): return None
    q=C.Pattern.from_coords(list(on|cat),R)
    hits=[]; Z=zoi(cat)
    for g in range(T):
        s=S(q)
        if any(((c in s)!=(c in cat)) for c in Z): hits.append(g)
        q.step(1)
    s=S(q)
    if not hits: return None
    if not all(((c in s)==(c in cat)) for c in Z): return None
    if not all(((c in s)==(c in st)) for c in zoi(st)): return None   # all drawn stills (incl. left block) intact, D empty
    if s & zoi(D) - st - cat: return None
    rest=s-st-cat
    return {"first":hits[0],"last":hits[-1],"rest":len(rest)}
if __name__=="__main__":
    L=open(sys.argv[1]).read().split("\n"); ok=0; tot=0
    for i,l in enumerate(L):
        if l.startswith("Completed:") and L[i+1].strip():
            tot+=1
            full=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+L[i+1])); bx,by=min(on)
            for fx,fy in full:
                dx,dy=fx-bx,fy-by
                if {(x+dx,y+dy) for x,y in on}<=full: break
            cat={(x-dx,y-dy) for x,y in full}-on
            r=judge(cat)
            if r: ok+=1; print(r, C.Pattern.from_coords(list(on|cat),R).to_rle().split("\n",1)[1].replace("\n",""))
    print("checked",tot,"passing",ok)
