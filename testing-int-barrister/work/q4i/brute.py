import sys, os, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qcheck import *
plain=C.Pattern.from_coords(list(on),R)
# generation at which the glider first disturbs the left block in the plain reaction
q=C.Pattern.from_coords(list(on),R); hitE=None
for g in range(60):
    s=S(q)
    if any(((c in s)!=(c in E)) for c in zoi(E)) : hitE=g; break
    q.step(1)
def judge2(cat,T=200):
    q=C.Pattern.from_coords(list(cat),R); q.step(1)
    if S(q)!=set(cat): return None
    q=C.Pattern.from_coords(list(on|cat),R); Z=zoi(cat); hits=[]; Ehit=False
    for g in range(T):
        s=S(q)
        if any(((c in s)!=(c in cat)) for c in Z): hits.append(g)
        if any(((c in s)!=(c in E)) for c in zoi(E)): Ehit=True
        q.step(1)
    s=S(q)
    if not hits or not Ehit or hits[0] < hitE: return None
    if not all(((c in s)==(c in cat)) for c in Z): return None
    if not all(((c in s)==(c in st)) for c in zoi(st)): return None
    if (s & zoi(D)) - st - cat: return None
    rest=s-st-cat
    return (hits[0],hits[-1],len(rest))
shapes={"block":"2o$2o","beehive":"b2o$o2bo$b2o","loaf":"b2o$o2bo$bobo$2bo","boat":"2o$obo$bo","ship":"2o$obo$b2o",
 "tub":"bo$obo$bo","pond":"b2o$o2bo$o2bo$b2o","eater1":"2o$obo$2bo$2b2o","snake":"2obo$ob2o","carrier":"2o$o2bo$2b2o",
 "longboat":"2o$obo$bobo$2bo","barge":"bo$obo$bobo$2bo","mango":"b2o$o2bo$bo2bo$2b2o","longbarge":"bo$obo$bobo$2bobo$3bo"}
def orients(cells):
    out=set()
    for t in range(8):
        c=cells
        if t&4: c={(y,x) for x,y in c}
        if t&1: c={(-x,y) for x,y in c}
        if t&2: c={(x,-y) for x,y in c}
        mx=min(x for x,y in c); my=min(y for x,y in c)
        out.add(frozenset((x-mx,y-my) for x,y in c))
    return out
jobs=[]
for name,r in shapes.items():
    base=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+r+"!"))
    q=C.Pattern.from_coords(list(base),R); q.step(1)
    if S(q)!=base: print("not stable in rule:",name); continue
    for o in orients(base):
        for x in range(9,27):
            for y in range(-8,10):
                cat={(a+x,b+y) for a,b in o}
                if cat & on: continue
                if any(max(abs(cx-a),abs(cy-b))<=2 for cx,cy in cat for a,b in st): continue
                jobs.append((name,frozenset(cat)))
def work(j):
    r=judge2(set(j[1])); return (j,r) if r else None
if __name__=="__main__":
    print("plain reaction hits block at gen",hitE,"; jobs",len(jobs))
    with mp.Pool(8) as p:
        res=[r for r in p.imap_unordered(work,jobs,chunksize=50) if r]
    for (name,cat),(f,l,rest) in sorted(res,key=lambda t:(t[1][2],len(t[0][1]))):
        print(name,"first",f,"last",l,"leftover cells",rest, C.Pattern.from_coords(list(on|set(cat)),R).to_rle().split("\n",1)[1].replace("\n",""))
    print("passing",len(res))
