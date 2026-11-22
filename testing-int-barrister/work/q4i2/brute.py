import multiprocessing as mp, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
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
    if S(q)!=base: print("not stable:",name); continue
    for o in orients(base):
        for x in range(-12,24):
            for y in range(-18,16):
                cat={(a+x,b+y) for a,b in o}
                if all(allowed(*c) for c in cat): jobs.append((name,frozenset(cat)))
def work(j):
    r=judge(set(j[1])); return (j,r) if r else None
if __name__=="__main__":
    q=C.Pattern.from_coords(list(on),R); q.step(250); print("no catalyst: leftover", len(S(q)-st), "; E hit at gen", hitE, "; jobs", len(jobs))
    with mp.Pool(10) as p: res=[r for r in p.imap_unordered(work,jobs,chunksize=50) if r]
    res.sort(key=lambda t:(not t[1][3], t[1][2], len(t[0][1])))
    for (name,cat),(f,l,rest,cok) in res[:25]:
        print(f"{name:9s} first {f:3d} last {l:3d} leftover {rest:3d} catalyst-intact@250 {cok}  ", C.Pattern.from_coords(list(on|set(cat)),R).to_rle().split("\n",1)[1].replace("\n",""))
    print("passing",len(res))
