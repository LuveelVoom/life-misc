# usage: check.py barrister.out [...]  -> passing results grouped by leftover debris
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
groups=collections.defaultdict(list); tot=0; seen=set()
for f in sys.argv[1:]:
    L=open(f).read().split("\n")
    for i,l in enumerate(L):
        if not (l.startswith("Completed:") and L[i+1].strip()): continue
        tot+=1
        full=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+L[i+1])); bx,by=min(on)
        for fx,fy in full:
            dx,dy=fx-bx,fy-by
            if {(x+dx,y+dy) for x,y in on}<=full: break
        else: continue
        cat=frozenset((x-dx,y-dy) for x,y in full)-frozenset(on)
        if cat in seen: continue
        seen.add(cat)
        r=judge(set(cat))
        if not r: continue
        q=C.Pattern.from_coords(list(on|cat),R); q.step(250); rest=frozenset(S(q)-st-cat)
        groups[rest].append((f,cat,r))
print("completions",tot,"distinct catalysts",len(seen),"passing",sum(len(v) for v in groups.values()),"distinct leftovers",len(groups))
for rest,v in sorted(groups.items(), key=lambda kv:len(kv[0]))[:15]:
    f,cat,r=min(v,key=lambda t:len(t[1]))
    onlane=any(-3<=x-y<=2 for x,y in rest)
    print(f"leftover {len(rest):2d} {'ON LANE ' if onlane else ''}x{len(v):<3} smallest cat {len(cat)} ({os.path.basename(f)}, first {r[0]}, intact@250 {r[3]}): ",
          C.Pattern.from_coords(list(on|cat),R).to_rle().split("\n",1)[1].replace("\n",""))
