# usage: check.py barrister.out [...]
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
        a=frozenset(leftover(set(cat),300)); b=leftover(set(cat),301)
        groups[(a, a!=b)].append((f,cat,r))
print("completions",tot,"distinct catalysts",len(seen),"passing",sum(len(v) for v in groups.values()),"distinct leftovers",len(groups))
for (a,changing),v in sorted(groups.items(), key=lambda kv:(bool(zoi(kv[0][0])&lane), len(kv[0][0]))):
    f,cat,r=min(v,key=lambda t:len(t[1]))
    print(f"leftover {len(a):2d} {'(changing) ' if changing else ''}{'TOUCHES LANE' if zoi(a)&lane else 'lane clear'} x{len(v):<3} smallest cat {len(cat)} first {r[0]}: ",
          C.Pattern.from_coords(list(on|cat),R).to_rle().split("\n",1)[1].replace("\n",""))
