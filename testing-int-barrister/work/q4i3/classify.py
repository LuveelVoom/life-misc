import sys,os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
def classify(cat):
    a=leftover(cat,300); b=leftover(cat,304)
    onlane=bool(zoi(a)&lane) and bool(a)
    moving = a!=b and len(a)==len(b)
    still=leftover(cat,302)==a
    kind = "nothing" if not a else ("glider/moving" if moving and not still and {(x,y) for x,y in a}!={(x,y) for x,y in leftover(cat,302)} and len(a)==5 else ("blinker/osc" if not still else "still"))
    return a,onlane,kind
if __name__=="__main__":
    for line in open(sys.argv[1]):
        if "catalyst-intact" not in line: continue
        r=line.split()[-1]
        full=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+r))
        bx,by=min(on)
        for fx,fy in full:
            dx,dy=fx-bx,fy-by
            if {(x+dx,y+dy) for x,y in on}<=full: break
        cat={(x-dx,y-dy) for x,y in full}-on
        a,onlane,kind=classify(cat)
        print(f"{line.split()[0]:9s} cat {len(cat):2d} leftover {len(a)} {kind:13s} {'TOUCHES LANE' if onlane else 'lane clear'}  leftover cells {sorted(a)[:6]}")
