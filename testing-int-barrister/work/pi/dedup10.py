# usage: dedup.py barrister.out base.rle  -> groups results by (interacting catalyst cells, outcome)
import sys, os, collections; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
R="B36/S23"
def S(q): return set(map(tuple,q.coords() if callable(q.coords) else q.coords))
def run(s,g):
    q=C.Pattern.from_coords(list(s),R); q.step(g); return S(q)
def zoi(s,r=1): return {(x+a,y+b) for x,y in s for a in range(-r,r+1) for b in range(-r,r+1)}
base=S(C.parse_rle(open(sys.argv[2]).read()))
gl={(x,y) for x,y in base if y<=2} if len(base)>7 else set()
stills=base-gl if gl else set()
L=open(sys.argv[1]).read().split("\n")
groups=collections.OrderedDict(); bad=collections.Counter()
for i,l in enumerate(L):
    if not l.startswith("Completed:") or not L[i+1].strip() or L[i+1].strip()=="!": continue
    full=S(C.parse_rle("x = 0, y = 0, rule = "+R+"\n"+L[i+1]))
    bx,by=min(base)
    for fx,fy in full:
        dx,dy=fx-bx,fy-by
        if {(x+dx,y+dy) for x,y in base}<=full: break
    else: bad["noalign"]+=1; continue
    cat={(x-dx,y-dy) for x,y in full}-base
    if run(cat,1)!=cat: bad["cat unstable"]+=1; continue
    # keep if, after the first interaction, the catalyst is intact for >= 10 consecutive gens
    pat=base|cat; Z=zoi(cat); E=set(); s=pat; hit=False; streak=0; okend=None
    for g in range(300):
        d={c for c in Z if (c in s)!=(c in cat)}
        if d: hit=True; streak=0; E|=d
        elif hit:
            streak+=1
            if streak>=10: okend=g; break
        s=run(s,1)
    if okend is None: bad["no 10-gen recovery"]+=1; continue
    if stills:
        s2=run(pat, max(okend+1, 45))   # stills are disturbed by the pi-forming reaction itself until ~gen 42
        if not all(((c in s2)==(c in stills)) for c in zoi(stills)): bad["stills damaged @max(recovery,45)"]+=1; continue
        s=s2
    inter=frozenset(cat & zoi(E,2))
    out=frozenset(s-cat)
    key=(inter,out)
    groups.setdefault(key,[]).append(cat)
print(sys.argv[1], "valid", sum(len(v) for v in groups.values()), "distinct reactions", len(groups), dict(bad))
for (inter,out),v in sorted(groups.items(), key=lambda kv:(len(kv[0][0]),len(kv[0][1])))[:int(os.environ.get("SHOW","12"))]:
    small=min(v,key=len)
    q=C.Pattern.from_coords(list(small),R)
    print(f"  x{len(v):<4} interacting cells {len(inter):2d}, debris pop {len(out):3d}, smallest cat pop {len(small)} : {q.to_rle(pos=True).split(chr(10),1)[1].replace(chr(10),'')}")

# optional: write a summary RLE (pi/base + smallest catalyst per distinct reaction), 8 per row
if len(sys.argv)>3:
    cells=set(); i=0
    for (inter,out),v in sorted(groups.items(), key=lambda kv:(len(min(kv[1],key=len)),len(kv[0][1]))):
        small=min(v,key=len); p=base|small
        mx=min(x for x,y in p); my=min(y for x,y in p)
        ox=(i%8)*40; oy=(i//8)*40
        cells|={(x-mx+ox,y-my+oy) for x,y in p}; i+=1
    if cells:
        open(sys.argv[3],"w").write(C.Pattern.from_coords(list(cells),R).to_rle())
