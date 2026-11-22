# usage: regrid.py file.rle [...]   Re-lays out a summary written as 8 per row on a 50-cell grid
# into a roughly square grid with spacing = largest entry + 12. Entry order (row-major) is kept,
# so the "#C n:..." index lines stay valid. Files with more than 400 entries are split into
# file_part1.rle, file_part2.rle, ... (400 entries each, numbering continues across parts).
import sys, os, math, re
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../casim"))
import casim as C
PER=400
G=int(os.environ.get("GRID","50"))
for f in sys.argv[1:]:
    txt=open(f).read()
    comments=[l for l in txt.split("\n") if l.startswith("#C")]
    p=C.parse_rle(txt); R=p.rule if isinstance(p.rule,str) else (p.rule.name if hasattr(p.rule,"name") else str(p.rule))
    m=re.search(r"rule = (\S+)", txt); R=m.group(1) if m else R
    cells=set(map(tuple,p.coords() if callable(p.coords) else p.coords))
    ents={}
    for x,y in cells: ents.setdefault((y//G)*8+x//G,set()).add((x,y))
    idx=sorted(ents); n=len(idx)
    norm=[]
    for i in idx:
        e=ents[i]; mx=min(x for x,y in e); my=min(y for x,y in e)
        norm.append({(x-mx,y-my) for x,y in e})
    W=max(max(x for x,y in e) for e in norm)+1+12; H=max(max(y for x,y in e) for e in norm)+1+12
    parts=[norm[k:k+PER] for k in range(0,n,PER)]
    for pi,part in enumerate(parts):
        cols=max(1,math.ceil(math.sqrt(len(part)*H/W)))
        out=set()
        for j,e in enumerate(part):
            ox=(j%cols)*W; oy=(j//cols)*H
            out|={(x+ox,y+oy) for x,y in e}
        first=pi*PER+1
        hdr=[]
        for l in comments:
            l=re.sub(r"(Grid: )?8 per row, \d+-cell spacing(, row-major)?\.?", f"{cols} per row, spacing {W}x{H}, row-major.", l)
            hdr.append(l)
        if not any("per row, spacing" in l for l in hdr):
            hdr.insert(0, f"#C {cols} per row, spacing {W}x{H}, row-major.")
        if len(parts)>1: hdr.insert(1, f"#C part {pi+1} of {len(parts)}: entries {first}-{first+len(part)-1}")
        name=f if len(parts)==1 else f.replace(".rle", f"_part{pi+1}.rle")
        open(name,"w").write("\n".join(hdr)+"\n"+C.Pattern.from_coords(list(out),R).to_rle())
        print(f"{name}: {len(part)} entries, {cols} per row, cell {W}x{H}, size {os.path.getsize(name)} bytes")
    if len(parts)>1: os.remove(f)
