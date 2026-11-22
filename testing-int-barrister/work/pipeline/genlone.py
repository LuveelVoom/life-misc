# usage: genlone.py name rule rle_body x0 x1 y0 y1 [maxact]   (region: box minus cells within 2 of the object)
import sys, os
name,rule,body=sys.argv[1:4]; x0,x1,y0,y1=map(int,sys.argv[4:8]); maxact=sys.argv[8] if len(sys.argv)>8 else "6"
cells=set(); x=y=n=0
for ch in body:
    if ch.isdigit(): n=n*10+int(ch); continue
    k=n or 1; n=0
    if ch=="o": cells|={(x+i,y) for i in range(k)}; x+=k
    elif ch=="b": x+=k
    elif ch=="$": y+=k; x=0
    elif ch=="!": break
grid={c:"A" for c in cells}
for xx in range(x0,x1):
    for yy in range(y0,y1):
        if all(max(abs(xx-a),abs(yy-b))>2 for a,b in cells): grid[(xx,yy)]="B"
xs=[a for a,b in grid]; ys=[b for a,b in grid]; X0,Y0,X1,Y1=min(xs),min(ys),max(xs),max(ys); rows=[]
for yy in range(Y0,Y1+1):
    row="".join(grid.get((xx,yy),".") for xx in range(X0,X1+1)).rstrip(".")
    out="";i=0
    while i<len(row):
        j=i
        while j<len(row) and row[j]==row[i]: j+=1
        out+=(str(j-i) if j-i>1 else "")+row[i]; i=j
    rows.append(out)
unsafe = "unsafe-allow-unrepresented-survival = true\n"
open(name+".toml","w").write(f'''# {rule}: catalysts for a lone {body} (region {name}; behind is skipped).
rule = "{rule}"
{unsafe}first-active-range = [0, 30]
active-window-range = [0, 30]
min-stable-interval = 10
max-active-cells = {maxact}
stabilise-results = true
pattern = \'\'\'
x = {X1-X0+1}, y = {Y1-Y0+1}, rule = LifeHistory
{"$".join(rows)}!
\'\'\'
''')
open(name+".base.rle","w").write(f"x = 0, y = 0, rule = {rule}\n{body}\n")
print(name, sum(1 for v in grid.values() if v=="B"))
