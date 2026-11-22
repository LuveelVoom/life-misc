# usage: gen.py [--body RLE] NAME [envelope_gens] [margin] [maxact]   (default body: the full pattern)
# B2i3/S23-c: catalysts for 2bo$b3o$2o2bo2$4b2o$4b2o$4b2o (no filter: copies are found by copies.py).
# Region (state B): cells within `margin` of the uncatalysed reaction (gens 0..envelope_gens),
# minus cells within 2 of the gen-0 pattern (no catalyst touches it at the start).
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
BODY="2bo$b3o$2o2bo2$4b2o$4b2o$4b2o!"
if "--body" in sys.argv:
    k=sys.argv.index("--body"); BODY=sys.argv[k+1]; del sys.argv[k:k+2]
name=sys.argv[1]; T=int(sys.argv[2]) if len(sys.argv)>2 else 60
M=int(sys.argv[3]) if len(sys.argv)>3 else 6; maxact=int(sys.argv[4]) if len(sys.argv)>4 else 8
RULE="B2i3/S23-c"
p=C.parse_rle(f"x=0,y=0,rule={RULE}\n{BODY}"); pat=set(p.coords()); ever=set()
for g in range(T+1):
    ever|=set(p.coords()); p.step(1)
def near(c,S,d): return any(max(abs(c[0]-a),abs(c[1]-b))<=d for a,b in S)
xs=[x for x,y in ever]; ys=[y for x,y in ever]
grid={c:"A" for c in pat}; nB=0
for y in range(min(ys)-M,max(ys)+M+1):
    for x in range(min(xs)-M,max(xs)+M+1):
        c=(x,y)
        if c in grid or near(c,pat,2) or not near(c,ever,M): continue
        grid[c]="B"; nB+=1
X0=min(x for x,y in grid); Y0=min(y for x,y in grid); X1=max(x for x,y in grid); Y1=max(y for x,y in grid)
assert X1-X0<58 and Y1-Y0<58, "too big for the 64x64 torus"
rows=[]
for y in range(Y0,Y1+1):
    row="".join(grid.get((x,y),".") for x in range(X0,X1+1)).rstrip(".")
    out="";i=0
    while i<len(row):
        j=i
        while j<len(row) and row[j]==row[i]: j+=1
        out+=(str(j-i) if j-i>1 else "")+row[i]; i=j
    rows.append(out)
open(name+".toml","w").write(f'''# {RULE}: catalysts for {BODY} (region: within {M} of the plain reaction up to gen {T},
# >= 3 cells from the gen-0 pattern; {nB} cells). Copies of the target are looked for afterwards (copies.py).
rule = "{RULE}"
first-active-range = [1, {T}]
active-window-range = [0, 40]
min-stable-interval = 14
max-active-cells = {maxact}
stabilise-results = true
minimise-results = true
pattern-center = [{-X0}, {-Y0}]
pattern = \'\'\'
x = {X1-X0+1}, y = {Y1-Y0+1}, rule = LifeHistory
{"$".join(rows)}!
\'\'\'
''')
print(name, nB, "region cells, box", X1-X0+1, Y1-Y0+1)
