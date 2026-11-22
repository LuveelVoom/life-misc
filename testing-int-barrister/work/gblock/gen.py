# usage: gen.py name T [margin] [maxact]
# B3/S23-q4i glider + block; only filter: block restored.
# Region (state B): cells within `margin` of the uncatalysed reaction (gens 0-90), minus
#  - cells within 2 of anything ON before generation T (so nothing touches the glider or the
#    block until the collision is T gens old: no glider eaters / eater-2 style catches),
#  - cells within 2 of the glider's incoming lane.
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
name=sys.argv[1]; T=int(sys.argv[2]); M=int(sys.argv[3]) if len(sys.argv)>3 else 6
maxact=int(sys.argv[4]) if len(sys.argv)>4 else 8
RULE="B3/S23-q4i"
glider=[(1,0),(2,1),(0,2),(1,2),(2,2)]; block=[(6,4),(7,4),(6,5),(7,5)]
p=C.Pattern.from_coords(glider+block, RULE)
early=set(); ever=set()
for g in range(91):
    cs=set(p.coords())
    if g<T: early|=cs
    ever|=cs; p.step(1)
lane={(x-k,y-k) for k in range(0,40) for x,y in glider}
def near(c,S,d): return any(max(abs(c[0]-a),abs(c[1]-b))<=d for a,b in S)
xs=[x for x,y in ever]; ys=[y for x,y in ever]
grid={c:"A" for c in glider+block}
nB=0
for y in range(min(ys)-M,max(ys)+M+1):
    for x in range(min(xs)-M,max(xs)+M+1):
        c=(x,y)
        if c in grid or not near(c,ever,M) or near(c,early,2) or near(c,lane,2): continue
        grid[c]="B"; nB+=1
X0=min(x for x,y in grid); Y0=min(y for x,y in grid); X1=max(x for x,y in grid); Y1=max(y for x,y in grid)
assert X1-X0<56 and Y1-Y0<56, "too big for the 64x64 torus"
rows=[]
for y in range(Y0,Y1+1):
    row="".join(grid.get((x,y),".") for x in range(X0,X1+1)).rstrip(".")
    out="";i=0
    while i<len(row):
        j=i
        while j<len(row) and row[j]==row[i]: j+=1
        out+=(str(j-i) if j-i>1 else "")+row[i]; i=j
    rows.append(out)
# pattern-center = block's top-left cell in pattern coordinates; filter is the block with an OFF ring
cx,cy=6-X0,4-Y0
open(name+".toml","w").write(f'''# {RULE}: glider + block, catalyst must restore the block. Region: within {M} of the
# plain reaction, >= 3 from everything ON before gen {T}, >= 3 from the glider's lane. {nB} cells.
rule = "{RULE}"
first-active-range = [{T}, 90]
active-window-range = [0, 60]
min-stable-interval = 14
max-active-cells = {maxact}
stabilise-results = true
minimise-results = true
pattern-center = [{cx}, {cy}]
pattern = \'\'\'
x = {X1-X0+1}, y = {Y1-Y0+1}, rule = LifeHistory
{"$".join(rows)}!
\'\'\'

# The block recovers:
[[filter]]
filter-gen = 200
filter-type = "EVER"
filter-pos = [-1, -1]
filter = \'\'\'
x = 4, y = 4, rule = LifeHistory
4D$D2CD$D2CD$4D!
\'\'\'
''')
print(name, nB, "cells, box", X1-X0+1, Y1-Y0+1)
