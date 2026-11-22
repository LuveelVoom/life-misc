import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../casim"))
import casim as C
here=os.path.dirname(os.path.abspath(__file__))
co=set(map(tuple,(lambda p: p.coords() if callable(p.coords) else p.coords)(C.parse_rle(open(os.path.join(here,"rx2.rle")).read()))))
st={(x,y) for x,y in co if y>2}
def rle(cells):
    xs=[x for x,y in cells]; ys=[y for x,y in cells]
    x0,y0,x1,y1=min(xs),min(ys),max(xs),max(ys); rows=[]
    for y in range(y0,y1+1):
        row="".join(cells.get((x,y),".") for x in range(x0,x1+1)).rstrip(".")
        out="";i=0
        while i<len(row):
            j=i
            while j<len(row) and row[j]==row[i]: j+=1
            out+=(str(j-i) if j-i>1 else "")+row[i]; i=j
        rows.append(out)
    return f"x = {x1-x0+1}, y = {y1-y0+1}, rule = LifeHistory\n"+"$".join(rows)+"!"
def make(name, side, maxact, xr, yr):
    cells={c:"A" for c in co}; n=0
    for x in range(*xr):
        for y in range(*yr):
            d=x-y
            ok=(side in("L","B") and d<=-7) or (side in("R","B") and d>=6)
            if not ok or (x,y) in co: continue
            if any(max(abs(x-a),abs(y-b))<=2 for a,b in st): continue
            cells[(x,y)]="B"; n+=1
    open(os.path.join(here,name+".toml"),"w").write(f'''# B36/S23 pi catalyst search, eater removed (crude INT port of Barrister)
# A = glider + still lifes. B = catalyst region, side {side}: >= 3 cells from the glider lane
# (lane x-y in [-2,1]; west x-y <= -7, east x-y >= 6), >= 3 cells from existing stills.
rule = "B36/S23"
first-active-range = [20, 70]
active-window-range = [0, 30]
min-stable-interval = 8
max-active-cells = {maxact}
stabilise-results = true
pattern = \'\'\'
{rle(cells)}
\'\'\'
''')
    print(name, n, "region cells")
make("L6","L",6,(-12,12),(-6,22))
make("L10","L",10,(-12,12),(-6,22))
make("R6","R",6,(4,24),(-8,14))
make("B8","B",8,(-12,24),(-8,22))
