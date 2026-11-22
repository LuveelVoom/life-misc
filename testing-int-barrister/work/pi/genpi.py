import os
here=os.path.dirname(os.path.abspath(__file__))
pi={(0,0),(1,0),(2,0),(0,1),(2,1),(0,2),(2,2)}
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
def make(name,xr,yr):
    cells={c:"A" for c in pi}; n=0
    for x in range(*xr):
        for y in range(*yr):
            if any(max(abs(x-a),abs(y-b))<=2 for a,b in pi): continue
            cells[(x,y)]="B"; n+=1
    open(os.path.join(here,name+".toml"),"w").write(f'''# B36/S23: catalysts for a lone pi (3o$obo$obo, moving north). Region {name}.
rule = "B36/S23"
first-active-range = [0, 30]
active-window-range = [0, 30]
min-stable-interval = 8
max-active-cells = 6
stabilise-results = true
pattern = \'\'\'
{rle(cells)}
\'\'\'
''')
    print(name,n)
make("piN",(-9,12),(-12,-2))
make("piE",(4,14),(-10,9))
make("piS",(-8,11),(4,12))
