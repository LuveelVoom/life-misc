import sys,os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
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
def make(name,xr,yr,maxact=8):
    cells={c:"A" for c in on}; n=0
    for x in range(*xr):
        for y in range(*yr):
            if allowed(x,y): cells[(x,y)]="B"; n+=1
    open(os.path.join(here,name+".toml"),"w").write(f'''# B3/S23-q4i: catalyst for the glider + block reaction; post-check (common.py) wants the E block and all stills restored.
rule = "B3/S23-q4i"
unsafe-allow-unrepresented-survival = true
first-active-range = [{hitE}, 70]
active-window-range = [0, 40]
min-stable-interval = 10
max-active-cells = {maxact}
stabilise-results = true
pattern = \'\'\'
{rle(cells)}
\'\'\'
''')
    print(name,n)
if len(sys.argv)>1:
    for a in sys.argv[1:]:
        n,x0,x1,y0,y1=a.split(","); make(n,(int(x0),int(x1)),(int(y0),int(y1)))
else:
    make("N",(-12,24),(-18,-1)); make("W",(-12,-1),(-1,16)); make("Ea",(10,24),(-1,16)); make("So",(-1,10),(9,16))
