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
for a in sys.argv[1:]:
    n,x0,x1,y0,y1,ma=a.split(",")
    cells={c:"A" for c in on}; k=0
    for x in range(int(x0),int(x1)):
        for y in range(int(y0),int(y1)):
            if allowed(x,y): cells[(x,y)]="B"; k+=1
    open(os.path.join(here,n+".toml"),"w").write(f'''# B3/S23-q4i: catalyst for the left-going debris (blinker) of the glider reaction.
# A = pattern as given. B = region: off the glider lane, >= 3 cells from existing still lifes.
# First interaction must be after the glider hits (gen {hit}), so incoming gliders pass untouched.
rule = "B3/S23-q4i"
unsafe-allow-unrepresented-survival = true
first-active-range = [{hit}, 80]
active-window-range = [0, 40]
min-stable-interval = 10
max-active-cells = {ma}
stabilise-results = true
pattern = \'\'\'
{rle(cells)}
\'\'\'
''')
    print(n,k)
