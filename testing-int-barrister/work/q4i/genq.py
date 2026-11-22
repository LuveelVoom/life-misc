import os
here=os.path.dirname(os.path.abspath(__file__))
on=set(eval(open(os.path.join(here,"on.txt")).read()))
gl={(x,y) for x,y in on if y<=2}; st=on-gl
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
def make(name,maxact,xr,yr):
    cells={c:"A" for c in on}; n=0
    for x in range(*xr):
        for y in range(*yr):
            if (x,y) in on: continue
            if any(max(abs(x-a),abs(y-b))<=2 for a,b in st|gl): continue
            cells[(x,y)]="B"; n+=1
    open(os.path.join(here,name+".toml"),"w").write(f'''# B3/S23-q4i: catalyst on the right; post-check wants the left block restored and no block on D (11-12,9-10).
# A = glider + all still lifes as drawn (D cells empty). B = catalyst region, >= 3 cells from stills and the glider's start.
rule = "B3/S23-q4i"
unsafe-allow-unrepresented-survival = true
first-active-range = [22, 60]
active-window-range = [0, 30]
min-stable-interval = 10
max-active-cells = {maxact}
stabilise-results = true
pattern = \'\'\'
{rle(cells)}
\'\'\'
''')
    print(name,n)
make("q6",6,(9,27),(-8,10))
make("q10",10,(9,27),(-8,10))
