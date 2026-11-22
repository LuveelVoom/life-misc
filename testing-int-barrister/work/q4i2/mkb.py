import sys
exec(open("gen.py").read().split("if len(sys.argv)")[0].replace("__file__",'"gen.py"'))
boat={(8, 10), (9, 9), (9, 11), (10, 10), (10, 11)}
name=sys.argv[1]; x0,x1,y0,y1=map(int,sys.argv[2:6])
cells={c:"A" for c in on}
for x in range(x0,x1):
    for y in range(y0,y1):
        if allowed(x,y): cells[(x,y)]="B"
for c in zoi(boat)-boat:
    if allowed(*c): cells[c]="B"
for c in boat: cells[c]="C"
open(name+".toml","w").write(open("g3.toml").read().split("pattern =")[0]+"exempt-existing = false\npattern = '''\n"+rle(cells)+"\n'''\n")
