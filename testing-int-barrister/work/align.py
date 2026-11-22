# Shared helper: offset of the starting pattern in a Barrister result, taken from the "A" cells
# (the active starting pattern) of the Winner block printed before each "Completed:" block.
# Matching the object's shape inside the completed pattern is unreliable (a glider shape can also
# occur inside a catalyst).
def winner_A(lines, i):
    j=i
    while j>=0 and not lines[j].startswith("Winner:"): j-=1
    if j<0: return None
    body=lines[j+2]; cells=set(); x=y=n=0
    for ch in body:
        if ch.isdigit(): n=n*10+int(ch); continue
        k=n or 1; n=0
        if ch=="$": y+=k; x=0
        elif ch=="!": break
        else:
            if ch=="A": cells|={(x+t,y) for t in range(k)}
            x+=k
    return cells
def offset(lines, i, base):
    A=winner_A(lines, i)
    if not A: return None
    ax,ay=min(A); bx,by=min(base)
    dx,dy=ax-bx,ay-by
    if {(x+dx,y+dy) for x,y in base}!=A: return None
    return dx,dy
