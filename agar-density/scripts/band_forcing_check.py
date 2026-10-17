# Is the row directly above the band FORCED (period-2 rules, no periodicity assumed)?
# Fix band rows 0,1 over 30 columns; row -1 free; require both transitions valid for row-0 cells in columns 3..26.
import itertools
from pysat.solvers import Cadical153
from pysat.formula import IDPool
def life(c,n): return 1 if (n==3 or (c and n==2)) else 0
TPL=[(b, life(b[0],sum(b[1:]))) for b in itertools.product((0,1),repeat=9)]
band=["XXSYYS","XX.YY."]; code={'.':0,'X':1,'Y':2,'S':3}; W=30
pool=IDPool(); v=lambda t,i,j: pool.id((t%2,i,j)); cls=[]
for i,row in enumerate(band):
    for j in range(W):
        s=code[row[j%6]]; cls+=[[v(0,i,j)] if s&1 else [-v(0,i,j)], [v(1,i,j)] if s&2 else [-v(1,i,j)]]
for t in (0,1):
    for j in range(3,W-3):
        lits=[v(t,0,j)]+[v(t,a,j+b) for a in(-1,0,1) for b in(-1,0,1) if (a,b)!=(0,0)]
        for bits,r in TPL: cls.append([-l if bb else l for l,bb in zip(lits,bits)]+[v(t+1,0,j) if r else -v(t+1,0,j)])
S=Cadical153(bootstrap_with=cls); assert S.solve()
forced=[]
for j in range(4,W-4):
    st=[]
    for val in range(4):
        a=[v(0,-1,j) if val&1 else -v(0,-1,j), v(1,-1,j) if val&2 else -v(1,-1,j)]
        if S.solve(assumptions=a): st.append('.XYS'[val])
    forced.append(''.join(st))
print("possible states of each row -1 cell (columns 4..25):",forced)
