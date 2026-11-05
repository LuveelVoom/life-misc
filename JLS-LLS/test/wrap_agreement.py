import sys
sys.path.insert(0,sys.argv[1]+'/lls')
from src.formatting import jdf_wrap_source,JDF_TRANSLATIONS
from src.SearchPattern import SearchPattern
import src.rules
w,h=5,5
ok=0
for tr,lst in JDF_TRANSLATIONS.items():
  for R,D in [(0,0),(2,-1),(-1,3)]:
    g=[[[f"c{t}_{y}_{x}" for x in range(w)] for y in range(h)] for t in range(3)]
    sp=SearchPattern(g,ignore_transition=[[[False]*w for _ in range(h)] for _ in range(3)],background_grid=[[["0"]]],background_ignore_transition=[[[False]]],rule=src.rules.rule_from_rulestring("B3/S23"))
    pairs=sp.cell_pairs_from_transformation([(lst or ["RO0"])[0],R,D,2])
    want=set()
    for y in range(h):
      for x in range(w):
        sx,sy=jdf_wrap_source(x,y,w,h,R,D,tr)
        if 0<=sx<w and 0<=sy<h: want.add((f"c2_{y}_{x}",f"c0_{sy}_{sx}"))
    got={(a,b) for a,b in pairs if a.startswith("c2_") and b.startswith("c0_")}
    got={p for p in got if p[0].split('_',1)[0]=='c2'}
    assert got==want,(tr,R,D,sorted(got^want)[:4]); ok+=1
print("wrap agreement",ok,"cases")
