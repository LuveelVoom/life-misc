# Cross-check rulespace.py Hensel classes against AlephAlpha/ca-rules ntlife.rs (path argv[1]).
import re, sys, itertools
import rulespace as R
src = open(sys.argv[1]).read()
src = src[src.index("parse_bs!"):src.index("parse_rule!()")]
ca = {}
for blk in re.finditer(r"'(\d)' => \{(.*?)\n        \}", src, re.S):
    n = int(blk.group(1))
    for m in re.finditer(r"'([a-z])' => \[([^\]]*)\]", blk.group(2)):
        ca[(n, m.group(1))] = {int(v, 16) for v in m.group(2).split(",")}
def to9(b):   # ca-rules: 8 bits row-major MSB=NW, centre excluded -> MAP 9-bit, centre 0
    return ((b & 0xE0) << 1) | ((b & 0x10) << 1) | (b & 0x0F)
mine = {k: v for k, v in R.CLASSES.items()}
theirs = {k: {to9(b) for b in v} for k, v in ca.items()}
assert set(mine) == set(theirs), set(mine) ^ set(theirs)
bad = [k for k in mine if mine[k] != theirs[k]]
print("classes:", len(mine), "mismatches:", bad)
# partition check: every centre-0 index in exactly one class with right count
cover = {}
for k, v in mine.items():
    for x in v:
        assert x & 0x10 == 0 and bin(x).count("1") == k[0], (k, x)
        assert x not in cover; cover[x] = k
assert len(cover) == 256
print("transitions:", R.TNAMES)
print("id 0      ->", R.rule_string(0), R.rule_id(R.MIN_RULE))
print("id 0xffff ->", R.rule_string(0xffff), hex(R.rule_id(R.MAX_RULE)))
# Life sanity: B3/S23 table == totalistic
t = R.table_from_sets(*R.parse_hensel("B3/S23"))
ok = all(t[i] == (1 if (bin(i & 0x1ef).count("1") == 3 or (i & 0x10 and bin(i & 0x1ef).count("1") == 2)) else 0) for i in range(512))
print("B3/S23 totalistic:", ok, R.map_string(t))
# round trip string<->id for random ids
import random
for rid in random.Random(0).sample(range(65536), 2000):
    assert R.rule_id(R.rule_string(rid)) == rid
print("round trip ok")
