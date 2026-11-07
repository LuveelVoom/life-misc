#!/usr/bin/env python3
"""period.py FILE.rle [RULE]: period and displacement of the pattern (rule from the header by default)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import minirule as M

cells, rule = M.parse_rle(open(sys.argv[1]).read())
rule = sys.argv[2] if len(sys.argv) > 2 else rule
p = M.period(cells, M.parse(rule), 1000)
print("no period within 1000 gens" if p is None else f"{p['kind']}, period {p['period']}, dx {p['dx']}, dy {p['dy']}")
