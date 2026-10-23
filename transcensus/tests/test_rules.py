#!/usr/bin/env python3
"""Rule parser check: transcensus-int --rule-info must give casim's table (casim follows Golly's liferules.cpp)."""
import os, random, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "casim"))
import numpy as np
import casim as C

EXE = os.path.join(HERE, "..", "transcensus-int")

def ours(rule):
    p = subprocess.run([EXE, "-r", rule, "--rule-info"], capture_output=True, text=True)
    return p.stdout.split("\n")[0] if p.returncode == 0 else None

def ref(rule):
    try:
        return C.map_string(C.parse_rule(rule))
    except C.RuleError:
        return None

rng = random.Random(7)
rules = ["B3/S23", "b3s23", "23/3", "S23/B3", "B3/S", "B/S", "B36/S23", "B2-a3/S23", "B3-cnqy4z/S23-a4i",
         "B2e3ai4n/S12k3-q4ent5", "B2ci3ai4c8/S02ae3eijkq4iz5ar6i7e", "B2-a3_S23", "b2-a3_s23", "B3/S2-i34q",
         "3/23", "B2a2/S", "B4-/S", "B3-/S23", "B3 /S 23", "B0/S8", "B0123478/S01234678",
         # rejected by both
         "B9/S23", "B3x/S23", "B3a-c/S", "B1e/S", "B8e/S", "B3/S23:T100,100", "R1,C2", "B3/S23V", "",
         "B2-a3/S23_", "B3_S23", "-3/S23", "BB3/S23", "B3/S23/S"]
for _ in range(150):
    t = np.zeros(512, np.uint8)
    for orb in {C._orbit(x) for x in range(512)}:
        if rng.random() < 0.3:
            for x in orb: t[x] = 1
    t[0] = 0
    rules.append(C.rule_name(t))
for _ in range(20):
    t = np.array([rng.random() < 0.3 for _ in range(512)], np.uint8)
    rules.append(C.map_string(t))
    rules.append(C.map_string(t) + "==")

bad = 0
for r in rules:
    a, b = ours(r), ref(r)
    if a != b:
        bad += 1
        print("MISMATCH", repr(r), "ours:", a, "casim:", b)
print(f"{len(rules)} rule strings, {bad} mismatches")
sys.exit(1 if bad else 0)
