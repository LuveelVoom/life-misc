#!/usr/bin/env python3
"""Rule space between B3-qr5ce6c/S23-a4er5-k and B3-q4ik5cejkr6ce7/S23-a4ejrz56aek7c
(inclusive, 2^16 rules), plus an explosiveness classifier driver for ./soupclass.

Rule ID convention: bit i of the ID adds TRANSITIONS[i] to the minimum rule.
  ID 0x0000 = B3-qr5ce6c/S23-a4er5-k, ID 0xFFFF = B3-q4ik5cejkr6ce7/S23-a4ejrz56aek7c.
Table convention passed to soupclass: 512 entries, index = MAP order (row-major from NW,
weights NW=256 N=128 NE=64 W=32 C=16 E=8 SW=4 S=2 SE=1), sent as 128 hex digits,
entry 0 first, 4 entries per hex digit with the lowest index in the digit's MSB.

Hensel letters: generated the way Golly's liferules.cpp does (representative 9-bit
neighbourhoods for counts 1-4, D4 orbit, counts 5-7 by complementing count 8-n).
Cross-checked against AlephAlpha/ca-rules src/rules/ntlife.rs by test_rulespace.py.
"""
import argparse, base64, csv, math, os, random, statistics, subprocess, sys, time

# ---- Hensel classes (Golly liferules.cpp: entry0..entry3, rule_letters) -------------
_LETTERS = ["ce", "ceaikn", "ceaiknjqry", "ceaiknjqrytwz"]
_ENTRIES = [[1, 2], [5, 10, 3, 40, 33, 68], [69, 42, 11, 7, 98, 13, 14, 70, 41, 97],
            [325, 170, 15, 45, 99, 71, 106, 102, 43, 101, 105, 78, 108]]

def _rot(x):   # Golly rotateBits90Clockwise
    return (((x & 0x4) << 6) | ((x & 0x20) << 2) | ((x & 0x100) >> 2) | ((x & 0x2) << 4)
            | (x & 0x10) | ((x & 0x80) >> 4) | ((x & 0x1) << 2) | ((x & 0x8) >> 2)
            | ((x & 0x40) >> 6))

def _flip(x):  # Golly flipBits
    return ((x & 0x07) << 6) | ((x & 0x1c0) >> 6) | (x & 0x38)

def _orbit(x):
    out, y = set(), x
    for _ in range(4):
        out.add(y); y = _rot(y)
    y = _flip(y)
    for _ in range(4):
        out.add(y); y = _rot(y)
    return out

def letters_for(n):
    if n in (0, 8):
        return "c"
    return "".join(sorted(_LETTERS[min(n, 8 - n) - 1]))

# CLASSES[(n, letter)] = set of 9-bit indices with centre 0
CLASSES = {(0, "c"): {0}, (8, "c"): {0x1ef}}
for n in range(1, 8):
    k = min(n, 8 - n)
    for li, L in enumerate(_LETTERS[k - 1]):
        x = _ENTRIES[k - 1][li] ^ (0x1ef if n > 4 else 0)
        CLASSES[(n, L)] = _orbit(x)

def parse_hensel(s):
    """'B3-qr5ce6c/S23-a4er5-k' -> (birth set, survival set) of (count, letter)."""
    s = s.strip()
    b, sv = s.split("/")
    assert b[0] in "Bb" and sv[0] in "Ss", s
    def part(t):
        out, i = set(), 0
        while i < len(t):
            n = int(t[i]); i += 1
            neg = False
            if i < len(t) and t[i] == "-":
                neg = True; i += 1
            ls = ""
            while i < len(t) and t[i].isalpha():
                ls += t[i]; i += 1
            valid = letters_for(n)
            for c in ls:
                if c not in valid:
                    raise ValueError(f"letter {c} invalid for {n} in {s}")
            if not ls:
                chosen = valid
            elif neg:
                chosen = [c for c in valid if c not in ls]
            else:
                chosen = ls
            out |= {(n, c) for c in chosen}
        return out
    return part(b[1:]), part(sv[1:])

def hensel_string(B, S):
    """Canonical-ish: letters alphabetical; '-' form only if strictly shorter."""
    def part(X):
        out = ""
        for n in range(9):
            valid = letters_for(n)
            have = [c for c in valid if (n, c) in X]
            if not have:
                continue
            if len(have) == len(valid):
                out += str(n)
            else:
                miss = [c for c in valid if c not in have]
                out += str(n) + ("-" + "".join(miss) if len(miss) < len(have) else "".join(have))
        return out
    return "B" + part(B) + "/S" + part(S)

def table_from_sets(B, S):
    t = [0] * 512
    for (n, c) in B:
        for x in CLASSES[(n, c)]:
            t[x] = 1
    for (n, c) in S:
        for x in CLASSES[(n, c)]:
            t[x | 0x10] = 1
    return t

def table_hex(t):
    return "".join("%x" % (t[i] << 3 | t[i + 1] << 2 | t[i + 2] << 1 | t[i + 3])
                   for i in range(0, 512, 4))

def map_string(t):
    """MAP string: 512 bits, index order as above, base64 (Golly/LifeWiki 'MAP' rules)."""
    by = bytes(int("".join(map(str, t[i:i + 8])), 2) for i in range(0, 512, 8))
    return "MAP" + base64.b64encode(by).decode().rstrip("=")

MIN_RULE = "B3-qr5ce6c/S23-a4er5-k"
MAX_RULE = "B3-q4ik5cejkr6ce7/S23-a4ejrz56aek7c"
_BMIN, _SMIN = parse_hensel(MIN_RULE)
_BMAX, _SMAX = parse_hensel(MAX_RULE)
assert _BMIN <= _BMAX and _SMIN <= _SMAX
TRANSITIONS = sorted([("B", x) for x in _BMAX - _BMIN]) + sorted([("S", x) for x in _SMAX - _SMIN])
assert len(TRANSITIONS) == 16, TRANSITIONS
TNAMES = [f"{bs}{n}{c}" for bs, (n, c) in TRANSITIONS]

def rule_sets(rid):
    B, S = set(_BMIN), set(_SMIN)
    for i, (bs, x) in enumerate(TRANSITIONS):
        if rid >> i & 1:
            (B if bs == "B" else S).add(x)
    return B, S

def rule_string(rid):
    return hensel_string(*rule_sets(rid))

def rule_table(rid):
    return table_from_sets(*rule_sets(rid))

def rule_id(s):
    B, S = parse_hensel(s)
    if not (_BMIN <= B <= _BMAX and _SMIN <= S <= _SMAX):
        raise ValueError("rule outside the space: " + s)
    rid = 0
    for i, (bs, x) in enumerate(TRANSITIONS):
        if x in (B if bs == "B" else S):
            rid |= 1 << i
    return rid

# ---- classification ---------------------------------------------------------------
CLASS_ORDER = ["EXPLOSIVE", "MOSTLY EXPLOSIVE", "SOMETIMES EXPLOSIVE", "RARELY EXPLOSIVE",
               "GROWING", "STABLE (long-lived)", "STABLE"]

def classify(r):
    n, p = r["soups"], r["p_explode"]
    if p >= 0.9:
        return "EXPLOSIVE"
    if p >= 0.5:
        return "MOSTLY EXPLOSIVE"
    if p >= 0.1:
        return "SOMETIMES EXPLOSIVE"
    if r["explode"] > 0:
        return "RARELY EXPLOSIVE"
    if r["grow"] / n >= 0.1:
        return "GROWING"
    if r["lifespan_med"] is not None and r["lifespan_med"] >= 500:
        return "STABLE (long-lived)"
    return "STABLE"

def score(r):
    """Single number for ranking, 0..100.  Explosion probability dominates (0..90);
    for rules with no explosions, growth fraction and median lifespan fill 0..10."""
    if r["explode"] > 0:
        return 10 + 90 * r["p_explode"]
    life = min(r["lifespan_med"] or 0, r["maxgen"]) / r["maxgen"]
    grow = min(1.0, (r["grow"] / r["soups"]) / 0.5)     # 50%+ growing soups saturates
    return 5 * grow + 5 * life

def wilson(k, n, z=1.96):
    if n == 0:
        return (0, 1)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)

def run_engine(ids, a):
    exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "soupclass")
    cmd = [exe, "-n", str(a.soups), "-g", str(a.gens), "-s", str(a.seed), "-d", str(a.density),
           "-w", str(a.soup_size), "-N", str(a.grid), "-T", str(a.explode_tiles)]
    if a.threads:
        cmd += ["-t", str(a.threads)]
    inp = "".join(f"{rid:04x} {table_hex(rule_table(rid))}\n" for rid in ids)
    p = subprocess.run(cmd, input=inp, capture_output=True, text=True, check=True)
    res = []
    for line in p.stdout.splitlines():
        f = dict(kv.split("=") for kv in line.split()[1:])
        rid = int(line.split()[0], 16)
        def num(k, typ=float):
            return None if f[k] == "na" else typ(f[k])
        r = dict(id=rid, soups=int(f["n"]), explode=int(f["explode"]), xtrend=int(f["xtrend"]),
                 grow=int(f["grow"]),
                 settled=int(f["settled"]), died=int(f["died"]), active=int(f["active"]),
                 explode_gen_med=num("tx"), speed_med=num("v"), lifespan_med=num("life"),
                 settle_gen_med=num("ts"), maxgen=a.gens)
        r["p_explode"] = r["explode"] / r["soups"]
        r["ci"] = wilson(r["explode"], r["soups"])
        r["class"] = classify(r)
        r["score"] = score(r)
        r["rule"] = rule_string(rid)
        res.append(r)
    return res

def fmt(x, f="{:.0f}"):
    return "–" if x is None else f.format(x)

def write_reports(res, a, out_prefix, elapsed):
    res = sorted(res, key=lambda r: (-r["score"], -(r["speed_med"] or 0), r["id"]))
    with open(out_prefix + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "rule", "class", "score", "p_explode", "ci_lo", "ci_hi", "explode",
                    "explode_by_trend", "grow", "settled", "died", "active", "soups", "explode_gen_med",
                    "speed_med_c", "lifespan_med", "settle_gen_med", "map"])
        for r in res:
            w.writerow([f"0x{r['id']:04x}", r["rule"], r["class"], f"{r['score']:.1f}",
                        f"{r['p_explode']:.3f}", f"{r['ci'][0]:.3f}", f"{r['ci'][1]:.3f}",
                        r["explode"], r["xtrend"], r["grow"], r["settled"], r["died"], r["active"], r["soups"],
                        fmt(r["explode_gen_med"]), fmt(r["speed_med"], "{:.3f}"),
                        fmt(r["lifespan_med"]), fmt(r["settle_gen_med"]),
                        map_string(rule_table(r["id"]))])
    L = []
    L.append(f"# Explosiveness of {len(res)} rules\n")
    L.append(f"Space: `{MIN_RULE}` (ID 0x0000) to `{MAX_RULE}` (ID 0xFFFF), inclusive. "
             f"Bit i of the ID adds transition i: "
             + ", ".join(f"{i}:{t}" for i, t in enumerate(TNAMES)) + ".\n")
    if getattr(a, "exclude", None):
        L.append("Excluded: rules containing any of " + ", ".join(a.exclude) + ".\n")
    L.append(f"Settings: {a.soups} soups/rule, {a.soup_size}×{a.soup_size} at {a.density}%, "
             f"max {a.gens} gens, grid {a.grid}², explosion = {a.explode_tiles} occupied "
             f"16×16 tiles (or ≥ ¼ of that with near-quadratic growth over two time-doublings), "
             f"seed {a.seed}. Run time {elapsed:.0f} s.\n")
    L.append("Scale: score = 10 + 90·P(explode) for rules with any exploding soup; for rules with "
             "none, 0–10 from the growing fraction and median lifespan. Ties broken by speed.\n")
    L.append("## Summary\n")
    L.append("| class | rules |\n|---|---:|")
    for c in CLASS_ORDER:
        k = sum(r["class"] == c for r in res)
        if k:
            L.append(f"| {c} | {k} |")
    L.append("")
    if len(res) >= 20:
        L.append("## Marginal effect of each transition (this sample)\n")
        L.append("Mean P(explode) of sampled rules with vs without the transition. "
                 "Correlational, sample-dependent.\n")
        L.append("| bit | transition | with | without | Δ |\n|---:|---|---:|---:|---:|")
        rows = []
        for i, t in enumerate(TNAMES):
            on = [r["p_explode"] for r in res if r["id"] >> i & 1]
            off = [r["p_explode"] for r in res if not r["id"] >> i & 1]
            if on and off:
                rows.append((statistics.mean(on) - statistics.mean(off), i, t,
                             statistics.mean(on), statistics.mean(off)))
        for d, i, t, m1, m0 in sorted(rows, reverse=True):
            L.append(f"| {i} | {t} | {m1:.2f} | {m0:.2f} | {d:+.2f} |")
        L.append("")
    L.append("## All rules (most explosive first)\n")
    L.append("Columns: P = fraction of soups that exploded (95% Wilson interval); t_x = median "
             "generation at which an explosion was detected; v = median approximate "
             "expansion speed of exploding soups (c, from occupied-area growth); grow = soups still "
             "growing at the end without exploding; settled = became periodic or died; life = "
             "median lifespan of non-exploding soups.\n")
    L.append("| ID | rule | class | score | P (95% CI) | t_x | v | grow | settled | life |")
    L.append("|---|---|---|---:|---|---:|---:|---:|---:|---:|")
    for r in res:
        n = r["soups"]
        L.append(f"| 0x{r['id']:04x} | `{r['rule']}` | {r['class']} | {r['score']:.0f} | "
                 f"{r['p_explode']:.2f} ({r['ci'][0]:.2f}–{r['ci'][1]:.2f}) | "
                 f"{fmt(r['explode_gen_med'])} | {fmt(r['speed_med'], '{:.2f}')} | "
                 f"{r['grow']}/{n} | {r['settled']}/{n} | {fmt(r['lifespan_med'])} |")
    with open(out_prefix + ".md", "w") as fh:
        fh.write("\n".join(L) + "\n")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--random", type=int, metavar="K", help="K distinct random rule IDs")
    g.add_argument("--ids", nargs="+", help="rule IDs (hex like 0x1a2b) or Hensel strings")
    g.add_argument("--all", action="store_true", help="all 65536 rules")
    g.add_argument("--info", nargs="+", help="print string/ID/MAP for IDs or strings and exit")
    ap.add_argument("--sample-seed", type=int, default=1)
    ap.add_argument("--soups", type=int, default=64)
    ap.add_argument("--gens", type=int, default=6000)
    ap.add_argument("--density", type=int, default=50)
    ap.add_argument("--soup-size", type=int, default=16)
    ap.add_argument("--grid", type=int, default=1024)
    ap.add_argument("--explode-tiles", type=int, default=512)
    ap.add_argument("--seed", type=int, default=12345, help="soup seed (same soups for every rule)")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=256, help="rules per engine call")
    ap.add_argument("--out", default="explosiveness")
    ap.add_argument("--exclude", nargs="+", default=[], metavar="T",
                    help="skip rules containing any of these transitions, e.g. B3r B4i")
    a = ap.parse_args()
    bad = [t for t in a.exclude if t not in TNAMES]
    if bad:
        ap.error(f"unknown transition(s) {bad}; choose from {TNAMES}")
    xmask = sum(1 << TNAMES.index(t) for t in a.exclude)

    def to_id(s):
        return int(s, 16) if s.lower().startswith("0x") or all(c in "0123456789abcdefABCDEF" for c in s) and len(s) <= 4 else rule_id(s)
    if a.info:
        for s in a.info:
            rid = to_id(s)
            print(f"0x{rid:04x}  {rule_string(rid)}  {map_string(rule_table(rid))}")
        return
    pool = [r for r in range(65536) if not r & xmask]
    if a.random:
        ids = random.Random(a.sample_seed).sample(pool, a.random)
    elif a.all:
        ids = pool
    else:
        ids = [r for r in (to_id(s) for s in a.ids) if not r & xmask]
    t0, res = time.time(), []
    for i in range(0, len(ids), a.chunk):
        res += run_engine(ids[i:i + a.chunk], a)
        # reports rewritten after every chunk, so an interrupted run keeps what it has done
        write_reports(res, a, a.out, time.time() - t0)
        print(f"{len(res)}/{len(ids)} rules, {time.time() - t0:.0f} s", file=sys.stderr)
    print(f"wrote {a.out}.md and {a.out}.csv", file=sys.stderr)

if __name__ == "__main__":
    main()
