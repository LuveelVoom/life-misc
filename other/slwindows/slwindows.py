#!/usr/bin/env python3
"""Small still lifes containing every 3x3 window a still life can contain.

Works for Life-like and INT rules; default B3/S23. B0 and non-isotropic rules are
rejected: a dead plane must be stable, and classes need the square symmetries.

  python3 slwindows.py classify                 # which of the 512 windows occur; writes windows_<rule>.json
  python3 slwindows.py find 20 20 --mode oriented --time 600 --maxpop 150 --out best.rle
  python3 slwindows.py find 13 13 --mode orbit --rule B3/S23-a5 --incremental

Modes: oriented = every window separately; orbit = at least one window from each
rotation/reflection class. Patterns are stable, not necessarily connected ("strict").
Every window position on the dead plane counts, including those overlapping the box edge.

Window index: Golly MAP, NW 256, N 128, NE 64, W 32, C 16, E 8, SW 4, S 2, SE 1;
x right, y down. SAT encoding is CNF run through kissat (several seeds in parallel).
"""
import argparse, json, os, subprocess, sys, tempfile, time
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType

OFFS = [(dx, dy) for dy in (-1, 0, 1) for dx in (-1, 0, 1)]
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import minirule as M  # noqa: E402


def load_rule(text):
    """Golly-canonical name and the 512 'cell stays as it is' flags."""
    try:
        table = M.parse(text)
    except M.RuleError as e:
        sys.exit(f"{e} (MAP rules aren't supported)")
    name = M.name(table)
    if M.b0(table):
        sys.exit(f"{name}: B0 rule, a dead plane is not stable; not supported")
    return name, [table[i] == ((i >> 4) & 1) for i in range(512)]


def windows_file(name):
    return os.path.join(HERE, "windows_" + name.replace("/", "_") + ".json")


def bits(i):
    return [(i >> (8 - k)) & 1 for k in range(9)]


def orbit_id():
    def g(i): return [[(i >> (8 - 3 * r - c)) & 1 for c in range(3)] for r in range(3)]
    def ix(a): return sum(a[r][c] << (8 - 3 * r - c) for r in range(3) for c in range(3))
    def rot(a): return [[a[2 - c][r] for c in range(3)] for r in range(3)]
    o = [-1] * 512
    for i in range(512):
        if o[i] < 0:
            a = g(i)
            for _ in range(4):
                a = rot(a)
                o[ix(a)] = i
                o[ix([r[::-1] for r in a])] = i
    return o


def stability_clauses(var, centres, stable):
    """var((x, y)) -> int literal, or None for a cell fixed dead."""
    cl = set()
    for x, y in centres:
        lits = [var((x + dx, y + dy)) for dx, dy in OFFS]
        for i in range(512):
            if stable[i]:
                continue
            pat = list(enumerate(bits(i)))
            if any(b and lits[k] is None for k, b in pat):
                continue
            cl.add(tuple(sorted(-lits[k] if b else lits[k] for k, b in pat if lits[k] is not None)))
    return [list(c) for c in cl if c]


def region(r):
    return {(x, y): 1 + (y + r) * (2 * r + 1) + (x + r) for y in range(-r, r + 1) for x in range(-r, r + 1)}


def window_assumps(v, w):
    return [v[(dx, dy)] if b else -v[(dx, dy)] for (dx, dy), b in zip(OFFS, bits(w))]


def classify(a, stable):
    res, wit, proofs = {}, {}, {}
    box = Cadical153(bootstrap_with=None)
    r = a.radius
    v = region(r)
    for c in stability_clauses(v.get, [(x, y) for y in range(-r - 1, r + 2) for x in range(-r - 1, r + 2)], stable):
        box.add_clause(c)
    local = {}
    for k in range(1, 5):                 # stability required only within radius k; outside free
        s = Cadical153()
        vv = region(k + 1)
        for c in stability_clauses(vv.get, [(x, y) for y in range(-k, k + 1) for x in range(-k, k + 1)], stable):
            s.add_clause(c)
        local[k] = (s, vv)
    for w in range(512):
        if not stable[w]:
            res[w] = 0
            proofs[w] = 0
            continue
        if box.solve(assumptions=window_assumps(v, w)):
            m = set(l for l in box.get_model() if l > 0)
            wit[w] = [c for c, l in v.items() if l in m]
            res[w] = 1
            continue
        for k in local:
            s, vv = local[k]
            if not s.solve(assumptions=window_assumps(vv, w)):
                res[w], proofs[w] = 0, k
                break
        else:
            res[w] = None                 # no proof either way
    return res, wit, proofs


def check(cells, rule):
    """Independent of the CNF: stability by simulation, windows by direct scan."""
    per = M.period(cells, M.parse(rule), 4) or {"kind": "no period", "period": "?"}
    live = set(cells)
    seen = {0}
    for x, y in {(x + i, y + j) for x, y in live for i in range(-2, 3) for j in range(-2, 3)}:
        i = 0
        for dx, dy in OFFS:
            i = 2 * i + ((x + dx, y + dy) in live)
        seen.add(i)
    return per, seen


def to_rle(cells, w, h, rule):
    import itertools
    live = set(cells)
    rows = []
    for y in range(h):
        row = "".join("o" if (x, y) in live else "b" for x in range(w)).rstrip("b")
        rows.append("".join(f"{n if n > 1 else ''}{k}" for k, n in
                            ((k, len(list(g))) for k, g in itertools.groupby(row))))
    return f"x = {w}, y = {h}, rule = {rule}\n" + "$".join(rows) + "!\n"


def groups_all(groups):
    return [w for g in groups.values() for w in g]


def incremental(a, cnf, cell, groups, z_of, W, H, rule):
    """Clauses for covered targets are added as we go; phases carry the pattern along."""
    n_cov = len(cnf) - sum(1 for g in groups.values() if 0 not in g)
    base, cover = cnf[:n_cov], cnf[n_cov:]
    s = Cadical153(bootstrap_with=base)
    t0 = time.time()
    done = 0
    for c in sorted(cover, key=len):      # rarest placements first
        if not c:
            print("a target has no possible position")
            return
        s.add_clause(c)
        if not s.solve():
            print(f"UNSAT with {done + 1}/{len(cover)} targets in {time.time() - t0:.1f}s")
            return
        done += 1
        m = set(l for l in s.get_model() if l > 0)
        cells = [k for k, l in cell.items() if l in m]
        print(f"[{time.time() - t0:7.1f}s] {done}/{len(cover)} pop {len(cells)}", flush=True)
        if a.out:
            open(a.out, "w").write(to_rle(cells, W, H, rule))
    per, seen = check(cells, rule)
    miss = sum(1 for g in groups.values() if not (set(g) & seen))
    print(f"DONE pop {len(cells)}, {per['kind']} p{per['period']}, missing {miss}, "
          f"oriented windows {len(seen & set(groups_all(groups)))}/{len(groups_all(groups))}")
    print(to_rle(cells, W, H, rule), end="")


def window_at(live, x, y):
    i = 0
    for dx, dy in OFFS:
        i = 2 * i + ((x + dx, y + dy) in live)
    return i


def repatch(live, free, W, H, allowed, groups, bound, rng, shuffle, budget):
    """Best-effort: new contents for the cells in `free`, everything else fixed, with at
    most `bound` live cells in total and all targets still covered. None if not found."""
    from pysat.solvers import Glucose4
    var = {c: n + 1 for n, c in enumerate(sorted(free))}
    nv = len(var)
    near = {(x + dx, y + dy) for x, y in free for dx, dy in OFFS}
    positions = {c for c in near if -1 <= c[0] <= W and -1 <= c[1] <= H}
    seen = {0}
    for y in range(-1, H + 1):
        for x in range(-1, W + 1):
            if (x, y) not in positions:
                seen.add(window_at(live, x, y))
    need = [g for g in groups if not (set(g) & seen)]
    cnf, z_of = [], {}
    for x, y in positions:
        one = []
        terms = [(var.get((x + dx, y + dy)), (x + dx, y + dy) in live) for dx, dy in OFFS]
        for w in allowed:
            bs = bits(w)
            if any(v is None and b != f for (v, f), b in zip(terms, bs)):
                continue
            nv += 1
            one.append(nv)
            z_of.setdefault(w, []).append(nv)
            cnf += [[-nv, v if b else -v] for (v, _), b in zip(terms, bs) if v is not None]
        if not one:
            return None
        cnf.append(one)
    for g in need:
        c = [z for w in g for z in z_of.get(w, [])]
        if not c:
            return None
        cnf.append(c)
    room = bound - len(live - free)
    if room < 0:
        return None
    if room < len(var):
        enc = CardEnc.atmost(lits=list(var.values()), bound=room, top_id=nv, encoding=EncType.seqcounter)
        cnf += enc.clauses
    with Glucose4(bootstrap_with=cnf) as s:
        if shuffle:
            s.set_phases([v * rng.choice((1, -1)) for v in var.values()])
        s.conf_budget(budget)
        if s.solve_limited() is not True:
            return None
        m = set(l for l in s.get_model() if l > 0)
    return (live - free) | {c for c, v in var.items() if v in m}


def polish_whole(a, rule, live, W, H, rng):
    """Whole box free: the pattern only seeds the solver's phases. Each success tightens
    the cap by --drop. Exact for the box if UNSAT comes back. One thread."""
    from pysat.card import ITotalizer
    from pysat.solvers import Glucose4
    a.width, a.height, a.sym = W, H, False
    cnf, cell, groups, nv = build(a, rule)[:4]
    tot = ITotalizer(lits=list(cell.values()), ubound=len(live) + 1, top_id=nv)
    s = Glucose4(bootstrap_with=cnf + tot.cnf.clauses)
    t0 = time.time()
    pop = len(live)
    print(f"{rule} whole {W}x{H} from pop {pop}, drop {a.drop}, noise {a.noise}", flush=True)
    while time.time() - t0 < a.time:
        cap = pop - a.drop
        s.set_phases([(l if ((c in live) != (rng.random() < a.noise)) else -l) for c, l in cell.items()])
        s.conf_budget(a.budget)
        r = s.solve_limited(assumptions=[-tot.rhs[cap]])
        if r is False:
            print(f"[{time.time() - t0:7.1f}s] UNSAT: no pattern of at most {cap} in this box", flush=True)
            return
        if r is None:
            continue
        m = set(l for l in s.get_model() if l > 0)
        live = {c for c, l in cell.items() if l in m}
        pop = len(live)
        per, seen = check(sorted(live), rule)
        miss = sum(1 for g in groups.values() if not (set(g) & seen))
        print(f"[{time.time() - t0:7.1f}s] pop {pop}; {per['kind']} p{per['period']}, missing {miss}", flush=True)
        if a.out and per["kind"] == "still life" and miss == 0:
            open(a.out, "w").write(to_rle(sorted(live), W, H, rule))
    print(f"END pop {pop}")


def polish(a):
    """Large-neighbourhood search: free a random patch, try to lose a cell; accept
    equal-population changes so the search can move. Heuristic, one thread."""
    import random
    rng = random.Random(a.seed)
    start, _ = M.parse_rle(open(a.start).read())
    rule, _ = load_rule(a.rule)
    d = json.load(open(windows_file(rule)))
    allowed = sorted(d["realizable"] + d["open"])
    orb = orbit_id()
    grouped = {}
    for w in d["realizable"]:
        grouped.setdefault(w if a.mode == "oriented" else orb[w], []).append(w)
    groups = list(grouped.values())
    x0, y0 = min(x for x, _ in start), min(y for _, y in start)
    w, h = max(x for x, _ in start) - x0 + 1, max(y for _, y in start) - y0 + 1
    live = {(x - x0 + a.pad, y - y0 + a.pad) for x, y in start}
    W, H = w + 2 * a.pad, h + 2 * a.pad
    if a.whole:
        return polish_whole(a, rule, live, W, H, rng)
    t0 = time.time()
    best = len(live)
    it = moves = 0
    print(f"{rule} start pop {best}, box {W}x{H}, {len(groups)} targets", flush=True)
    while time.time() - t0 < a.time:
        it += 1
        n = rng.randint(a.minpatch, a.maxpatch)
        px, py = rng.randint(0, max(0, W - n)), rng.randint(0, max(0, H - n))
        free = {(x, y) for y in range(py, min(H, py + n)) for x in range(px, min(W, px + n))}
        new = repatch(live, free, W, H, allowed, groups, len(live) - 1, rng, False, a.budget)
        if new is None:
            new = repatch(live, free, W, H, allowed, groups, len(live), rng, True, a.budget)
            if new is not None and new != live:
                moves += 1
                live = new
        else:
            live = new
            if len(live) < best:
                best = len(live)
                per, seen = check(sorted(live), rule)
                miss = sum(1 for g in groups if not (set(g) & seen))
                print(f"[{time.time() - t0:7.1f}s] pop {best} after {it} patches, {moves} sideways; "
                      f"{per['kind']} p{per['period']}, missing {miss}", flush=True)
                if a.out and per["kind"] == "still life" and miss == 0:
                    open(a.out, "w").write(to_rle(sorted(live), W, H, rule))
        if it % 500 == 0:
            print(f"[{time.time() - t0:7.1f}s] {it} patches, {moves} sideways, current {len(live)}, best {best}", flush=True)
    print(f"END best pop {best}, {it} patches, {moves} sideways")


def build(a, rule):
    path = windows_file(rule)
    if not os.path.exists(path):
        sys.exit(f"run: slwindows.py classify --rule {rule}")
    d = json.load(open(path))
    real = d["realizable"]
    allowed = sorted(real + d["open"])    # open windows may occur, so they stay legal
    orb = orbit_id()
    groups = {}
    for w in real:
        groups.setdefault(w if a.mode == "oriented" else orb[w], []).append(w)
    W, H = a.width, a.height
    def canon(x, y):
        if not a.sym:
            return x, y
        x, y = min(x, W - 1 - x), min(y, H - 1 - y)
        return min(x, y), max(x, y)
    if a.sym and W != H:
        sys.exit("--sym needs a square box")
    cell = {(x, y): 1 + canon(x, y)[1] * W + canon(x, y)[0] for y in range(H) for x in range(W)}
    nv = W * H
    cnf = []
    z_of = {}
    for y in range(-1, H + 1):
        for x in range(-1, W + 1):
            if canon(x, y) != (x, y):
                continue                  # the rest follows by symmetry
            lits = [cell.get((x + dx, y + dy)) for dx, dy in OFFS]
            one = []
            for i in allowed:
                bs = bits(i)
                if any(b and l is None for b, l in zip(bs, lits)):
                    continue
                nv += 1
                one.append(nv)
                z_of.setdefault(i, []).append(nv)
                for b, l in zip(bs, lits):
                    if l is not None:
                        cnf.append([-nv, l if b else -l])
            cnf.append(one)               # stability: each position shows a stable window
    for g in groups.values():
        if 0 in g:
            continue
        cnf.append([z for w in g for z in z_of.get(w, [])])
    return cnf, cell, groups, nv, real, allowed, z_of


def find(a):
    rule, stable = load_rule(a.rule)
    cnf, cell, groups, nv, real, allowed, z_of = build(a, rule)
    W, H = a.width, a.height
    if a.maxpop:
        enc = CardEnc.atmost(lits=list(cell.values()), bound=a.maxpop, top_id=nv, encoding=EncType.seqcounter)
        nv = max(nv, enc.nv)
        cnf = enc.clauses + cnf
    if a.incremental:
        return incremental(a, cnf, cell, groups, z_of, W, H, rule)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".cnf", delete=False, dir=a.tmp)
    tmp.write(f"p cnf {nv} {len(cnf)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cnf))
    tmp.close()
    print(f"{rule} {W}x{H} {a.mode}: {len(groups)} targets, {nv} vars, {len(cnf)} clauses", flush=True)
    t0 = time.time()
    outs = [tempfile.TemporaryFile("w+", dir=a.tmp) for _ in range(a.threads)]   # a pipe would fill and hang
    procs = [subprocess.Popen(["kissat", "-q", f"--seed={s}", f"--time={int(a.time)}", tmp.name],
                              stdout=outs[s], text=True) for s in range(a.threads)]
    out = None
    try:
        while out is None and any(p.poll() is None for p in procs):
            time.sleep(0.5)
            for p, f in zip(procs, outs):
                if p.poll() in (10, 20):
                    f.seek(0)
                    out = (p.returncode, f.read())
                    break
        if out is None:
            for p, f in zip(procs, outs):
                if p.poll() in (10, 20):
                    f.seek(0)
                    out = (p.returncode, f.read())
                    break
    finally:
        for p in procs:
            if p.poll() is None:
                p.kill()
        os.unlink(tmp.name)
    if out is None:
        print(f"UNKNOWN after {time.time() - t0:.0f}s")
        return
    if out[0] == 20:
        print(f"UNSAT in {time.time() - t0:.1f}s")
        return
    model = set(int(t) for l in out[1].splitlines() if l.startswith("v") for t in l.split()[1:])
    cells = [c for c, l in cell.items() if l in model]
    per, seen = check(cells, rule)
    miss = sum(1 for g in groups.values() if not (set(g) & seen))
    extra = seen - set(allowed)
    print(f"SAT in {time.time() - t0:.1f}s: pop {len(cells)}, {per['kind']} p{per['period']}, "
          f"targets missing {miss}, oriented windows {len(seen & set(real))}/{len(real)}, "
          f"windows outside realizable set {len(extra)}")
    rle = to_rle(cells, W, H, rule)
    print(rle, end="")
    if a.out:
        open(a.out, "w").write(rle)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("classify")
    c.add_argument("--radius", type=int, default=7, help="witness search box half-width")
    c.add_argument("--rule", default="B3/S23")
    f = sub.add_parser("find")
    f.add_argument("--rule", default="B3/S23")
    f.add_argument("width", type=int)
    f.add_argument("height", type=int)
    f.add_argument("--mode", choices=["oriented", "orbit"], default="oriented")
    f.add_argument("--time", type=float, default=300)
    f.add_argument("--threads", type=int, default=10)
    f.add_argument("--maxpop", type=int)
    f.add_argument("--incremental", action="store_true",
                   help="add targets one at a time in one CaDiCaL (SIGTERM to stop); needs --out to keep partial work")
    f.add_argument("--sym", action="store_true",
                   help="force D4 symmetry; with --mode orbit this covers every oriented window too")
    f.add_argument("--out")
    f.add_argument("--tmp", default=tempfile.gettempdir())
    q = sub.add_parser("polish", help="reduce the population of a pattern by patch-wise SAT")
    q.add_argument("start", help="RLE file")
    q.add_argument("--rule", default="B3/S23")
    q.add_argument("--mode", choices=["oriented", "orbit"], default="orbit")
    q.add_argument("--time", type=float, default=600)
    q.add_argument("--pad", type=int, default=1, help="dead margin added around the pattern")
    q.add_argument("--minpatch", type=int, default=5)
    q.add_argument("--maxpatch", type=int, default=9)
    q.add_argument("--budget", type=int, default=5000, help="conflicts per attempt")
    q.add_argument("--seed", type=int, default=1)
    q.add_argument("--whole", action="store_true", help="free the whole box instead of patches")
    q.add_argument("--drop", type=int, default=1, help="--whole: cells to lose per step")
    q.add_argument("--noise", type=float, default=0.0, help="--whole: fraction of seed phases flipped")
    q.add_argument("--out")
    a = p.parse_args()
    if a.cmd == "classify":
        rule, stable = load_rule(a.rule)
        res, wit, proofs = classify(a, stable)
        real = [w for w in range(512) if res[w] == 1]
        orb = orbit_id()
        hist = {}
        for w, k in proofs.items():
            hist[k] = hist.get(k, 0) + 1
        centre = sum(stable)
        print(f"{rule}: centre stable {centre}; realizable {len(real)} ({len({orb[w] for w in real})} D4 classes); "
              f"impossible by proof radius {sorted(hist.items())}; open {[w for w in res if res[w] is None]}")
        json.dump({"rule": rule, "realizable": real, "open": [w for w in res if res[w] is None],
                   "witness": {str(w): wit[w] for w in real}}, open(windows_file(rule), "w"))
    elif a.cmd == "polish":
        polish(a)
    else:
        find(a)


if __name__ == "__main__":
    main()
