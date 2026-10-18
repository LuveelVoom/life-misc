#!/usr/bin/env python3
"""Brute-force crawler search.

A crawler (LifeWiki) is a reaction in which an active pattern collides with an object and
comes back unchanged but moved, while the object is left moved backwards or forwards (and
possibly rephased). Repeated on a line of objects, it travels forever (e.g. the 17c/45
pi heptomino crawler over blinkers, or the 31c/240 Herschel crawler through blocks).

The active patterns are either given (--active RLE, tried in all 8 orientations at every
position within --distance of the object), or every pattern of the cells in a box (--region).
A pattern A is reported if after t <= --max-period generations the board is exactly
    A moved by d   plus   the object (in any phase) moved by e
with A's move d different from 0 and from e (otherwise nothing travels along the objects).
With --allow-debris the reaction may also leave other cells behind (like the pi crawler's
exhaust); A then only has to come back moved, with nothing else within 2 cells of it.
Each find must then repeat on a line of four objects, spaced so that the reaction can repeat:
after each of three cycles A is back (moved) and the objects ahead are untouched.

Line mode (--active with --line DX DY) is for crawlers that only start repeating after a
startup phase (e.g. a bare pi heptomino over blinkers): the pattern is run in front of a line
of objects at each spacing in --spacing (and each phase step between oscillator objects), and
a crawler is reported when the board near the front of the reaction repeats, moved along the
line, twice in a row with the same period.

Examples (a block at 0 0; the last finds the 17c/45 pi crawler in about 90 s):
  python3 crawler_search.py --rule B3/S23-q4i --object '2o$2o!' --active 'o$3o$obo$2bo!' --max-period 60
  python3 crawler_search.py --rule B3/S23-q4i --object '2o$2o!' --region -6 -1 5 4 --max-period 30
  python3 crawler_search.py --object '3o!' --object-pos 13 18 --active '3o$obo$obo!' --positions 9 20 9 30 \\
      --line 0 -1 --spacing 17 17 --max-period 50 --max-gen 300 --window 12

Coordinates: x to the right, y down; the object's top-left cell is at --object-pos (default 0 0).
"""
import argparse, itertools, json, math, os, sys, time
import numpy as np

# Neighbour offsets and their bit in the 8-bit neighbourhood mask (Hensel table order).
NB = {(-1, -1): 0x80, (0, -1): 0x40, (1, -1): 0x20, (-1, 0): 0x10,
      (1, 0): 0x08, (-1, 1): 0x04, (0, 1): 0x02, (1, 1): 0x01}


def rule_table(rule):
    """512-entry next-state table, index = alive * 256 + neighbourhood mask.
    Accepts Hensel/INT notation such as B3/S23-q4i, B2-a/S12, B36/S23."""
    hensel = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hensel.json')))
    birth, surv = [0] * 256, [0] * 256
    mode, i, r = -1, 0, rule.strip()
    while i < len(r):
        ch = r[i]
        if ch == '/':
            i += 1; continue
        if ch in 'Bb':
            mode = 0; i += 1; continue
        if ch in 'Ss':
            mode = 1; i += 1; continue
        if not ch.isdigit() or mode < 0:
            sys.exit(f"Cannot read rule {rule!r} at {ch!r}")
        i += 1
        neg = i < len(r) and r[i] == '-'
        if neg:
            i += 1
        letters = ''
        while i < len(r) and r[i].isalpha() and r[i] not in 'BbSs':
            letters += r[i].lower(); i += 1
        dest = birth if mode == 0 else surv
        for n, letter, masks in hensel:
            if n == ch and (not letters or ((letter not in letters) if neg else (letter in letters))):
                for m in masks:
                    dest[m] = 1
    if birth[0]:
        sys.exit("B0 rules are not supported")
    return np.array(birth + surv, dtype=np.uint8)


def step(a, table):
    """One generation of a batch of boards a[batch, y, x] (uint8 0/1); cells beyond the edge are dead."""
    b, h, w = a.shape
    p = np.zeros((b, h + 2, w + 2), dtype=np.uint16)
    p[:, 1:-1, 1:-1] = a
    idx = a.astype(np.uint16) << 8
    for (dx, dy), bit in NB.items():
        idx |= p[:, 1 + dy:1 + dy + h, 1 + dx:1 + dx + w] * bit
    return table[idx]


def step_live(a, table):
    """step() on the part of the canvas that can change: the box around all live cells of the batch."""
    occ = a.any(axis=0)
    rows, cols = np.nonzero(occ.any(axis=1))[0], np.nonzero(occ.any(axis=0))[0]
    out = np.zeros_like(a)
    if len(rows):
        y0, y1 = max(rows[0] - 1, 0), min(rows[-1] + 2, a.shape[1])
        x0, x1 = max(cols[0] - 1, 0), min(cols[-1] + 2, a.shape[2])
        out[:, y0:y1, x0:x1] = step(a[:, y0:y1, x0:x1], table)
    return out


def run(cells, gens, table):
    """Evolve a set of (x, y) cells on the plane; returns the list of sets for generations 0..gens.
    The board grows only when something comes near its edge."""
    if not cells:
        return [set()] * (gens + 1)
    pad = 8
    x0 = min(x for x, y in cells) - pad; y0 = min(y for x, y in cells) - pad
    w = max(x for x, y in cells) - x0 + pad + 1; h = max(y for x, y in cells) - y0 + pad + 1
    a = np.zeros((1, h, w), dtype=np.uint8)
    for x, y in cells:
        a[0, y - y0, x - x0] = 1
    out = [set(cells)]
    for _ in range(gens):
        if a[0, :2].any() or a[0, -2:].any() or a[0, :, :2].any() or a[0, :, -2:].any():
            a = np.pad(a, ((0, 0), (pad, pad), (pad, pad)))
            x0 -= pad; y0 -= pad
        a = step_live(a, table)
        ys, xs = np.nonzero(a[0])
        out.append({(int(x) + x0, int(y) + y0) for x, y in zip(xs, ys)})
    return out


def move(cells, d):
    return {(x + d[0], y + d[1]) for x, y in cells}


def decode_rle(t):
    cells, x, y, n = set(), 0, 0, 0
    for ch in t:
        if ch.isdigit():
            n = n * 10 + int(ch)
        elif ch in 'ob$.A':
            n = n or 1
            if ch in 'oA':
                cells.update((x + k, y) for k in range(n)); x += n
            elif ch in 'b.':
                x += n
            else:
                y += n; x = 0
            n = 0
        elif ch == '!':
            break
    return cells


def encode_rle(cells):
    """(rle, width, height) of a set of (x, y) cells."""
    if not cells:
        return '!', 0, 0
    x0 = min(x for x, y in cells); y0 = min(y for x, y in cells)
    x1 = max(x for x, y in cells); y1 = max(y for x, y in cells)
    rows = []
    for y in range(y0, y1 + 1):
        row = ''
        for alive, run_ in itertools.groupby(((x, y) in cells for x in range(x0, x1 + 1))):
            n = len(list(run_))
            row += (str(n) if n > 1 else '') + ('o' if alive else 'b')
        rows.append(row.rstrip('0123456789b') if row.endswith('b') else row)  # drop trailing dead run
    return '$'.join(rows) + '!', x1 - x0 + 1, y1 - y0 + 1


def normal(cells):
    """Shape (moved to the origin) and the offset it was moved by."""
    x0 = min(x for x, y in cells); y0 = min(y for x, y in cells)
    return frozenset((x - x0, y - y0) for x, y in cells), (x0, y0)


class ObjectInfo:
    """The object's phases over one period, its period and its own movement per period
    (non-zero for a spaceship, e.g. gliders in a climber)."""

    def __init__(self, cells, table, maxgen):
        ev = run(cells, maxgen, table)
        self.period, self.shift = 0, (0, 0)
        s0, o0 = normal(cells)
        for t in range(1, maxgen + 1):
            if not ev[t]:
                break
            s, o = normal(ev[t])
            if s == s0:
                self.period, self.shift = t, (o[0] - o0[0], o[1] - o0[1])
                break
        if not self.period:
            sys.exit(f'The object does not repeat within {maxgen} generations '
                     '(it must be a still life, oscillator or spaceship in this rule)')
        self.phases = ev[:self.period]           # absolute cells, phase k at generation k

    def at(self, t):
        """The object's free evolution at generation t (t may be negative)."""
        k, n = t % self.period, t // self.period
        return move(self.phases[k], (n * self.shift[0], n * self.shift[1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rule', default='B3/S23')
    ap.add_argument('--object', required=True, help='RLE of the object the pattern crawls over')
    ap.add_argument('--object-pos', type=int, nargs=2, default=[0, 0], metavar=('X', 'Y'))
    ap.add_argument('--active', action='append', metavar='RLE',
                    help='an active pattern to try at every position near the object (repeatable)')
    ap.add_argument('--distance', type=int, default=8,
                    help='--active: try positions up to this many cells from the object (default 8)')
    ap.add_argument('--positions', type=int, nargs=4, metavar=('X', 'Y', 'W', 'H'),
                    help='--active: try only these top-left positions instead (a box)')
    ap.add_argument('--as-given', action='store_true',
                    help='--active: use the patterns only as given (default: all 8 orientations)')
    ap.add_argument('--region', type=int, nargs=4, metavar=('X', 'Y', 'W', 'H'),
                    help='instead of --active: try every pattern of the cells in this box (object cells excluded)')
    ap.add_argument('--max-period', type=int, default=30, help='longest period tried (default 30)')
    ap.add_argument('--max-cells', type=int, default=0,
                    help='only patterns with at most this many live cells (0 = no limit)')
    ap.add_argument('--allow-debris', action='store_true',
                    help='the reaction may leave more than the moved object behind (debris, gliders...), as long '
                         'as it still repeats on a line of objects; slower')
    ap.add_argument('--line', type=int, nargs=2, metavar=('DX', 'DY'),
                    help='--active: line mode. Objects are placed in a line in this direction (e.g. 0 -1 = '
                         'upwards); finds crawlers that only start repeating after a startup phase')
    ap.add_argument('--spacing', type=int, nargs=2, default=[3, 20], metavar=('MIN', 'MAX'),
                    help='--line: object spacings tried, in steps of the direction (default 3 20)')
    ap.add_argument('--max-gen', type=int, default=0,
                    help='--line: generations simulated (default 4 x max-period)')
    ap.add_argument('--window', type=int, default=0,
                    help='--line: cells behind the front of the reaction that must repeat (default 2 x spacing)')
    ap.add_argument('--check-cycles', type=int, default=10,
                    help='cycles a find must repeat for on a line of objects (default 10)')
    ap.add_argument('--batch', type=int, default=2048, help='patterns simulated at once')
    ap.add_argument('--out', default='crawlers.rle')
    args = ap.parse_args()

    table = rule_table(args.rule)
    P = args.max_period
    global CYCLES
    CYCLES = args.check_cycles
    if bool(args.active) == bool(args.region):
        sys.exit('Give either --active (patterns to place) or --region (box to enumerate)')
    obj = move(decode_rle(args.object), args.object_pos)
    if not obj:
        sys.exit('--object has no live cells')
    oi = ObjectInfo(obj, table, P)
    kind = 'still life' if oi.period == 1 else (f'oscillator, period {oi.period}' if oi.shift == (0, 0)
                                                else f'spaceship, ({oi.shift[0]},{oi.shift[1]})c/{oi.period}')

    if args.active:
        # every orientation of every pattern, at every position that does not overlap the object
        shapes = set()
        for rle in args.active:
            cells = decode_rle(rle)
            if not cells:
                sys.exit(f'--active {rle} has no live cells')
            for sx, sy, sw in itertools.product((1, -1), (1, -1), (False, True)):
                if args.as_given and (sx, sy, sw) != (1, 1, False):
                    continue
                shapes.add(normal({((y if sw else x) * sx, (x if sw else y) * sy) for x, y in cells})[0])
        ox0, oy0 = min(x for x, y in obj), min(y for x, y in obj)
        ox1, oy1 = max(x for x, y in obj), max(y for x, y in obj)
        placements = set()
        for shape in shapes:
            w = max(x for x, y in shape) + 1; h = max(y for x, y in shape) + 1
            if args.positions:
                px, py, pw, ph = args.positions
                xr, yr = range(px, px + pw), range(py, py + ph)
            else:
                D = args.distance
                xr, yr = range(ox0 - D - w, ox1 + D + 2), range(oy0 - D - h, oy1 + D + 2)
            for x in xr:
                for y in yr:
                    placed = frozenset(move(shape, (x, y)))
                    if placed.isdisjoint(obj):
                        placements.add(placed)
        placements = sorted(placements, key=sorted)
        total = len(placements)
        span = set().union(*placements)
        print(f'rule {args.rule}, object {kind}, {len(shapes)} pattern orientation(s), {total} placements, '
              f'periods up to {P}')
        what = f'active {" ".join(args.active)}'
        if args.line:
            return line_search(args, table, obj, oi, placements)
    elif args.line:
        sys.exit('--line needs --active')
    else:
        gx, gy, gw, gh = args.region
        span = {(gx, gy), (gx + gw - 1, gy + gh - 1)}
        what = f'region {gx} {gy} {gw} {gh}'

    # Canvas: patterns and object plus a margin of P + 2 (nothing spreads faster than one cell per
    # generation, so nothing can reach the edge within P generations).
    xs_all = [x for x, y in obj | span]; ys_all = [y for x, y in obj | span]
    m = P + 2
    X0, Y0 = min(xs_all) - m, min(ys_all) - m
    W, H = max(xs_all) - X0 + m + 1, max(ys_all) - Y0 + m + 1
    base = np.zeros((H, W), dtype=np.uint8)
    for x, y in obj:
        base[y - Y0, x - X0] = 1

    if args.active:
        def boards():
            for i in range(0, total, args.batch):
                chunk = placements[i:i + args.batch]
                g = np.repeat(base[None], len(chunk), axis=0)
                for k, cells in enumerate(chunk):
                    for x, y in cells:
                        g[k, y - Y0, x - X0] = 1
                yield g
    else:
        free = [(y, x) for y in range(gy, gy + gh) for x in range(gx, gx + gw) if (x, y) not in obj]
        n = len(free)
        total = sum(math.comb(n, c) for c in range(1, args.max_cells + 1)) if args.max_cells else 2 ** n - 1
        print(f'rule {args.rule}, object {kind}, region {gw}x{gh} with {n} free cells, {total} patterns, '
              f'periods up to {P}')
        if total > 50_000_000:
            sys.exit('Too many patterns. Use a smaller region or --max-cells.')
        fy = np.array([y - Y0 for y, x in free]); fx = np.array([x - X0 for y, x in free])

        def boards():
            def fill(bits):
                g = np.repeat(base[None], len(bits), axis=0)
                g[:, fy, fx] = bits
                return g
            if args.max_cells:
                buf = []
                for c in range(1, args.max_cells + 1):
                    for combo in itertools.combinations(range(n), c):
                        v = np.zeros(n, dtype=np.uint8); v[list(combo)] = 1
                        buf.append(v)
                        if len(buf) == args.batch:
                            yield fill(np.array(buf)); buf = []
                if buf:
                    yield fill(np.array(buf))
            else:
                bits = np.arange(n, dtype=np.int64)
                for start in range(1, 2 ** n, args.batch):
                    ids = np.arange(start, min(start + args.batch, 2 ** n), dtype=np.int64)
                    yield fill(((ids[:, None] >> bits) & 1).astype(np.uint8))

    # The object's phases as small arrays, for finding where a phase sits in a board.
    phase_shapes = []
    for k in range(oi.period):
        s, o = normal(oi.phases[k])
        phase_shapes.append((k, s, o, len(s)))
    obj_pops = {len(s) for k, s, o, p in phase_shapes}

    results = []
    tried, t0 = 0, time.time()
    for g in boards():
        B = len(g)
        popA = (g ^ base).reshape(B, -1).sum(axis=1).astype(np.int64)
        idx = np.arange(B)
        cur, prev1, prev2 = g, None, None
        for t in range(1, P + 1):
            prev2, prev1 = prev1, cur
            cur = step_live(cur, table)
            # Edges absorb what reaches them (escaping spaceships). Without --allow-debris such a
            # reaction is not clean, so it is dropped.
            edge = cur[:, :2, :].any((1, 2)) | cur[:, -2:, :].any((1, 2)) | cur[:, :, :2].any((1, 2)) | cur[:, :, -2:].any((1, 2))
            cur[:, :2, :] = 0; cur[:, -2:, :] = 0; cur[:, :, :2] = 0; cur[:, :, -2:] = 0
            pop = cur.reshape(len(idx), -1).sum(axis=1)
            keep = (pop > 0) & (~edge if not args.allow_debris else True)
            # a board that is now still or period 2 stays that way: it cannot become a crawler
            flat = cur.reshape(len(idx), -1)
            settled = (flat == prev1.reshape(len(idx), -1)).all(axis=1)
            if prev2 is not None:
                settled |= (flat == prev2.reshape(len(idx), -1)).all(axis=1)
            cand = pop > popA[idx] if args.allow_debris else np.zeros(len(idx), dtype=bool)
            for p in obj_pops:
                cand |= pop == popA[idx] + p
            for i in np.nonzero(cand & keep)[0]:
                j = idx[i]
                res = decompose(cur[i], g[j] ^ base, t, phase_shapes, oi, X0, Y0, table, args)
                if res and repeats(res, oi, table):
                    keep[i] = False
                    record(res, oi, results)
            keep &= ~settled
            idx, cur, prev1 = idx[keep], cur[keep], prev1[keep]
            if not len(idx):
                break
        tried += B
        if time.time() - t0 > 10:
            print(f'  {tried} of {total} starting patterns tried, {len(results)} found so far', flush=True)
            t0 = time.time()

    print(f'{tried} starting patterns tried: {len(results)} result(s)')
    if results:
        with open(args.out, 'w') as f:
            f.write(f'#C crawler_search.py: {len(results)} result(s) in {args.rule}; object {args.object} at '
                    f'{args.object_pos[0]},{args.object_pos[1]}, {what}\n'
                    f'#C Each row: the pattern and a line of {CYCLES + 1} objects, spaced so that the reaction repeats '
                    f'(checked for {CYCLES} cycles).\n')
            for k, r in enumerate(results, 1):
                f.write(f'#C {k}: {r["desc"]}\n')
            cells, yoff = set(), 0
            for r in results:
                c = r['stream']
                y0 = min(y for x, y in c); x0 = min(x for x, y in c)
                cells |= {(x - x0, y - y0 + yoff) for x, y in c}
                yoff += max(y for x, y in c) - y0 + 1 + 10
            body, w, h = encode_rle(cells)
            f.write(f'x = {w}, y = {h}, rule = {args.rule}\n')
            f.write('\n'.join(body[i:i + 70] for i in range(0, len(body), 70)) + '\n')
        print(f'Written to {args.out}')


def decompose(board, d0, t, phase_shapes, oi, X0, Y0, table, args):
    """Is board (generation t) = the starting pattern moved by d + an object phase moved by e?"""
    ys, xs = np.nonzero(d0)
    A = {(int(x) + X0, int(y) + Y0) for x, y in zip(xs, ys)}
    shapeA, offA = normal(A)
    ys, xs = np.nonzero(board)
    S = {(int(x) + X0, int(y) + Y0) for x, y in zip(xs, ys)}
    separate = None   # the pattern and the object evolved apart (computed when first needed)
    if args.allow_debris:
        return decompose_debris(A, S, t, oi, table)
    for k, shapeO, offO, popO in phase_shapes:
        if len(S) != len(A) + popO:
            continue
        # where can this phase sit? anchor its first cell on each live cell
        ax, ay = min(shapeO)
        for sx, sy in S:
            q = (sx - ax, sy - ay)
            placed = move(shapeO, q)
            if not placed <= S:
                continue
            rest = S - placed
            shapeR, offR = normal(rest)
            if shapeR != shapeA:
                continue
            d = (offR[0] - offA[0], offR[1] - offA[1])
            # e: where phase k is, relative to where the object's free evolution has phase k at
            # the latest generation <= t; actual: how far it is from its starting position.
            tk = t - (t - k) % oi.period
            e = (q[0] - normal(oi.at(tk))[1][0], q[1] - normal(oi.at(tk))[1][1])
            start = normal(oi.phases[k])[1]
            actual = (q[0] - start[0], q[1] - start[1])
            if d == (0, 0) or d == actual or max(abs(d[0]), abs(d[1])) > t:
                continue   # did not move, moved along with the object, or faster than light (cells swapped roles)
            if separate is None:
                separate = run(A, t, table)[t] | oi.at(t)
            if S == separate:
                return None    # no reaction yet (e.g. a spaceship still on its way to the object)
            return {'A': A, 'd': d, 'e': e, 'k': k, 't': t, 'objfinal': placed}
    return None


def decompose_debris(A, S, t, oi, table):
    """--allow-debris: the pattern is back, moved by d, with nothing else within 2 cells of it;
    everything else on the board counts as the remains of the object."""
    shapeA, offA = normal(A)
    ax, ay = min(shapeA)
    separate = None
    for sx, sy in S:
        q = (sx - ax, sy - ay)
        placed = move(shapeA, q)
        if (sx, sy) != min(placed) or not placed <= S:
            continue
        d = (q[0] - offA[0], q[1] - offA[1])
        if d == (0, 0) or max(abs(d[0]), abs(d[1])) > t:
            continue
        near = {(x + i, y + j) for x, y in placed for i in range(-2, 3) for j in range(-2, 3)}
        if (S & near) != placed:
            continue           # something else is touching it
        if separate is None:
            separate = run(A, t, table)[t] | oi.at(t)
        if S == separate:
            return None        # no reaction yet
        return {'A': A, 'd': d, 'e': None, 'k': 0, 't': t, 'objfinal': S - placed}
    return None


# ---------------------------------------------------------------- line mode
def line_search(args, table, obj, oi, placements):
    """Line mode: the pattern is placed in front of a line of objects (spacing n along the
    direction, object j in phase j*phi mod period) and simulated for max-gen generations. A
    crawler is reported when the board near the front of the reaction (from --window cells behind
    the frontmost changed cell to one spacing ahead) repeats moved along the line, twice in a row
    with the same period: the reaction may take any number of generations to settle first."""
    dx, dy = args.line
    if (dx, dy) == (0, 0) or max(abs(dx), abs(dy)) != 1:
        sys.exit('--line DX DY: a direction with entries -1, 0 or 1, e.g. 0 -1 (up) or 1 1 (down-right)')
    P = args.max_period
    G = args.max_gen or 4 * P
    q = dx * dx + dy * dy                      # projection units per cell along the direction
    proj = lambda c: c[0] * dx + c[1] * dy
    results, seen = [], set()
    nmin, nmax = args.spacing
    combos = [(n, phi) for n in range(nmin, nmax + 1) for phi in range(oi.period)]
    print(f'line mode: direction ({dx},{dy}), spacings {nmin}-{nmax}, {oi.period} phase step(s), '
          f'{len(placements)} placements, {G} generations, periods up to {P}')
    t_start = time.time()
    for ci, (n, phi) in enumerate(combos):
        s = (n * dx, n * dy)
        depth = (args.window or 2 * n) * q
        # objects needed: the front moves at most 1 cell/gen relative to the objects, which may move too
        K = (2 * G + 2 * n) // n + 3
        objs = [move(oi.at(j * phi), (j * s[0], j * s[1])) for j in range(K)]

        def line_at(t):
            return set().union(*[move(oi.at(t + j * phi), (j * s[0], j * s[1])) for j in range(K)])
        line = set().union(*objs)
        if any(len(o & p) for i, o in enumerate(objs) for p in objs[i + 1:]):
            continue
        # objects must not interact with each other
        free = run(line, oi.period, table)
        if any(free[t] != line_at(t) for t in range(1, oi.period + 1)):
            continue
        p1 = proj(next(iter(objs[1])))
        pl = [A for A in placements if A.isdisjoint(line) and max(proj(c) for c in A) < max(proj(c) for c in objs[1])]
        if not pl:
            continue
        # canvas: the line (wherever it is during the run) and placements, margin P + 2.
        # Its edges absorb whatever reaches them (escaping spaceships, debris).
        Ft = [line_at(t) for t in range(G + 2)]
        cells = set().union(*Ft) | set().union(*pl)
        m = P + 2
        X0 = min(x for x, y in cells) - m; Y0 = min(y for x, y in cells) - m
        W = max(x for x, y in cells) - X0 + m + 1; H = max(y for x, y in cells) - Y0 + m + 1
        F = np.zeros((G + 2, H, W), dtype=np.uint8)      # the line on its own, every generation
        for t in range(G + 2):
            for x, y in Ft[t]:
                F[t, y - Y0, x - X0] = 1
        yy, xx = np.mgrid[0:H, 0:W]
        Pm = ((xx + X0) * dx + (yy + Y0) * dy).astype(np.int64)
        unit = n * q                                     # projection change per spacing
        bs = min(args.batch, 256)
        for b0 in range(0, len(pl), bs):
            chunk = pl[b0:b0 + bs]
            B = len(chunk)
            # D = difference from the free line (F); the board is F ^ D. Only the box around D is
            # stepped and examined. Edges absorb whatever reaches them (escaping spaceships).
            D = np.zeros((B, H, W), dtype=np.uint8)
            for k, A in enumerate(chunk):
                for x, y in A:
                    D[k, y - Y0, x - X0] = 1
            ids = np.arange(B)                     # boards still being simulated
            rejected = set()                       # (board, period, move) already checked and refused
            Hh = np.zeros((B, G + 1), dtype=np.uint64)
            Ox = np.zeros((B, G + 1), dtype=np.int64); Oy = np.zeros((B, G + 1), dtype=np.int64)
            occ = D.any(axis=0)
            rows, cols = np.nonzero(occ.any(axis=1))[0], np.nonzero(occ.any(axis=0))[0]
            for t in range(G + 1):
                if t:
                    y0, y1 = max(rows[0] - 2, 0), min(rows[-1] + 3, H)
                    x0, x1 = max(cols[0] - 2, 0), min(cols[-1] + 3, W)
                    stepped = step(F[t - 1][y0:y1, x0:x1] ^ D[:, y0:y1, x0:x1], table)
                    # the box's own edge cells were computed with dead cells outside: keep the inside
                    iy0, iy1 = (1 if y0 > 0 else 0), (y1 - y0 - 1 if y1 < H else y1 - y0)
                    ix0, ix1 = (1 if x0 > 0 else 0), (x1 - x0 - 1 if x1 < W else x1 - x0)
                    D[:, y0:y1, x0:x1] = 0
                    D[:, y0 + iy0:y0 + iy1, x0 + ix0:x0 + ix1] = (stepped[:, iy0:iy1, ix0:ix1]
                                                                 ^ F[t][y0 + iy0:y0 + iy1, x0 + ix0:x0 + ix1])
                    D[:, :2, :] = 0; D[:, -2:, :] = 0; D[:, :, :2] = 0; D[:, :, -2:] = 0
                    occ = D[:, y0:y1, x0:x1].any(axis=0)
                    rr, cc = np.nonzero(occ.any(axis=1))[0], np.nonzero(occ.any(axis=0))[0]
                    if not len(rr):
                        break
                    rows, cols = rr + y0, cc + x0
                # examine the changes plus one spacing ahead
                r0, r1 = max(rows[0] - n - 2, 0), min(rows[-1] + n + 3, H)
                c0, c1 = max(cols[0] - n - 2, 0), min(cols[-1] + n + 3, W)
                Ds = D[:, r0:r1, c0:c1]
                gs = Ds ^ F[t][r0:r1, c0:c1]
                Ps = Pm[r0:r1, c0:c1]
                Bn = len(ids)
                alive = Ds.reshape(Bn, -1).any(axis=1)
                front = np.where(Ds.astype(bool), Ps[None], np.iinfo(np.int64).min // 2).reshape(Bn, -1).max(axis=1)
                win = gs.astype(bool) & (Ps[None] >= (front - depth)[:, None, None]) & (Ps[None] <= (front + unit)[:, None, None])
                bi, wy, wx = np.nonzero(win)
                # position-independent hash of the window (the line itself may be moving)
                big = np.iinfo(np.int64).max
                mx = np.full(Bn, big, dtype=np.int64); my = np.full(Bn, big, dtype=np.int64)
                np.minimum.at(mx, bi, wx); np.minimum.at(my, bi, wy)
                xr = (wx - mx[bi]).astype(np.int64); yr = (wy - my[bi]).astype(np.int64)
                v = (xr * 0x9E3779B1 + yr * 0x85EBCA77 + 0x165667B1).astype(np.uint64)
                v ^= v >> np.uint64(15); v *= np.uint64(0x2C1B3C6D297A1A2B); v ^= v >> np.uint64(29)
                h = np.zeros(Bn, dtype=np.uint64)
                np.add.at(h, bi, v)
                Hh[ids, t] = h
                Ox[ids, t] = mx + c0; Oy[ids, t] = my + r0
                done = np.zeros(Bn, dtype=bool)
                Tm = min(P, t // 2)
                cols_ = np.arange(t - 2 * Tm, t + 1) if Tm else np.arange(t, t + 1)
                hx, hy, hh = Ox[ids][:, cols_], Oy[ids][:, cols_], Hh[ids][:, cols_]   # generations t-2Tm..t
                base_ = t - 2 * Tm
                for T in range(1, Tm + 1):
                    a0, a1, a2 = t - base_, t - T - base_, t - 2 * T - base_
                    ddx, ddy = hx[:, a0] - hx[:, a1], hy[:, a0] - hy[:, a1]
                    hit = (alive & ~done & (hh[:, a0] == hh[:, a1]) & (hh[:, a1] == hh[:, a2])
                           & (ddx == hx[:, a1] - hx[:, a2]) & (ddy == hy[:, a1] - hy[:, a2])
                           & ((ddx != 0) | (ddy != 0)))
                    if not hit.any():
                        continue
                    for i in np.nonzero(hit)[0]:
                        key = (int(ids[i]), T, int(ddx[i]), int(ddy[i]))
                        if key in rejected or passes(T, (int(ddx[i]), int(ddy[i])), s, phi, oi) is None:
                            continue
                        if line_record(chunk[ids[i]], objs, n, phi, s, int(T), t - 2 * T, depth, unit,
                                       proj, table, results, seen, oi):
                            done[i] = True
                        else:
                            rejected.add(key)
                keep = alive & ~done
                if not keep.all():
                    if not keep.any():
                        break
                    ids, D = ids[keep], D[keep]
                    occ = D.any(axis=0)
                    rows, cols = np.nonzero(occ.any(axis=1))[0], np.nonzero(occ.any(axis=0))[0]
        print(f'  spacing {n}, phase step {phi} done ({ci + 1}/{len(combos)}, {time.time() - t_start:.0f} s), '
              f'{len(results)} found so far', flush=True)

    print(f'{len(results)} result(s)')
    if results:
        with open(args.out, 'w') as f:
            f.write(f'#C crawler_search.py line mode: {len(results)} result(s) in {args.rule}; object {args.object}, '
                    f'active {" ".join(args.active)}, direction ({dx},{dy})\n'
                    f'#C Each row: the starting pattern and the line of objects it crawls along (checked for {CYCLES} cycles).\n')
            for k, r in enumerate(results, 1):
                f.write(f'#C {k}: {r["desc"]}\n')
            cells, yoff = set(), 0
            for r in results:
                c = r['stream']
                y0 = min(y for x, y in c); x0 = min(x for x, y in c)
                cells |= {(x - x0, y - y0 + yoff) for x, y in c}
                yoff += max(y for x, y in c) - y0 + 1 + 10
            body, w, h = encode_rle(cells)
            f.write(f'x = {w}, y = {h}, rule = {args.rule}\n')
            f.write('\n'.join(body[i:i + 70] for i in range(0, len(body), 70)) + '\n')
        print(f'Written to {args.out}')


def passes(T, d, s, phi, oi):
    """Objects passed per period T if the reaction moves d (object j at t+T is object j-k at t,
    moved by d: T + k*phi is a whole number of object periods and d - k*s the objects' own
    movement), or None."""
    return next((k for k in range(1, 2 * T + 1) if (T + k * phi) % oi.period == 0 and
                 (d[0] - k * s[0], d[1] - k * s[1]) == ((T + k * phi) // oi.period * oi.shift[0],
                                                        (T + k * phi) // oi.period * oi.shift[1])), None)


def line_record(A, objs, n, phi, s, T, t0, depth, unit, proj, table, results, seen, oi):
    """Confirm a line-mode find exactly (no hashing, unbounded board) and record it. Returns
    False if it is not a crawler after all."""
    # the exact check runs longer than the search, so it needs a longer line
    K2 = (2 * (t0 + CYCLES * T) + 2 * n) // n + 3
    objs = objs + [move(oi.at(j * phi), (j * s[0], j * s[1])) for j in range(len(objs), K2)]
    ev = run(set(A) | set().union(*objs), t0 + CYCLES * T, table)

    def window(t, dep=depth):
        free = set().union(*[move(oi.at(t + j * phi), (j * s[0], j * s[1])) for j in range(len(objs))])
        D = ev[t] ^ free
        if not D:
            return None, None
        front = max(proj(c) for c in D)
        cells = {c for c in ev[t] if front - dep <= proj(c) <= front + unit}
        return normal(cells) if cells else (None, None)
    wins = [window(t0 + i * T) for i in range(CYCLES + 1)]
    w0, o0 = wins[0]
    if w0 is None or wins[1][0] is None:
        return False
    d = (wins[1][1][0] - o0[0], wins[1][1][1] - o0[1])
    if d == (0, 0) or any(w != w0 or o is None or (o[0] - o0[0], o[1] - o0[1]) != (i * d[0], i * d[1])
                          for i, (w, o) in enumerate(wins)):
        return False          # does not keep repeating for CYCLES cycles
    # how many objects does it pass per period? object j at t+T must be object j-k at t, moved by d:
    # T + k*phi must be a whole number of object periods, and d - k*s the objects' own movement
    k = passes(T, d, s, phi, oi)
    if k is None:
        return False          # repeats, but does not get through any objects (e.g. travels beside them)
    # smallest period of the front itself (one spacing behind to one ahead; the wake further
    # behind may repeat more slowly)
    for T2 in range(1, T):
        if T % T2 or (k * T2) % T:
            continue
        d2 = (d[0] * T2 // T, d[1] * T2 // T)
        w_a, o_a = window(t0, unit); w_b, o_b = window(t0 + T2, unit)
        if w_b == w_a and o_b is not None and (o_b[0] - o_a[0], o_b[1] - o_a[1]) == d2:
            T, k, d = T2, k * T2 // T, d2
            break
    # Fake crawler: the repeating front is just something (e.g. a spaceship) that moves the same way
    # without the objects. Real crawlers may still emit spaceships; those are left behind as debris.
    free0 = set().union(*[move(oi.at(t0 + j * phi), (j * s[0], j * s[1])) for j in range(len(objs))])
    D0 = ev[t0] ^ free0
    front = max(proj(c) for c in D0)
    part = {c for c in D0 if c in ev[t0] and proj(c) >= front - depth}
    if part:
        shape, off = normal(part)
        alone = run(part, T, table)[T]
        if alone and normal(alone) == (shape, (off[0] + d[0], off[1] + d[1])):
            return False
    # the same reaction found from another start (or another phase of it) counts once
    key = (frozenset(window(t0 + i, unit)[0] for i in range(T)), T, k)
    if key in seen:
        return True
    seen.add(key)
    # crawler or fuse: are the objects it passed still there behind it (moved or not)?
    last = ev[t0 + CYCLES * T]
    free_last = set().union(*[move(oi.at(t0 + CYCLES * T + j * phi), (j * s[0], j * s[1])) for j in range(len(objs))])
    behind = {c for c in last if proj(c) < max(proj(c) for c in last ^ free_last) - depth}
    copies = count_objects(behind, oi)
    what = ('crawler' if copies >= k * (CYCLES - 1) else
            'fuse' if copies == 0 else f'crawler or fuse ({copies} objects left behind)')
    speed = f'({abs(d[0])},{abs(d[1])})c/{T}' if d[0] and d[1] else f'{max(abs(d[0]), abs(d[1]))}c/{T}'
    desc = (f'{what}; period {T}, passes {k} object(s) per period, moves ({d[0]},{d[1]}) = {speed}'
            + f'; objects every {n} cells' + (f', each {phi} generation(s) on from the last' if oi.period > 1 else '')
            + f'; repeating from generation {t0}; pattern {encode_rle(set(A))[0]}')
    shown = set(A) | set().union(*objs[:max(10, (t0 + CYCLES * T) * k // T + 3)])
    results.append({'desc': desc, 'stream': shown})
    print(f'  found: {desc}', flush=True)
    return True


def count_objects(cells, oi):
    """Isolated copies of the object (any phase) among cells."""
    n = 0
    for k in range(oi.period):
        shape = normal(oi.phases[k])[0]
        ax, ay = min(shape)
        for sx, sy in cells:
            placed = move(shape, (sx - ax, sy - ay))
            if (sx, sy) == min(placed) and placed <= cells:
                near = {(x + a, y + b) for x, y in placed for a in range(-1, 2) for b in range(-1, 2)}
                if (cells & near) == placed:
                    n += 1
    return n // oi.period if oi.period > 1 and all(len(normal(p)[0]) for p in oi.phases) and False else n


CYCLES = 3   # cycles checked on a line of CYCLES + 1 objects


def repeats(res, oi, table):
    """Does the reaction repeat on a line of objects? Object j must reach the pattern's position
    + j*d in phase 0 exactly when the pattern arrives (generation j*t), so it starts as the
    object's free evolution at generation -j*t, moved by j*d. At the end of each of CYCLES
    cycles, the pattern must be back (moved j*d, nothing else within 2 cells) and every object
    still ahead untouched. What happens behind (remains, debris) is not checked further."""
    A, d, t = res['A'], res['d'], res['t']
    line = set(A)
    for j in range(CYCLES + 1):
        line |= move(oi.at(-j * t), (j * d[0], j * d[1]))
    res['line'] = line
    ev = run(line, CYCLES * t, table)
    for j in range(1, CYCLES + 1):
        S = ev[j * t]
        placed = move(A, (j * d[0], j * d[1]))
        near = {(x + a, y + b) for x, y in placed for a in range(-2, 3) for b in range(-2, 3)}
        if (S & near) != placed:
            return False
        for m in range(j + 1, CYCLES + 1):
            if not move(oi.at((j - m) * t), (m * d[0], m * d[1])) <= S:
                return False
    return True


def find_object(cells, oi):
    """Movement of a copy of the object (any phase) among cells, relative to the same phase at
    the start, or None."""
    for k in range(oi.period):
        shape, start = normal(oi.phases[k])
        ax, ay = min(shape)
        for sx, sy in cells:
            q = (sx - ax, sy - ay)
            if move(shape, q) <= cells:
                return (q[0] - start[0], q[1] - start[1])
    return None


def record(res, oi, results):
    A, d, t, e = res['A'], res['d'], res['t'], res['e']
    if e is None:   # --allow-debris: is the object still there (moved), or used up?
        moved = find_object(res['objfinal'], oi)
        extra = len(res['objfinal']) - len(oi.phases[0])
        what = (f'crawler or fuse: a copy of the object is among the remains, ({moved[0]},{moved[1]}) from where it '
                f'was, plus {extra} other cells' if moved
                else f'fuse (the object is used up; {len(res["objfinal"])} cells left behind)')
    else:
        what = f'crawler, object moved ({e[0]},{e[1]})' + (f' and rephased to phase {res["k"]}' if res['k'] else '')
    speed = f'({abs(d[0])},{abs(d[1])})c/{t}' if d[0] and d[1] else f'{max(abs(d[0]), abs(d[1]))}c/{t}'
    desc = (f'{what}; period {t}, pattern moves ({d[0]},{d[1]}) = {speed}'
            + f'; pattern {encode_rle(A)[0]}')
    results.append({'desc': desc, 'stream': res['line']})
    print(f'  found: {desc}', flush=True)


if __name__ == '__main__':
    main()
