#!/usr/bin/env python3
"""Pull the rare and easy-to-miss objects of one rule out of a rulesearch soup census + probe.

Runs `rulesearch CFG probe RULE N` once (soups, then catalyst / wick / fossil probes) with a discoveries log,
and lays the interesting finds out as one RLE pattern collection in the rule (one band of rows per group, rarest
first; a #C line per object gives what it is and its top-left position, x right, y down):

  ships         seen in at most --rare of the soups, and every oblique ship
  oscillators   period >= --min-period, seen in at most --rare of the soups
  soups         puffers, rakes, guns and never-settling active regions (the RLE is the soup)
  catalysts     ships fired at still lifes that come back: eaters, reflectors, converters, ...
  signals       signals and 1D CA behaviour along wicks (fuses only with --fuses)
  fossils       large rare objects left in the ash
  replicators   patterns that make exact copies of themselves, most controllable first

Ships and oscillators are re-simulated with casim (an independent simulator) when it can be imported; it is optional. Objects are
spaced --gap cells apart, but ships fly off and soups grow when the collection is run, so look before you run.

Example:
  python3 rarefind.py B2-a3/S23 -n 20000 -o rare.rle
  python3 rarefind.py B3/S23-a5 --rare 0.001 --set "soup 20 20 0.4"
"""
import argparse, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CASIM = os.path.join(HERE, '..', '..', 'casim')


def parse_log(text):
    """Log entries: one or more '#' lines, then an RLE ending in '!'."""
    entries, head, body = [], [], []
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith('#') and not body:
            head.append(line)
        elif head:
            body.append(line)
            if line.rstrip().endswith('!'):
                entries.append((head, '\n'.join(body)))
                head, body = [], []
    return entries


def classify(head, soups):
    """-> (section, sort key, description, count) or None."""
    h = head[0]
    m = re.match(r'# (?:NEW BEHAVIOUR|new shape) (.*), (\d+) cells, (\d+) in (\d+) soups', h)
    if m:
        desc, cells, cnt = m.group(1), int(m.group(2)), int(m.group(3))
        sec = 'ships' if 'ship' in desc else 'oscillators'
        return sec, cnt, f'{desc}, {cells} cells, in {cnt} of {soups} soups', cnt
    m = re.match(r'# (?:NEW SOUP BEHAVIOUR|soup behaviour) (.*), (\d+) in (\d+) soups', h)
    if m:
        return 'soups', int(m.group(2)), f'{m.group(1)}; in {m.group(2)} of {soups} soups; {head[1][2:]}', int(m.group(2))
    m = re.match(r'# (?:NEW )?catalyst: (.*), (\d+) setups', h)
    if m:
        return 'catalysts', -int(m.group(2)), f'{m.group(1)} ({m.group(2)} setups)', None
    m = re.match(r'# (?:NEW )?wick: (.*), (\d+) setups', h)
    if m:
        kind = 'fuses' if m.group(1).startswith('fuse') else 'signals'
        return kind, -int(m.group(2)), f'{m.group(1)} ({m.group(2)} setups)', None
    m = re.match(r'# (?:NEW )?replicator: (.*), (\d+) setups', h)
    if m:
        f = float(re.search(r'controllability ([\d.]+)', m.group(1)).group(1))
        return 'replicators', -f, f'replicator: {m.group(1)} (found in {m.group(2)} seeds/soups)', None
    m = re.match(r'# (?:NEW )?fossil: (.*), (\d+) in (\d+) soups', h)
    if m:
        return 'fossils', int(m.group(2)), f'{m.group(1)}; in {m.group(2)} of {soups} soups', int(m.group(2))
    return None


def rle_cells(rle):
    """Live cells (x, y) of a 2-state RLE; header and # lines skipped."""
    body = ''.join(l.strip() for l in rle.splitlines() if l.strip() and not l.startswith(('#', 'x ')))
    cells, x, y, n = [], 0, 0, ''
    for ch in body:
        if ch.isdigit():
            n += ch
            continue
        k = int(n) if n else 1
        n = ''
        if ch == '!':
            break
        if ch == '$':
            y += k; x = 0
        elif ch in 'bB.':
            x += k
        else:
            cells += [(x + i, y) for i in range(k)]
            x += k
    return cells


def rle_body(cells):
    """RLE body of a set of cells (top-left at 0,0), wrapped at 70 characters."""
    rows = {}
    for x, y in cells:
        rows.setdefault(y, set()).add(x)
    toks, py = [], 0
    for y in sorted(rows):
        d = y - py
        if d:
            toks.append(f'{d if d > 1 else ""}$')
        py, xs, x, i = y, sorted(rows[y]), 0, 0
        while i < len(xs):
            j = i
            while j + 1 < len(xs) and xs[j + 1] == xs[j] + 1:
                j += 1
            if xs[i] > x:
                toks.append(f'{xs[i] - x if xs[i] - x > 1 else ""}b')
            toks.append(f'{j - i + 1 if j > i else ""}o')
            x, i = xs[j] + 1, j + 1
    toks.append('!')
    lines, cur = [], ''
    for t in toks:
        if len(cur) + len(t) > 70:
            lines.append(cur); cur = ''
        cur += t
    return '\n'.join(lines + [cur])


def has_b0(rule):
    """True if the rule has birth on 0 neighbours (B/S, S/B or Golly's digits-only S/B form, or MAP bit 0)."""
    r = rule.strip().split(':')[0]
    if r[:3].upper() == 'MAP':
        b64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
        return len(r) > 3 and r[3] in b64 and bool(b64.index(r[3]) >> 5 & 1)   # table[0] = top bit of the first char
    parts = r.split('/')
    for q in parts:
        if q[:1] in 'Bb':
            return q[1:2] == '0'
    return len(parts) == 2 and parts[1][:1] == '0'   # digits only: survival/birth


def casim_check(C, rle, desc):
    """Re-simulate a ship / oscillator: True, False, or None when the description has no period."""
    m = re.search(r'\(p(\d+)\)|^p(\d+) |c/(\d+) oblique', desc)
    if not m:
        return None
    P = int(next(g for g in m.groups() if g))
    p = C.Pattern.from_rle(rle)
    s0 = sorted(p.coords())
    q = p.copy().step(P)
    s1 = sorted(q.coords())
    if len(s0) != len(s1) or not s0:
        return False
    dx, dy = s1[0][0] - s0[0][0], s1[0][1] - s0[0][1]
    same = [(x + dx, y + dy) for x, y in s0] == s1
    return same and ((dx, dy) != (0, 0)) == ('ship' in desc)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rule', help='rule in Golly notation (Life-like, INT or MAP; B0 rules are rejected)')
    ap.add_argument('-n', '--soups', type=int, default=10000, help='soups to run (default 10000)')
    ap.add_argument('--rare', type=float, default=0.005,
                    help='ships/oscillators count as rare when seen in at most this fraction of soups (default 0.005)')
    ap.add_argument('--min-period', type=int, default=3, help='smallest oscillator period kept (default 3)')
    ap.add_argument('--fuses', action='store_true', help='also list fuses along wicks')
    ap.add_argument('--seed', type=int, default=1, help='soup seed (default 1)')
    ap.add_argument('--set', action='append', default=[], metavar='"KEY VALUE"',
                    help='extra rulesearch config line, e.g. "soup 20 20 0.4" or "probe_stills 40" (repeatable)')
    ap.add_argument('-j', '--threads', type=int, help='threads (default: all cores)')
    ap.add_argument('-o', '--out', help='output RLE file (default rare_<rule>.rle)')
    ap.add_argument('--gap', type=int, default=24, help='empty cells between objects (default 24)')
    ap.add_argument('--width', type=int, default=600, help='start a new row past this width (default 600)')
    ap.add_argument('--rulesearch', default=os.path.join(HERE, 'rulesearch'), help='path to the rulesearch binary')
    ap.add_argument('--no-recheck', action='store_true', help='skip the casim re-simulation of ships and oscillators')
    a = ap.parse_args()

    if not os.path.exists(a.rulesearch):
        sys.exit(f'{a.rulesearch} not found: run make in {HERE} first')
    C = None
    if not a.no_recheck:
        try:
            sys.path.insert(0, CASIM)
            import casim as C
        except ImportError:
            print('casim not found; ships and oscillators will not be re-checked', file=sys.stderr)
    if has_b0(a.rule):
        sys.exit('B0 rules are not supported (rulesearch forces B0 off)')
    # rulesearch refuses a non-isotropic rule under symmetry isotropic: then use symmetry none
    sym = None
    with tempfile.TemporaryDirectory() as tmp:
        for sy in ('isotropic', 'none'):
            cfg = os.path.join(tmp, 'sym.cfg')
            open(cfg, 'w').write(f'rule {a.rule}\nsymmetry {sy}\nlog none\narchive none\n')
            r = subprocess.run([a.rulesearch, cfg, 'check'], capture_output=True, text=True)
            if r.returncode == 0:
                sym = sy
                break
            if 'not symmetric' not in r.stdout + r.stderr:
                sys.exit(f'rulesearch cannot use rule {a.rule}:\n{(r.stdout + r.stderr)[-1000:]}')
    if sym is None:
        sys.exit(f'rulesearch cannot use rule {a.rule}')
    out = a.out or 'rare_' + re.sub(r'[^A-Za-z0-9-]+', '_', a.rule) + '.rle'

    with tempfile.TemporaryDirectory() as tmp:
        cfg, log = os.path.join(tmp, 'rare.cfg'), os.path.join(tmp, 'log.txt')
        with open(cfg, 'w') as f:
            f.write(f'rule {a.rule}\nsymmetry {sym}\nlog {log}\narchive none\n' + ''.join(s + '\n' for s in a.set))
        env = dict(os.environ)
        if a.threads:
            env['OMP_NUM_THREADS'] = str(a.threads)
        print(f'running {a.soups} soups + probes in {a.rule} ...', file=sys.stderr)
        r = subprocess.run([a.rulesearch, cfg, 'probe', a.rule, str(a.soups), str(a.seed)],
                           capture_output=True, text=True, env=env)
        if r.returncode != 0:
            sys.exit(f'rulesearch failed:\n{r.stdout[-2000:]}{r.stderr[-2000:]}')
        census = next((l.strip() for l in r.stdout.splitlines() if l.startswith('census')), '')
        text = open(log).read() if os.path.exists(log) else ''

    sections = {k: [] for k in ('replicators', 'ships', 'oscillators', 'soups', 'catalysts', 'signals', 'fuses', 'fossils')}
    lim = a.rare * a.soups
    for head, rle in parse_log(text):
        c = classify(head, a.soups)
        if not c:
            continue
        sec, key, desc, cnt = c
        if sec == 'ships' and cnt > lim and 'oblique' not in desc:
            continue
        if sec == 'oscillators':
            P = int(re.match(r'p(\d+)', desc).group(1))
            if P < a.min_period or cnt > lim:
                continue
        if sec == 'fuses' and not a.fuses:
            continue
        if sec in ('ships', 'oscillators') and C is not None:
            ok = casim_check(C, rle, desc)
            if ok is False:
                desc += '  [casim check FAILED]'
        sections[sec].append((key, desc, rle))

    titles = {'replicators': 'Replicators (most controllable first; they fill space when run)', 'ships': 'Rare ships (and all oblique ships)', 'oscillators': f'Rare oscillators (p{a.min_period}+)',
              'soups': 'Soups with puffers, rakes, guns or never-settling regions (RLE = the soup)',
              'catalysts': 'Catalysts (ship + still life, still life restored)',
              'signals': 'Signals and 1D CA along wicks', 'fuses': 'Fuses along wicks', 'fossils': 'Rare large objects in the ash'}
    # layout: each section starts a new band; objects left to right, wrapping past --width
    order = ['ships', 'oscillators', 'catalysts', 'signals', 'fuses', 'fossils', 'replicators', 'soups']   # growing things last
    comments, live, X, Y, rowh, maxw = [], [], 0, 0, 0, 0
    for sec in order:
        items = [it for it in sorted(sections[sec], key=lambda x: x[0]) if rle_cells(it[2])]
        if not items:
            continue
        if X:
            Y += rowh + 2 * a.gap; X = rowh = 0
        comments.append(f'===== {titles[sec]}: {len(items)} =====')
        for key, desc, rle in items:
            c = rle_cells(rle)
            if not c:
                continue
            x0, y0 = min(x for x, y in c), min(y for x, y in c)
            w, h = max(x for x, y in c) - x0 + 1, max(y for x, y in c) - y0 + 1
            if X and X + w > a.width:
                Y += rowh + a.gap; X = rowh = 0
            live += [(x - x0 + X, y - y0 + Y) for x, y in c]
            comments.append(f'({X},{Y}) {desc}')
            X += w + a.gap; rowh = max(rowh, h); maxw = max(maxw, X - a.gap)
    W, H = maxw, Y + rowh
    with open(out, 'w') as f:
        f.write(f'#N rarefind {a.rule}\n#C {a.soups} soups, seed {a.seed}; rare = seen in at most {lim:g} soups\n#C {census}\n')
        f.write('#C (x,y) = top-left corner of each object, x right, y down from the top-left of the pattern\n')
        for cm in comments:
            f.write(f'#C {cm}\n')
        f.write(f'x = {W}, y = {H}, rule = {a.rule}\n' + (rle_body(live) if live else '!') + '\n')
    print(f'{out}: ' + ', '.join(f'{len(v)} {k}' for k, v in sections.items() if v or k != 'fuses'), file=sys.stderr)


if __name__ == '__main__':
    main()
