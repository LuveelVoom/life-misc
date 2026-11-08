"""Small self-contained INT (Hensel) rule engine for slwindows: parse a rule string into a
512-entry table (Golly MAP index order: NW 256, N 128, NE 64, W 32, C 16, E 8, SW 4, S 2, SE 1),
name it canonically, step a set of live cells, find a period, read RLE.

Parsing and naming follow Golly's liferules.cpp. other/slwindows/test_minirule.py (in life-misc) compares
everything here with casim, an independent simulator. Copied from other/slwindows; nothing
here needs casim.
"""
import re

LETTERS = "ceaiknjqrytwz"
# Hensel letter of each outer neighbourhood: index i in 0..511 with the centre bit (16) clear,
# in increasing order of i; '-' for counts 0 and 8, which have no letters.
_TABLE = ("-ceaccaieaeaknjaekejanaairerririccknncqnaijaqnwakykkqyqjrtjnzrqaekirkyrtejerkkjnekejjkrnej"
          "ecjyccanriqyzraariqjqajkjywkqkrnccqknccnkqccnnkqkqyykjaqjwinaarzjqtrnaccyyccyennkjyekenykkne"
          "jeirykrikeaqrznyirjwjqkkykaqrqajiarqcnnkccintrneriaanajekeajnkaeaeiaccaec-")
_OUTER = [i for i in range(512) if not i & 16]
LETTER = {i: _TABLE[n] for n, i in enumerate(_OUTER)}
COUNT = {i: bin(i).count("1") for i in _OUTER}
# VALID[k] = letters that exist for k live neighbours, in Golly's order
VALID = {k: [L for L in LETTERS if any(COUNT[i] == k and LETTER[i] == L for i in _OUTER)]
         for k in range(1, 8)}


class RuleError(ValueError):
    pass


def _apply(table, part, c):
    """Golly's setRuleFromString for one part (c = 0 birth, 1 survival)."""
    touched, sign = set(), {}
    i, n = 0, len(part)
    while i < n:
        ch = part[i]
        i += 1
        if not ch.isdigit() or int(ch) > 8:
            raise RuleError(f"unexpected {ch!r} in {part!r}")
        k = int(ch)
        nxt = part[i] if i < n else ""
        if not nxt or nxt not in LETTERS:
            if k not in touched:          # Golly skips a count already set
                touched.add(k)
                for o in _OUTER:
                    if COUNT[o] == k:
                        table[o | c << 4] = 1
        neg = nxt == "-"
        if neg:
            i += 1
        while i < n and part[i] in LETTERS:
            L = part[i]
            i += 1
            if L not in VALID.get(k, []):
                raise RuleError(f"letter {L!r} not valid after {k}")
            if sign.setdefault(k, neg) != neg:
                raise RuleError(f"count {k} mixes normal and negated letters")
            touched.add(k)
            for o in _OUTER:
                if COUNT[o] == k and LETTER[o] == L:
                    table[o | c << 4] = 0 if neg else 1


def parse(rule):
    """512-entry table of 0/1 for a B.../S... rule (Life-like or Hensel letters)."""
    t = rule.strip().lower().replace(" ", "")
    m = re.fullmatch(r"b([0-8a-z\-]*)/?s([0-8a-z\-]*)", t)
    if not m:
        raise RuleError(f"{rule!r}: expected B.../S... (Life-like or INT)")
    table = [0] * 512
    _apply(table, m.group(1), 0)
    _apply(table, m.group(2), 1)
    return table


def _part(table, c):
    out = ""
    for k in range(9):
        if k in (0, 8):
            if table[(0x1ef if k == 8 else 0) | c << 4]:
                out += str(k)
            continue
        on = [L for L in VALID[k] if any(table[o | c << 4] for o in _OUTER
                                          if COUNT[o] == k and LETTER[o] == L)]
        m = len(VALID[k])
        if not on:
            continue
        if len(on) == m:
            out += str(k)
        elif len(on) > m // 2 and not (len(on) == 7 and m == 13):
            out += str(k) + "-" + "".join(sorted(L for L in VALID[k] if L not in on))
        else:
            out += str(k) + "".join(sorted(on))
    return out


def name(table):
    """Golly's canonical name of an isotropic table."""
    return "B" + _part(table, 0) + "/S" + _part(table, 1)


def b0(table):
    return bool(table[0])


def step(live, table):
    """One generation of a set of (x, y) cells (y down)."""
    cand = {(x + dx, y + dy) for x, y in live for dx in (-1, 0, 1) for dy in (-1, 0, 1)}
    out = set()
    for x, y in cand:
        i = 0
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                i = 2 * i + ((x + dx, y + dy) in live)
        if table[i]:
            out.add((x, y))
    return out


def _key(live):
    if not live:
        return (), 0, 0
    x0 = min(x for x, _ in live)
    y0 = min(y for _, y in live)
    return tuple(sorted((x - x0, y - y0) for x, y in live)), x0, y0


def period(cells, table, max_gen=10000):
    """Same result as casim's Pattern.period for B0-free rules: None, or
    {period, dx, dy, start, kind} with kind 'dies', 'still life', 'oscillator', 'spaceship'."""
    live = set(cells)
    k, x0, y0 = _key(live)
    seen = {k: (0, x0, y0)}
    for g in range(1, max_gen + 1):
        if not live:
            return {"period": 1, "dx": 0, "dy": 0, "start": g - 1, "kind": "dies"}
        live = step(live, table)
        k, x, y = _key(live)
        hit = seen.get(k)
        if hit is not None:
            s, hx, hy = hit
            per, dx, dy = g - s, x - hx, y - hy
            kind = "spaceship" if dx or dy else ("still life" if per == 1 else "oscillator")
            return {"period": per, "dx": dx, "dy": dy, "start": s, "kind": kind}
        seen[k] = (g, x, y)
    return None


def parse_rle(text):
    """(cells, rule): live cells of a two-state RLE as (x, y), and the header's rule (or None)."""
    rule, body = None, []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s.lower().startswith("x"):
            m = re.search(r"rule\s*=\s*([^\s,]+)", s, re.I)
            rule = m.group(1) if m else None
            continue
        body.append(s)
    cells, x, y, num = [], 0, 0, ""
    for ch in "".join(body):
        if ch.isdigit():
            num += ch
            continue
        n = int(num) if num else 1
        num = ""
        if ch == "!":
            break
        if ch == "$":
            y += n
            x = 0
        elif ch in "bB.":
            x += n
        else:                 # o, A, or any other live state
            cells.extend((x + i, y) for i in range(n))
            x += n
    return cells, rule
