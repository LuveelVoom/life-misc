#!/usr/bin/env python3
"""Generate random isotropic non-totalistic (INT / Hensel) rules, optionally restricted to part of rule space.

Each of the 102 transitions (51 birth, 51 survival) is switched on independently with probability -p,
after applying the restrictions:
  --exclude  transitions that are never on        (e.g. B0,B1,B2a  -- B1 means every B1 letter)
  --include  transitions that are always on       (e.g. S2,S3)
  --within   a rule: only its transitions may be on (an upper bound)
  --contains a rule: all its transitions are on     (a lower bound)
Rules are printed one per line in Golly notation (B2-a3/S23 style; each count written in the shorter of
the positive or the "-" form).

Examples:
  python3 randrule.py -n 1000 --exclude B0,B1,B2a --seed 1 > rules.txt
  python3 randrule.py -n 20 --within B2345678/S012345678 --contains B3/S23 -p 0.2
  python3 randrule.py -n 250 --totalistic --exclude B0,B1,B2      # Life-like rules
"""
import argparse, random, re, sys

# Hensel letters per neighbour count (Golly's current letter set).
LETTERS = ['', 'ce', 'cekain', 'cekainyqjr', 'cekainyqjrtwz', 'cekainyqjr', 'cekain', 'ce', '']
TRANSITIONS = [(bs, n, l) for bs in 'BS' for n in range(9) for l in (LETTERS[n] or [''])]  # 102


def parse_tokens(spec):
    """'B0,B1,B2a S2-in' -> set of transitions. A token is B or S, then digits each optionally followed by
    letters (or '-' and the letters to leave out), as in a rule half."""
    out = set()
    for tok in re.split(r'[,\s]+', spec.strip()):
        if tok:
            out |= parse_half(tok[0].upper(), tok[1:], tok)
    return out


def parse_half(bs, s, ctx):
    if bs not in 'BS':
        sys.exit(f'bad transition list "{ctx}": must start with B or S')
    out, i = set(), 0
    while i < len(s):
        if not s[i].isdigit():
            sys.exit(f'bad transition list "{ctx}" at "{s[i:]}"')
        n = int(s[i]); i += 1
        neg = i < len(s) and s[i] == '-'
        if neg:
            i += 1
        j = i
        while j < len(s) and s[j].isalpha():
            j += 1
        given = s[i:j]; i = j
        for c in given:
            if c not in LETTERS[n]:
                sys.exit(f'"{ctx}": {bs}{n} has no letter {c} (letters: {LETTERS[n] or "none"})')
        if not given:
            letters = LETTERS[n] or ['']
        elif neg:
            letters = [c for c in LETTERS[n] if c not in given]
        else:
            letters = list(given)
        out |= {(bs, n, l) for l in letters}
    return out


def parse_rule(rule):
    m = re.fullmatch(r'\s*[Bb]([0-8a-z-]*)/[Ss]([0-8a-z-]*)\s*', rule)
    if not m:
        sys.exit(f'cannot read rule "{rule}" (use B.../S... Hensel notation)')
    return parse_half('B', m.group(1), rule) | parse_half('S', m.group(2), rule)


def rule_string(on):
    half = {}
    for bs in 'BS':
        s = ''
        for n in range(9):
            if not LETTERS[n]:
                s += str(n) if (bs, n, '') in on else ''
                continue
            have = [l for l in LETTERS[n] if (bs, n, l) in on]
            if not have:
                continue
            miss = [l for l in LETTERS[n] if l not in have]
            if not miss:
                s += str(n)
            elif len(miss) < len(have):
                s += f'{n}-{"".join(miss)}'
            else:
                s += f'{n}{"".join(have)}'
        half[bs] = s
    return f'B{half["B"]}/S{half["S"]}'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-n', type=int, default=10, help='number of rules (default 10)')
    ap.add_argument('-p', type=float, default=0.5, help='probability each free transition is on (default 0.5)')
    ap.add_argument('--exclude', default='', help='transitions never on, e.g. B0,B1,B2a')
    ap.add_argument('--include', default='', help='transitions always on, e.g. S23')
    ap.add_argument('--within', help='only transitions of this rule may be on')
    ap.add_argument('--contains', help='every transition of this rule is on')
    ap.add_argument('--seed', type=int, help='random seed (default: random)')
    ap.add_argument('--allow-duplicates', action='store_true', help='do not skip rules already printed')
    ap.add_argument('--totalistic', action='store_true',
                    help='outer-totalistic (Life-like) rules: each neighbour count is on or off as a whole')
    a = ap.parse_args()

    allowed = set(TRANSITIONS) - parse_tokens(a.exclude)
    if a.within:
        allowed &= parse_rule(a.within)
    forced = parse_tokens(a.include) | (parse_rule(a.contains) if a.contains else set())
    if forced - allowed:
        sys.exit(f'--include/--contains conflict with --exclude/--within: {rule_string(forced - allowed)}')
    free = sorted(allowed - forced, key=TRANSITIONS.index)
    if a.totalistic:   # units are whole counts: allowed only if every letter is allowed
        groups = {}
        for t in TRANSITIONS:
            groups.setdefault(t[:2], []).append(t)
        part = [g for g, ts in groups.items() if 0 < len(set(ts) & forced) < len(ts)]
        if part:
            sys.exit('--totalistic: --include/--contains force only part of ' + ', '.join(f'{b}{n}' for b, n in part))
        free = [frozenset(ts) for g, ts in groups.items() if set(ts) <= allowed and not set(ts) & forced]
    else:
        free = [frozenset([t]) for t in free]
    space = 2 ** len(free)
    if not a.allow_duplicates and a.n > space:
        sys.exit(f'only {space} rules fit these restrictions; asked for {a.n}')
    print(f'# {len(free)} free transitions, {len(forced)} forced on; {space:.3g} possible rules', file=sys.stderr)

    rng = random.Random(a.seed)
    seen, printed = set(), 0
    while printed < a.n:
        r = rule_string(forced.union(*[u for u in free if rng.random() < a.p]))
        if r in seen and not a.allow_duplicates:
            continue
        seen.add(r)
        print(r)
        printed += 1


if __name__ == '__main__':
    main()
