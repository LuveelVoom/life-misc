"""Edge stability of an agar: strip of H rows, wrapping in columns (period multiple),
for the agar and its 90-degree rotation. Reports the maximal self-forcing subset of the strip
per row (count forced per row)."""
import sys, sfcore as S
def rot(rows):  # rotate 90 clockwise
    a, b = len(rows), len(rows[0])
    return ["".join(rows[a - 1 - r][c] for r in range(a)) for c in range(b)]
def strip(rule, rows, H, reps=3, shift=0):
    a, b = len(rows), len(rows[0]); L = b * reps
    P = {(r, c): int(rows[(r + shift) % a][c % b] == 'o') for r in range(H) for c in range(L)}
    norm = lambda x: (x[0], x[1] % L)
    D = S.max_self_forcing(rule, P, set(P), log=None, norm=norm)
    return [sum((r, c) in D for c in range(L)) for r in range(H)], L
rule, rows, H = sys.argv[1], sys.argv[2].strip('/').split('/'), int(sys.argv[3])
for name, rr in [("rows (top/bottom edges)", rows), ("cols (left/right edges)", rot(rows))]:
    for sh in range(len(rr)):
        cnt, L = strip(rule, rr, H, shift=sh)
        print(f"{name} shift {sh}: width {L}, forced per row:", cnt, flush=True)
