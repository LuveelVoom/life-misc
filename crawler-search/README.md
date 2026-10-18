# crawler-search

Brute-force search for **crawlers**: an active pattern collides with an object and comes back
unchanged but moved, while the object is left moved or rephased (LifeWiki definition; see the
saved `Crawler - LifeWiki.html`). Repeated on a line of objects, the reaction travels forever.

```
python3 crawler_search.py --rule B3/S23-q4i --object '2o$2o!' --active 'o$3o$obo$2bo!' --max-period 60
python3 crawler_search.py --rule B3/S23-q4i --object '2o$2o!' --region -5 -1 4 4 --max-period 30
```

Two ways to choose the starting patterns:
- `--active RLE` (repeatable): your own active patterns, each tried in all 8 orientations at every
  position within `--distance` cells of the object (default 8), or only at the top-left positions
  in `--positions X Y W H`. `--as-given` skips the other orientations. Fast: a pattern near a
  blinker at distance 12 is about 3,500 placements (~40 s at max-period 60).
- `--region X Y W H`: every pattern of the cells in a box (blind enumeration; doubles per cell).

**Line mode** (`--active ... --line DX DY`): for crawlers that only start repeating after a
startup phase. The pattern is run in front of a line of objects in direction (DX, DY) (entries
-1/0/1; `0 -1` = objects above each other, crawling up), for every spacing in `--spacing MIN MAX`
(default 3 20) and, for oscillator objects, every phase step between neighbouring objects. A
crawler is reported when the board near the front of the reaction (`--window` cells behind the
frontmost change, default 2 spacings, up to one spacing ahead) repeats moved along the line, twice
in a row with the same period, within `--max-gen` generations (default 4 x max-period). Debris
left further behind is allowed. Each find is re-checked exactly before it is reported.
Example: a bare pi heptomino below a column of blinkers finds the 17c/45 crawler (repeating from
generation 10) in about 90 s:
```
python3 crawler_search.py --object '3o!' --object-pos 13 18 --active '3o$obo$obo!' \
    --positions 9 20 9 30 --line 0 -1 --spacing 17 17 --max-period 50 --max-gen 300 --window 12
```
Cost grows with placements x spacings x phase steps: about 45 s per spacing and phase step for
1,000 placements at max-gen 300, so narrow `--spacing` and `--positions` when you can.

The simulation board is not limited to these: it has a margin of max-period + 2 cells on every
side, so nothing can reach its edge, and the line check runs on an unbounded board.

Needs Python 3 with numpy. `hensel.json` (INT neighbourhood table) must sit next to the script.
Rules: Hensel/INT notation (`B3/S23-q4i`, `B36/S23`, `B2-a/S12`); B0 rules are not supported.

| Option | Meaning |
|---|---|
| `--object RLE` | the object crawled over (still life, oscillator or spaceship) |
| `--object-pos X Y` | its top-left cell (default `0 0`); x right, y down |
| `--active RLE` | an active pattern to place (repeatable); with `--distance N`, `--positions X Y W H`, `--as-given` |
| `--region X Y W H` | instead of `--active`: every pattern of the box's cells (object cells excluded) |
| `--max-period N` | longest cycle tried (default 30) |
| `--max-cells N` | only patterns with at most N live cells (for bigger regions) |
| `--line DX DY` | line mode (with `--active`); `--spacing MIN MAX`, `--max-gen N`, `--window N` |
| `--check-cycles N` | cycles a find must keep repeating for (default 10) |
| `--allow-debris` | the reaction may also leave debris/gliders behind (slower) |
| `--out FILE` | results (default `crawlers.rle`) |

**What counts.** After t generations the board must be exactly the starting pattern moved by d
plus the object (any phase) moved by e, with d not 0 and not the object's own movement. With
`--allow-debris`, the pattern only has to come back moved with nothing else within 2 cells;
finds are then labelled crawler (a copy of the object is among the remains) or fuse (the object
is used up). Every find must repeat on a line of four objects (spacing d, phases chosen so each
object is in its starting phase when the pattern arrives): after each of 3 cycles the pattern is
back and the objects ahead are untouched. Debris behind the pattern is not followed further.

**Output.** One row per result in the RLE file: the pattern plus the line of four objects, ready
to run in Golly; `#C` lines give period, speed and what happens to the object.

**Speed.** Every pattern is simulated in numpy batches (only the box around live cells is
stepped; patterns that die, settle into a still/period-2 board or reach the canvas edge are
dropped). About 65k patterns (a 4x4 region) take ~9 s at max-period 30. More than 50 million
patterns is refused: use a smaller region or `--max-cells`.

**Checked (2026-10-18).** The LifeWiki 17c/45 pi heptomino reaction (pi phase + one blinker)
is recognised with `--allow-debris` (period 45, moves 17) and passes the line test. Its own
pattern is 22 cells wide, so it is not reachable by brute force from a small region; the 7-cell
pi heptomino itself never comes back isolated in that cycle. Searches in Life near blinkers
find a 4c/5 blinker fuse (`3o$3o$bo!` on horizontal blinkers 4 apart; checked over 9 cycles).
The strict (no-debris) mode found nothing in 27 test searches (B3/S23-q4i, B36/S23, B3/S23-k;
block, blinker, beehive; 4x4 regions on three sides); clean crawlers are rare.

**Spaceships.** The object may be a spaceship (e.g. a glider stream in line mode, for climbers).
Crawlers that emit spaceships are kept: the board's edges absorb whatever escapes, and in line
mode only the front of the reaction has to repeat. Fake crawlers are rejected: a front that moves
the same way without the objects (e.g. a spaceship flying alongside the line), or one that never
gets through any objects.

**Labels (line mode).** Each find is reduced to the smallest period of its front and listed once
per distinct front, labelled `crawler` (the objects it passed are still there behind it, moved or
not) or `fuse` (they are used up).

**Limits.** One object type. Without `--line`, the starting pattern must itself come back (moved),
so crawlers with a startup phase (such as the bare pi heptomino) need line mode. The repeating
part must fit in `--window`. A spaceship fired forward along the line would be taken as the front.

**B3/S23-a5 test (2026-10-18), B-heptomino `b3o$o2bo$o2bo!`.** Common objects from a 300-soup
census: blinker, a p2 2x3 oscillator (`2o$2o$2o!`), loaf, boat, pond, beehive, tub, ship; glider
(p4) also tested. Basic mode (distance 10, max-period 60): a (12,2)c/49 fuse through ships that
holds for 9 cycles, then its growing debris catches up (rejected at the default 10 cycles).
Line mode (upward line, spacings 6-14, small position box, max-gen 150, ~2 min): only light-speed
fuses along loaves 6 apart (front periods 6, 12 and 24). No crawlers in this short test. Results that are the same crawler in a
different phase of its cycle are listed separately.
