# crawler-search

Finds **crawlers** and **fuses** in Life-like and INT rules by brute force.

- A **crawler** (LifeWiki definition) is a reaction where an active pattern hits an object and
  comes back unchanged but moved, while the object is left moved or rephased. Repeated on a line
  of objects, the pattern travels forever. Example: the 17c/45 pi heptomino crawler over blinkers.
- A **fuse** is the same, except that the objects are used up (burnt) instead of moved.

You give the script one object (the thing crawled over) and either your own active pattern(s) or
a box of cells to try every pattern in. It simulates every placement and reports the reactions
that keep repeating along a line of copies of the object.

Requirements: Python 3 with numpy. `hensel.json` must sit next to the script. Rules use Hensel/INT
notation (`B3/S23`, `B36/S23`, `B3/S23-q4i`, `B2-a/S12`); B0 rules are not supported.
Coordinates: x to the right, y down.

## Quick start

```
python3 crawler_search.py --rule B3/S23-a5 --object '2o$obo$b2o!' \
    --active 'b3o$o2bo$o2bo!' --distance 10 --max-period 60 --allow-debris
```

This places the B-heptomino in every orientation at every position within 10 cells of a ship,
runs each for up to 60 generations, and writes what it finds to `crawlers.rle`. It finds a
(12,2)c/49 fuse that burns through a line of ships (it lasts 9 cycles, so it is only reported
with `--check-cycles 9`).

Open `crawlers.rle` in Golly. Each result is one row: the pattern plus a line of objects, ready
to run. The `#C` lines at the top say what each row does (period, speed, crawler or fuse).

## Three ways to search

### 1. Give the active pattern (`--active`): the pattern must come back by itself

```
python3 crawler_search.py --rule B3/S23 --object '3o!' --active '3o$3o$bo!' --max-period 12 --allow-debris
```

Tries the pattern in all 8 orientations at every position within `--distance` cells of the
object (default 8), and reports it if, after some period, the board is the pattern moved plus
the object moved/rephased. This example finds a 4c/5 fuse that burns a line of blinkers.

- `--positions X Y W H` tries only top-left positions in that box instead of `--distance`.
- `--as-given` uses the pattern only as given (no other orientations).
- `--active` can be repeated to try several patterns in one run.

### 2. Try everything in a box (`--region`): when you do not know the pattern

```
python3 crawler_search.py --rule B3/S23-q4i --object '2o$2o!' --region -5 -1 4 4 --max-period 30
```

Tries every pattern of live cells in the 4x4 box at (-5,-1) next to a block (65,535 patterns,
about 9 s). The count doubles with each extra cell; above 50 million patterns the script refuses.
`--max-cells N` keeps only patterns with at most N live cells, which makes bigger boxes feasible.

### 3. Line mode (`--active ... --line DX DY`): crawlers with a startup phase

Some crawlers only settle into their cycle after a while; the pi heptomino never reappears as
itself in the 17c/45 crawler. Line mode puts the pattern in front of a whole line of objects and
waits for the front of the reaction to repeat:

```
python3 crawler_search.py --rule B3/S23 --object '3o!' --object-pos 13 18 --active '3o$obo$obo!' \
    --positions 9 20 9 30 --line 0 -1 --spacing 17 17 --max-period 50 --max-gen 300 --window 12
```

This finds the 17c/45 pi crawler from a bare pi heptomino in a few seconds.

- `--line DX DY`: direction from one object to the next (entries -1, 0, 1). `0 -1` = objects
  stacked upwards (crawling up), `1 0` = sideways, `1 -1` = diagonal. Put the pattern on the side
  the line starts from (below it for `0 -1`, left of it for `1 0`, below-left for `1 -1`).
- `--spacing MIN MAX`: object spacings to try, in steps of the direction (default 3 20).
  For oscillators and spaceships every phase step between neighbours is also tried.
- `--max-gen N`: generations to simulate (default 4 x max-period).
- `--window N`: how far behind the front must repeat (default 2 spacings).

Cost grows with placements x spacings x phase steps; narrow `--spacing` and `--positions` when
you can.

## Searching many objects at once

The script takes one object per run. Run several in parallel from the shell:

```
P='2bo$2ob2o$bobo$2bo!'
for o in '3o!' '2o$2o$2o!' '2bo$bobo$o2bo$b2o!' '2o$obo$bo!' 'b2o$o2bo$o2bo$b2o!' \
         'b2o$o2bo$b2o!' 'bo$obo$bo!' '2o$obo$b2o!' 'bo$2bo$3o!' '3o$bob!'; do
  name=$(echo "$o" | tr -d '$!')
  python3 crawler_search.py --rule B3/S23-a5 --object "$o" --active "$P" \
      --positions -8 3 12 12 --line 0 -1 --spacing 4 16 --max-period 60 --max-gen 200 \
      --out "up_$name.rle" > "up_$name.log" &
done; wait
```

Each run uses one core. With 10 objects this took about 9 minutes on a 14-core machine. Change
`--line` and `--positions` for the other directions:

| Line | Example `--positions` (object at 0 0) |
|---|---|
| `--line 0 -1` (up) | `-8 3 12 12` (below the object) |
| `--line 1 0` (right) | `-14 -6 12 12` (left of the object) |
| `--line 1 -1` (up-right) | `-12 2 12 12` (below-left) |

## Reading the results

Each `#C` line describes one result, for example:

```
#C 1: fuse; period 5, passes 1 object(s) per period, moves (0,-4) = 4c/5; objects every 4 cells; ...
```

- `crawler`: the objects it passed are still there behind it (moved or not).
- `fuse`: the objects are used up.
- `crawler or fuse (N objects left behind)`: some object-shaped remains; check in Golly.
- `objects every N cells, each K generation(s) on from the last`: the line it needs.

Every result has been checked to repeat for `--check-cycles` cycles (default 10) on a line of
objects. Results are reduced to the smallest period of the front and listed once per distinct
front (mirror images and the same reaction from other starting positions are merged).

## What counts as a result

- **Default (`--active` or `--region`)**: after t generations the board is exactly the pattern
  moved by d plus the object moved by e (any phase), with d not 0 and not the object's own
  movement.
- **`--allow-debris`**: the pattern only has to come back moved, with nothing else within 2 cells
  of it; the rest (debris, emitted spaceships) is left behind.
- **Line mode**: the front of the reaction repeats, moved along the line; anything left behind is
  allowed.

Rejected in all modes: reactions where the "crawler" is just something moving the same way
without the objects (for example a spaceship flying alongside the line), and reactions that never
get through any objects. Crawlers that emit spaceships are kept: the board edges absorb whatever
escapes.

## All options

| Option | Default | Meaning |
|---|---|---|
| `--rule RULE` | `B3/S23` | Hensel/INT rule |
| `--object RLE` | (required) | the object crawled over: still life, oscillator or spaceship |
| `--object-pos X Y` | `0 0` | top-left cell of the object |
| `--active RLE` | | an active pattern to place (repeatable) |
| `--distance N` | 8 | `--active`: positions up to N cells from the object |
| `--positions X Y W H` | | `--active`: only these top-left positions |
| `--as-given` | off | `--active`: no other orientations |
| `--region X Y W H` | | instead of `--active`: every pattern of the box's cells |
| `--max-cells N` | 0 (no limit) | `--region`: at most N live cells |
| `--max-period N` | 30 | longest cycle tried |
| `--allow-debris` | off | the reaction may leave debris or spaceships behind |
| `--line DX DY` | | line mode (needs `--active`) |
| `--spacing MIN MAX` | `3 20` | line mode: object spacings |
| `--max-gen N` | 4 x max-period | line mode: generations simulated |
| `--window N` | 2 spacings | line mode: depth of the repeating front |
| `--check-cycles N` | 10 | cycles a result must keep repeating for |
| `--batch N` | 2048 | patterns simulated at once (memory) |
| `--out FILE` | `crawlers.rle` | output file |

## Limits

- One object type per run.
- Without `--line`, the starting pattern itself must come back; use line mode for crawlers with a
  startup phase.
- In line mode the repeating part must fit in `--window`, and a spaceship fired forward along the
  line would be mistaken for the front.
- The labels crawler/fuse are a quick check (are there object-shaped remains?); look at
  interesting results in Golly.

## Results so far (`results/`)

All in B3/S23-a5, searched on 2026-10-18 against the common objects of that rule (blinker, p2
2x3 oscillator `2o$2o$2o!`, loaf, boat, pond, beehive, tub, ship) plus the glider and the c/4
ship `3o$bob!`.

- **B-heptomino** `b3o$o2bo$o2bo!`: a (12,2)c/49 fuse through ships (9 cycles only) and
  light-speed fuses along loaves.
- **`2bo$2ob2o$bobo$2bo!`** (`results/a5_7cell_results*.txt`, one-line RLEs):
  - a possible climber on gliders, (5,10)c/42;
  - oblique fuses through boats, (6,10)c/39 and (7,10)c/54;
  - a 4c/5 fuse that turns a column of beehives into a trail of beehive pairs;
  - oblique fuses along a row of c/4 ships, (4,1)c/4, (5,2)c/6 and (5,2)c/5;
  - 5c/7 and 5c/14 fuses along the p2 2x3 oscillator; 4c/5 fuses on tubs, loaves and ponds;
    3c/4 diagonal fuses on boats and ships; many light-speed fuses.
  - No true crawlers.
