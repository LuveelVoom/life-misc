# Handoff: Game of Life with light-speed delay

Paste this file into a new chat (and attach `life-lightcone.html`) to continue.

## What was asked (in order)
1. "Code a model of Conway's Game of Life viewed from the perspective of a virtual observer (which can be given specified speed or be stationary) and information about distant objects only propagates at c and is not instantly seen. Make the starting pattern have several unique oscillators and spaceships to help make the effect clear."
2. "Make the world bigger, add more objects with unique properties, and increase the maximum observer speed."
3. Export everything to continue on another account (this bundle).

## Current state: `life-lightcone.html`
A single self-contained HTML page (no build step; open in a browser). Fonts load from Google Fonts; everything else is inline.

- **Two side-by-side views:** the instantaneous lattice (left) and the observer's delayed view (right).
- **Delay model:** the observer at time T sees cell (x, y) as it was at generation ⌊T − d/c⌋, where d is the straight-line distance from the observer's current position. This is classical retardation in the lattice rest frame, with no Lorentz effects; that was a deliberate modelling choice because Life has a preferred frame.
- **World:** 300×210 toroidal grid (the edges wrap). Light is traced only along direct paths inside the rectangle, not around the wrap.
- **History:** 1280 generations stored bit-packed, about 10 MB.
- **Observer:**
  - set by click or drag
  - velocity vx, vy from −4 to 4 cells/gen
  - bounces off the walls, with sub-stepping at high speed
- **Signal speed:** c from 0.3 to 4.
- **Display options:** delay colour ramp (0 to 200+ gen), delay rings every 20 gen, "Stack views" toggle.
- **Presets:**
  - centre
  - far corner
  - chase (0.5c)
  - fly-by (0.9c)
  - outrun the light (2.5c)
  - slow light (c = 0.4)
- **Diehard respawn toggle:** respawns the diehard every 300 gen.

### Objects (layout lives in the `FIXED` / `MOVERS` / `DIEHARD` arrays)
- **Gun systems:** two Gosper guns firing toward each other in a V. Each stream is absorbed by an eater placed at gun + (90+k, 76+k), with k = 6.
- **Fixed oscillators:**
  - queen bee shuttle (p30)
  - pulsar ×3 (p3)
  - tumbler (p14)
  - figure eight (p8)
  - octagon 2 (p5)
  - pentadecathlon ×2 (p15)
  - Kok's galaxy ×2 (p8)
  - beacon, toad, blinker (p2)
- **Moving objects:**
  - LWSS ×3, MWSS, HWSS ×2 (c/2, two-way traffic in the top lanes)
  - weekender (2c/7)
  - loafer (c/7)
  - copperhead (c/10)

## Verification already done (headless Node)
- **Pattern behaviour:** every pattern's period or velocity was checked by simulation.
- **Bugs found and fixed:**
  - The first MWSS RLE moved diagonally and was replaced.
  - The first tumbler RLE did not oscillate and was replaced.
  - With the original 180×120 dead-boundary layout, border-erased ships left debris that eventually hit a gun. Switching to a wrap-around world fixed this.
- **Long-run stability:** the full world was run for 8000 generations.
  - Every fixed object stayed intact at its own period.
  - The populations of both ship-lane regions (top and bottom) matched from one full cycle to the next, so no collisions.
  - The diehard dies at generation 130, as expected.
- **Speed:** about 0.4 ms per generation in Node.

## Files in this bundle
- `life-lightcone.html`: the current deliverable (open in a browser).
- `dev/world.js`: the world block extracted verbatim from the page (between the `===== WORLD` markers), plus a `module.exports` line for Node.
- `dev/verify.js`: stability and period checks. Run `node verify.js 8000` inside `dev/`. It checks that every fixed object is intact at its own period, that ship-lane populations repeat, and each oscillator's period.
- The one-off scripts used earlier (eater placement search, ship velocity tests) were not preserved; their results are recorded above.

## Possible next steps (not requested yet)
- Proper relativistic aberration (an observer-frame transform).
- Light that propagates across the wrap-around edges.
- Pattern editor or RLE import.

## User's standing preferences
These were set in the original account's settings, not in the chat. Re-enter them in the new account under Settings → Profile if wanted:
- concise, no repetition
- no unfounded "deep" insights
- don't agree just to agree; evaluate corrections
- treat each chat fresh
- state confidence levels (high / moderate / low) and source types
