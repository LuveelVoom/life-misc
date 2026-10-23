# Life misc

Smaller Conway's Game of Life tools and research. Each folder has its own `HANDOFF.md` with full details.

| Folder | What it is |
|---|---|
| [`lightcone/`](lightcone) | Game of Life seen by a moving observer, with information travelling at a finite speed c. Open `life-lightcone.html` in a browser; `dev/verify.js` checks the world in Node. |
| [`agar-density/`](agar-density) | Computer-assisted work on the conjecture that every oscillating agar has density at most 1/2: LP bounds (≈0.5146 for 3×4 windows), an exclusion theorem with an exact certificate, and strip certificates. Write-up in `writeup/`, scripts in `scripts/`. All results still need independent verification. |
| [`crawler-search/`](crawler-search) | Brute-force search for crawlers (a pattern that crawls over a line of objects, moving them) in Life-like and INT rules. `python3 crawler_search.py --help`. |
| [`ruleclimb/`](ruleclimb) | Constrained search over range-1 rules (Life-like, isotropic, MAP): hill-climb to make chosen objects soup-common, or explore for new ships/oscillators, puffers/guns, catalysts, signals and fuses. `make`, then see `ruleclimb/README.md`; results in `ruleclimb/findings/`. |
| [`transcensus/`](transcensus) | Port of the 2016 `transcensus.c` active-region soup census (`tc.zip`) to Life-like, INT and isotropic MAP rules, with deterministic threads; identical output to the original for B3/S23. See `transcensus/README.md`. |
| [`randrule/`](randrule) | Random INT rule generator, with options to limit the rules to part of rule space (exclude, include, within, contains). `python3 randrule/randrule.py --help`. |
| [`solver-explorer/`](solver-explorer) | Interactive charts of a long search-solver run log (partitions, branch depths, timings). Open `solver_run_explorer.html`; it can load any `.jsonl` log of the same format. |
