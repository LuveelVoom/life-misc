# Life misc

Smaller cellular-automaton tools and research. Each folder stands on its own; results go in its `results/`.

| Folder | What it is |
|---|---|
| [`lightcone/`](lightcone) | Life seen by a moving observer, with information travelling at a finite speed c. Open `life-lightcone.html`. |
| [`crawler-search/`](crawler-search) | Brute-force crawler search (a pattern crawling along a line of objects) in Life-like and INT rules. `python3 crawler_search.py --help`. |
| [`transcensus/`](transcensus) | Port of the 2016 `transcensus.c` active-region soup census (`orig/tc.zip`) to Life-like, INT and isotropic MAP rules, with deterministic threads; output identical to the original for B3/S23. |
| [`JLS-LLS/`](JLS-LLS) | JavaLifeSearch 1.7 with a "Run in LLS (INT rule)" command. Needs `JLS-1.7.zip` next to the folder (not in the repo). |

## `other/`

One-off research and small tools.

| Folder | What it is |
|---|---|
| [`agar-density/`](other/agar-density) | Computer-assisted work on the conjecture that every oscillating agar has density at most 1/2: LP bounds (≈0.5146 for 3×4 windows), an exclusion theorem with an exact certificate, and strip certificates. Write-up in `writeup/`, scripts in `scripts/`. Not yet independently checked. |
| [`intcover/`](other/intcover) | Smallest pattern containing all 102 INT transitions (annealer plus CP-SAT). Best so far: 68 cells in 15×15, `results/best.rle`. |
| [`slwindows/`](other/slwindows) | Small Life still lifes containing every realizable 3×3 window (SAT): 259 oriented windows, 51 classes. |
| [`smallspaceexploration/`](other/smallspaceexploration) | Explosiveness survey of a 2^16-rule INT space around B3/S23-a4er variants (`rulespace.py` + the `soupclass` C engine). `smallspaceexploration.zip` is an earlier snapshot (23 Oct). |
| [`solver-explorer/`](other/solver-explorer) | Interactive charts of a long solver run log (partitions, branch depths, timings). Open `solver_run_explorer.html`; it loads any `.jsonl` log of the same format. |

`ruleclimb` and `randrule` moved to their own repo, [LuveelVoom/ruleclimb](https://github.com/LuveelVoom/ruleclimb).
