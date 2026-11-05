# Small rule-space exploration

Explosiveness survey of the 2^16 INT rules between `B3-qr5ce6c/S23-a4er5-k` (ID 0x0000)
and `B3-q4ik5cejkr6ce7/S23-a4ejrz56aek7c` (ID 0xFFFF); bit i of a rule ID adds one
transition (list in `rulespace.py`).

| Path | What it is |
|---|---|
| `rulespace.py` | Rule IDs ↔ Hensel strings ↔ MAP tables, and the driver for `soupclass`. `python3 rulespace.py --help` |
| `soupclass.c` | Soup engine that classifies each rule (built as `./soupclass`, found next to the script) |
| `test_rulespace.py`, `test_engine.py` | Tests |
| `results/` | Runs: `explosiveness_100` (100 random rules), `noB3rB4i` (16384 rules without B3r/B4i), `noB3rB4i_variety` |

Output goes to `--out NAME` (`NAME.md` + `NAME.csv`) in the current folder, so run from
`results/` (`python3 ../rulespace.py ...`) to keep it there.
