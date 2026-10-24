# randrule

Generates random isotropic non-totalistic (INT / Hensel) rules, one per line in Golly notation. You can limit
which part of rule space they come from. Python 3, no dependencies.

```
python3 randrule.py -n 1000 --exclude B0,B1,B2a --seed 1 > rules.txt
python3 randrule.py -n 20 --within B2345678/S012345678 --contains B3/S23 -p 0.2
```

Each of the 102 transitions (51 birth, 51 survival) is switched on independently with probability `-p`
(default 0.5), after the restrictions:

| Option | Meaning |
|---|---|
| `-n N` | number of rules (default 10); duplicates are skipped unless `--allow-duplicates` |
| `-p P` | chance each free transition is on. Lower values give sparser rules. |
| `--exclude LIST` | transitions that are never on, e.g. `B0,B1,B2a` (`B1` = every B1 letter) |
| `--include LIST` | transitions that are always on, e.g. `S23` |
| `--within RULE` | only transitions of RULE may be on (upper bound) |
| `--contains RULE` | every transition of RULE is on (lower bound) |
| `--totalistic` | outer-totalistic (Life-like) rules: each neighbour count is switched as a whole |
| `--seed S` | repeatable output |

Lists use rule-half syntax: `B2a`, `B2-a` (every B2 letter except a), `S2ce3`, comma or space separated.
How many free transitions and possible rules are left is printed to stderr.

Each rule is written with the shorter of the positive and the `-` form for each count (`B2-a` rather than
`B2cekin`) [1].

B0 rules can be generated (`B0` isn't excluded unless you say so), but many tools reject them.

## Checks

1. Checked against casim: 3000 random transition sets give the same rule tables as the fully spelled-out form, and generated rules respect `--exclude`, `--within` and `--contains`.
