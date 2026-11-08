# doublerule

`doublerule.py RULE` writes `Double<RULE>.rule`: a 4-state Golly rule tree holding two
independent copies of a Life-like or INT rule, laid out like DoubleB3S23
(viewtopic.php?p=22316#p22316). State 1 = layer A, 2 = layer B, 3 = both. The default name
is `Double` + Golly's canonical name with `/` removed (`B2-a3/S23` → `DoubleB2-a3S23`);
`-n` overrides it, `-o -` writes to stdout.

Tree neighbour order is Golly's `num_neighbors=8` order (NW, NE, SW, SE, N, W, E, S, C), taken
from lifelib's `genera/rulefiles/parsetree.py`. For B0 rules the tree is literal: Golly does not
emulate B0 for multi-state rules, so an empty background becomes state 3.

`minirule.py` is copied from `../slwindows`. `python3 test_doublerule.py` checks that B3/S23
reproduces the forum tree and that each layer of a few INT rules follows its own table.
