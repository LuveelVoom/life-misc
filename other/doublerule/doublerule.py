#!/usr/bin/env python3
"""Write a Golly .rule file for Double[rule]: two independent copies of a Life-like or INT
rule on one grid, as in DoubleB3S23 (viewtopic.php?p=22316#p22316).

State = bit 0 (layer A) | bit 1 (layer B): 0 empty, 1 A only, 2 B only, 3 both. Each layer
evolves under the rule while ignoring the other one.

    python3 doublerule.py B2-a3/S23            # writes DoubleB2-a3S23.rule
    python3 doublerule.py B3/S23 -o - -n MyName

Output is a @TREE. Neighbour order is Golly's for num_neighbors=8: NW, NE, SW, SE, N, W, E, S,
then C, the root splitting on NW (as read by lifelib's genera/rulefiles/parsetree.py).
"""
import argparse
import functools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import minirule  # noqa: E402

# MAP index weight of each tree variable, in tree order (NW NE SW SE N W E S C)
TREE_WEIGHTS = [256, 64, 4, 1, 128, 32, 8, 2, 16]


def build_tree(table, layers=2):
    """Rule tree for `layers` independent copies of a 512-entry table.
    Returns (num_states, list of node tuples); the root is the last node."""
    n_states = 1 << layers
    nodes, index = [], {}

    def node(t):
        if t not in index:
            index[t] = len(nodes)
            nodes.append(t)
        return index[t]

    @functools.lru_cache(maxsize=None)
    def recur(depth, idx):
        # idx[l] = partial MAP index of layer l from the variables fixed so far
        if depth == 8:  # last variable is C: emit a leaf of output states
            w = TREE_WEIGHTS[8]
            out = tuple(sum(table[idx[l] | (w if s >> l & 1 else 0)] << l for l in range(layers))
                        for s in range(n_states))
            return node((1,) + out)
        w = TREE_WEIGHTS[depth]
        kids = tuple(recur(depth + 1, tuple(idx[l] | (w if s >> l & 1 else 0)
                                            for l in range(layers)))
                     for s in range(n_states))
        return node((9 - depth,) + kids)

    recur(0, (0,) * layers)
    return n_states, nodes


def rule_file(name, source, n_states, nodes):
    lines = [f"@RULE {name}", "",
             f"Two independent copies of {source}: bit 1 = layer A, bit 2 = layer B.", "",
             "@TREE", "", f"num_states={n_states}", "num_neighbors=8", f"num_nodes={len(nodes)}"]
    lines += [" ".join(map(str, t)) for t in nodes]
    if n_states == 4:
        lines += ["", "@COLORS", "", "1 255 0 0", "2 0 255 0", "3 255 255 255"]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("rule", help="Life-like or INT rule, e.g. B2-a3/S23")
    ap.add_argument("-n", "--name", help="rule name (default Double + canonical name without '/')")
    ap.add_argument("-o", "--out", help="output file, '-' for stdout (default <name>.rule)")
    a = ap.parse_args()

    table = minirule.parse(a.rule)
    canon = minirule.name(table)
    if minirule.b0(table):
        print(f"warning: {canon} has B0. Golly runs multi-state rules literally, so an empty "
              "background turns into state 3 (no B0 emulation).", file=sys.stderr)
    name = a.name or "Double" + canon.replace("/", "")
    n_states, nodes = build_tree(table)
    text = rule_file(name, canon, n_states, nodes)
    out = a.out or name + ".rule"
    if out == "-":
        sys.stdout.write(text)
    else:
        with open(out, "w") as f:
            f.write(text)
        print(f"{out}: {len(nodes)} nodes", file=sys.stderr)


if __name__ == "__main__":
    main()
