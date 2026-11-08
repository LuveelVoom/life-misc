"""python3 test_doublerule.py  (about 2 s)"""
import itertools
import doublerule as D
import minirule as M


def ev(nodes, xs):
    n = nodes[-1]
    for x in xs[:-1]:
        n = nodes[n[1 + x]]
    return n[1 + xs[-1]]


def test_b3s23_matches_forum_tree():
    # DoubleB3S23 from viewtopic.php?p=22316 has 130 nodes and its root is "9 104 114 124 128"
    _, nodes = D.build_tree(M.parse("B3/S23"))
    assert len(nodes) == 130 and nodes[-1] == (9, 104, 114, 124, 128)


def test_layers_independent_int():
    for r in ["B2-a3/S23", "B2ce3ai/S1e2-k3", "B0/S8"]:
        t = M.parse(r)
        _, nodes = D.build_tree(t)
        for xs in itertools.product(range(4), repeat=9):
            s = ev(nodes, xs)
            for l in (0, 1):
                i = sum(w for w, x in zip(D.TREE_WEIGHTS, xs) if x >> l & 1)
                assert (s >> l & 1) == t[i], (r, xs, l)


if __name__ == "__main__":
    test_b3s23_matches_forum_tree()
    test_layers_independent_int()
    print("ok")
