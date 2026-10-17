# Which of the excluded window types occur in a given oscillator?  Usage: python check_windows.py
import pickle, numpy as np, life_agar_density as L
W = L.reachable_windows()                                   # same window order as the certify run
excluded = pickle.load(open("excluded_windows.pkl", "rb"))  # boolean mask, one entry per window
print("excluded window types:", int(excluded.sum()))
index = {int(k): i for i, k in enumerate((W * (4 ** np.arange(9))).sum(1))}
def windows_of(phases):
    """window index at every (cell, time) of an oscillator given as a list of 0/1 arrays (one per phase)"""
    p = len(phases); m, n = phases[0].shape; out = []
    for t in range(p):
        a, b = phases[t], phases[(t + 1) % p]
        for i in range(m):
            for j in range(n):
                code = sum(int(a[(i+r-1) % m, (j+c-1) % n] | (b[(i+r-1) % m, (j+c-1) % n] << 1)) * 4 ** (r*3 + c)
                           for r in range(3) for c in range(3))
                out.append(index[code])
    return np.array(out)
g = lambda rows, alive: np.array([[1 if ch in alive else 0 for ch in r] for r in rows])
band = ["XXSYYS","XX.YY.",".SS.SS","YY.XX.","YYSXXS","XX.YY.","XX.YY.","YYSXXS","YY.XX.","SS.SS.","XX.YY.","XXSYYS","YY.XX.","YY.XX."]
tests = {"band agar (6x14, period 2)": [g(band, "XS"), g(band, "YS")],
         "stripes (still life)": [g(["oooo", "....", "oooo", "...."], "o")]}
for name, phases in tests.items():
    w = windows_of(phases); hits = excluded[w]
    density = np.mean([ph.mean() for ph in phases])
    print(f"{name}: density {density:.4f}; excluded windows at {int(hits.sum())} of {len(w)} (cell, time) positions")
