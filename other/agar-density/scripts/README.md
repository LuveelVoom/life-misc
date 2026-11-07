# Oscillating-agar density scripts

Setup (fresh environment recommended):

    conda create -n agar python=3.12 -y && conda activate agar
    pip install numpy scipy highspy python-sat python-flint ortools

| Script | What it does | Typical command |
|---|---|---|
| `bound_16_31.py` | Standalone 3x3 all-period LP: optimum 16/31 and rigorous bound 0.516129050 | `python bound_16_31.py` |
| `life_agar_density.py` | 3x3 windows: peeling + exact rational certificate (main theorem, 769 excluded windows); band embedding tests | `python life_agar_density.py certify` |
| `big_window_lp.py` | Any window size h x w: solve, peel, rigorous certificates, mirror symmetry (`--sym`), MPS export | `python big_window_lp.py solve --h 3 --w 4 --sym --threads 16` |
| `strip_certificate.py` | Exact proof (or counterexample) for all period-p agars repeating every W columns | `python strip_certificate.py --width 6 --period 2 --jobs 16` |
| `band_forcing_check.py` | Solver check of the band-forcing lemma without assuming periodicity | `python band_forcing_check.py` |
| `check_windows.py` | Reports where the excluded window types occur in a given oscillator (edit the `tests` list to add your own) | `python check_windows.py` |
| `coverage_search2.py` | CP-SAT search for counterexamples to the charge-sharing ("coverage") rule; needs `certificate.pkl` from `life_agar_density.py certify` | `python coverage_search2.py --type 0 --radius 1 --decide` |

Notes
- `coverage_search2.py` imports `life_agar_density.py`; keep them in the same folder.
- `big_window_lp.py certify` reuses dual values saved by `solve`/`peel` (no re-solve).
- Expected outputs for the main commands are listed in the write-up (Section 10).
- The first coverage search script (`coverage_search.py`) is superseded by `coverage_search2.py` and not included.
