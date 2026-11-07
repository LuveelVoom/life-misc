# Explosiveness of 100 rules

Space: `B3-qr5ce6c/S23-a4er5-k` (ID 0x0000) to `B3-q4ik5cejkr6ce7/S23-a4ejrz56aek7c` (ID 0xFFFF), inclusive. Bit i of the ID adds transition i: 0:B3r, 1:B4i, 2:B4k, 3:B5j, 4:B5k, 5:B5r, 6:B6e, 7:B7c, 8:B7e, 9:S4j, 10:S4z, 11:S5k, 12:S6a, 13:S6e, 14:S6k, 15:S7c.

Settings: 64 soups/rule, 16×16 at 50%, max 6000 gens, grid 1024², explosion = 512 occupied 16×16 tiles (or ≥ ¼ of that with near-quadratic growth over two time-doublings), seed 12345. Run time 500 s.

Scale: score = 10 + 90·P(explode) for rules with any exploding soup; for rules with none, 0–10 from the growing fraction and median lifespan. Ties broken by speed.

## Summary

| class | rules |
|---|---:|
| EXPLOSIVE | 76 |
| MOSTLY EXPLOSIVE | 7 |
| SOMETIMES EXPLOSIVE | 2 |
| RARELY EXPLOSIVE | 1 |
| GROWING | 2 |
| STABLE | 12 |

## Marginal effect of each transition (this sample)

Mean P(explode) of sampled rules with vs without the transition. Correlational, sample-dependent.

| bit | transition | with | without | Δ |
|---:|---|---:|---:|---:|
| 0 | B3r | 0.99 | 0.63 | +0.36 |
| 1 | B4i | 0.98 | 0.67 | +0.31 |
| 4 | B5k | 0.88 | 0.75 | +0.13 |
| 13 | S6e | 0.88 | 0.75 | +0.12 |
| 2 | B4k | 0.87 | 0.75 | +0.12 |
| 10 | S4z | 0.87 | 0.76 | +0.11 |
| 3 | B5j | 0.86 | 0.78 | +0.09 |
| 11 | S5k | 0.86 | 0.78 | +0.08 |
| 5 | B5r | 0.86 | 0.79 | +0.07 |
| 15 | S7c | 0.85 | 0.78 | +0.07 |
| 9 | S4j | 0.83 | 0.80 | +0.03 |
| 12 | S6a | 0.83 | 0.80 | +0.03 |
| 14 | S6k | 0.83 | 0.80 | +0.02 |
| 7 | B7c | 0.83 | 0.81 | +0.02 |
| 8 | B7e | 0.82 | 0.81 | +0.01 |
| 6 | B6e | 0.82 | 0.82 | +0.01 |

## All rules (most explosive first)

Columns: P = fraction of soups that exploded (95% Wilson interval); t_x = median generation at which an explosion was detected; v = median approximate expansion speed of exploding soups (c, from occupied-area growth); grow = soups still growing at the end without exploding; settled = became periodic or died; life = median lifespan of non-exploding soups.

| ID | rule | class | score | P (95% CI) | t_x | v | grow | settled | life |
|---|---|---|---:|---|---:|---:|---:|---:|---:|
| 0x76df | `B3-q4ik5cejk6ce7c/S23-a4ejrz5-k6aek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 800 | 0.15 | 0/64 | 0/64 | – |
| 0xd42f | `B3-q4ik5cejr6c/S23-a4erz5-k6ak7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 800 | 0.14 | 0/64 | 0/64 | – |
| 0x5f2f | `B3-q4ik5cejr6c7e/S23-a4ejrz56ak` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 800 | 0.14 | 0/64 | 0/64 | – |
| 0xfddb | `B3-q4i5cejk6ce7/S23-a4erz56aek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 832 | 0.14 | 0/64 | 0/64 | – |
| 0x3c5f | `B3-q4ik5cejk6ce/S23-a4erz56ae` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 832 | 0.13 | 0/64 | 0/64 | – |
| 0xd81f | `B3-q4ik5cejk6c/S23-a4er56ak7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 864 | 0.13 | 0/64 | 0/64 | – |
| 0xb0ff | `B3-q4ik5cejkr6ce7c/S23-a4er5-k6ae7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 896 | 0.13 | 0/64 | 0/64 | – |
| 0x9df7 | `B3-q4ik5cekr6ce7/S23-a4erz56a7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 896 | 0.13 | 0/64 | 0/64 | – |
| 0xcefe | `B3-qr4ik5cejkr6ce7c/S23-a4ejrz56k7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 800 | 0.13 | 0/64 | 0/64 | – |
| 0xbbf7 | `B3-q4ik5cekr6ce7/S23-a4ejr56ae7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 928 | 0.13 | 0/64 | 0/64 | – |
| 0xbabe | `B3-qr4ik5cejkr6c7c/S23-a4ejr56ae7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 928 | 0.12 | 0/64 | 0/64 | – |
| 0xe0bf | `B3-q4ik5cejkr6c7c/S23-a4er5-k6ek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 896 | 0.12 | 0/64 | 0/64 | – |
| 0x6b7f | `B3-q4ik5cejkr6ce7e/S23-a4ejr56ek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 896 | 0.12 | 0/64 | 0/64 | – |
| 0xdd93 | `B3-q4i5cek6c7/S23-a4erz56ak7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 896 | 0.12 | 0/64 | 0/64 | – |
| 0xfda9 | `B3-q5cejr6c7/S23-a4erz56aek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 896 | 0.12 | 0/64 | 0/64 | – |
| 0x6627 | `B3-q4ik5cer6c/S23-a4ejrz5-k6ek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 992 | 0.11 | 0/64 | 0/64 | – |
| 0x758e | `B3-qr4ik5cej6c7/S23-a4erz5-k6aek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 928 | 0.11 | 0/64 | 0/64 | – |
| 0x3457 | `B3-q4ik5cek6ce/S23-a4erz5-k6ae` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 896 | 0.11 | 0/64 | 0/64 | – |
| 0xb6ad | `B3-q4k5cejr6c7c/S23-a4ejrz5-k6ae7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 864 | 0.11 | 0/64 | 0/64 | – |
| 0xffac | `B3-qr4k5cejr6c7/S23-a4ejrz56aek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 928 | 0.11 | 0/64 | 0/64 | – |
| 0xfab9 | `B3-q5cejkr6c7c/S23-a4ejr56aek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 992 | 0.11 | 0/64 | 0/64 | – |
| 0xe623 | `B3-q4i5cer6c/S23-a4ejrz5-k6ek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 960 | 0.11 | 0/64 | 0/64 | – |
| 0xd81e | `B3-qr4ik5cejk6c/S23-a4er56ak7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1008 | 0.11 | 0/64 | 0/64 | – |
| 0x945e | `B3-qr4ik5cejk6ce/S23-a4erz5-k6a7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 928 | 0.11 | 0/64 | 0/64 | – |
| 0xbdb5 | `B3-q4k5cekr6c7/S23-a4erz56ae7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 928 | 0.11 | 0/64 | 0/64 | – |
| 0x44cb | `B3-q4i5cej6ce7c/S23-a4erz5-k6k` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 912 | 0.11 | 0/64 | 0/64 | – |
| 0x7756 | `B3-qr4ik5cek6ce7e/S23-a4ejrz5-k6aek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 960 | 0.11 | 0/64 | 0/64 | – |
| 0xd515 | `B3-q4k5cek6c7e/S23-a4erz5-k6ak7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 928 | 0.10 | 0/64 | 0/64 | – |
| 0x0ede | `B3-qr4ik5cejk6ce7c/S23-a4ejrz5` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1024 | 0.10 | 0/64 | 0/64 | – |
| 0xd3a9 | `B3-q5cejr6c7/S23-a4ejr5-k6ak7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1056 | 0.10 | 0/64 | 0/64 | – |
| 0x89de | `B3-qr4ik5cejk6ce7/S23-a4er57c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1056 | 0.10 | 0/64 | 0/64 | – |
| 0xbfd9 | `B3-q5cejk6ce7/S23-a4ejrz56ae7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1008 | 0.10 | 0/64 | 0/64 | – |
| 0xaa57 | `B3-q4ik5cek6ce/S23-a4ejr56e7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1120 | 0.10 | 0/64 | 0/64 | – |
| 0x00cf | `B3-q4ik5cej6ce7c/S23-a4er5-k` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 976 | 0.10 | 0/64 | 0/64 | – |
| 0x0d07 | `B3-q4ik5ce6c7e/S23-a4erz5` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1088 | 0.10 | 0/64 | 0/64 | – |
| 0x917d | `B3-q4k5cejkr6ce7e/S23-a4er5-k6a7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1008 | 0.10 | 0/64 | 0/64 | – |
| 0x3de5 | `B3-q4k5cer6ce7/S23-a4erz56ae` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1008 | 0.09 | 0/64 | 0/64 | – |
| 0x7431 | `B3-q5cekr6c/S23-a4erz5-k6aek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1088 | 0.09 | 0/64 | 0/64 | – |
| 0x9b53 | `B3-q4i5cek6ce7e/S23-a4ejr56a7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1088 | 0.09 | 0/64 | 0/64 | – |
| 0x3333 | `B3-q4i5cekr6c7e/S23-a4ejr5-k6ae` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1088 | 0.09 | 0/64 | 0/64 | – |
| 0x1643 | `B3-q4i5ce6ce/S23-a4ejrz5-k6a` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1152 | 0.09 | 0/64 | 0/64 | – |
| 0x11ad | `B3-q4k5cejr6c7/S23-a4er5-k6a` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1088 | 0.09 | 0/64 | 0/64 | – |
| 0xc32d | `B3-q4k5cejr6c7e/S23-a4ejr5-k6k7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1024 | 0.09 | 0/64 | 0/64 | – |
| 0x204f | `B3-q4ik5cej6ce/S23-a4er5-k6e` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1168 | 0.09 | 0/64 | 0/64 | – |
| 0xeb13 | `B3-q4i5cek6c7e/S23-a4ejr56ek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1168 | 0.09 | 0/64 | 0/64 | – |
| 0x1cbc | `B3-qr4k5cejkr6c7c/S23-a4erz56a` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1120 | 0.09 | 0/64 | 0/64 | – |
| 0xf652 | `B3-qr4i5cek6ce/S23-a4ejrz5-k6aek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1120 | 0.09 | 0/64 | 0/64 | – |
| 0xc795 | `B3-q4k5cek6c7/S23-a4ejrz5-k6k7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1168 | 0.08 | 0/64 | 0/64 | – |
| 0x6ee6 | `B3-qr4ik5cer6ce7c/S23-a4ejrz56ek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1104 | 0.08 | 0/64 | 0/64 | – |
| 0x0e83 | `B3-q4i5ce6c7c/S23-a4ejrz5` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1152 | 0.08 | 0/64 | 0/64 | – |
| 0xf1ca | `B3-qr4i5cej6ce7/S23-a4er5-k6aek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1216 | 0.08 | 0/64 | 0/64 | – |
| 0x5738 | `B3-qr5cejkr6c7e/S23-a4ejrz5-k6ak` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1296 | 0.08 | 0/64 | 0/64 | – |
| 0xc472 | `B3-qr4i5cekr6ce/S23-a4erz5-k6k7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1360 | 0.08 | 0/64 | 0/64 | – |
| 0x7634 | `B3-qr4k5cekr6c/S23-a4ejrz5-k6aek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1248 | 0.08 | 0/64 | 0/64 | – |
| 0xe409 | `B3-q5cej6c/S23-a4erz5-k6ek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1264 | 0.08 | 0/64 | 0/64 | – |
| 0x6133 | `B3-q4i5cekr6c7e/S23-a4er5-k6ek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1392 | 0.07 | 0/64 | 0/64 | – |
| 0x04b6 | `B3-qr4ik5cekr6c7c/S23-a4erz5-k` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1504 | 0.07 | 0/64 | 0/64 | – |
| 0xf04a | `B3-qr4i5cej6ce/S23-a4er5-k6aek7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1408 | 0.07 | 0/64 | 0/64 | – |
| 0xc965 | `B3-q4k5cer6ce7e/S23-a4er56k7c` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1344 | 0.07 | 0/64 | 0/64 | – |
| 0x5651 | `B3-q5cek6ce/S23-a4ejrz5-k6ak` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1488 | 0.07 | 0/64 | 0/64 | – |
| 0x300e | `B3-qr4ik5cej6c/S23-a4er5-k6ae` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1488 | 0.07 | 0/64 | 0/64 | – |
| 0x6935 | `B3-q4k5cekr6c7e/S23-a4er56ek` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1712 | 0.06 | 0/64 | 0/64 | – |
| 0x0fa9 | `B3-q5cejr6c7/S23-a4ejrz5` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1744 | 0.06 | 0/64 | 0/64 | – |
| 0x0b6d | `B3-q4k5cejr6ce7e/S23-a4ejr5` | EXPLOSIVE | 100 | 1.00 (0.94–1.00) | 1712 | 0.05 | 0/64 | 0/64 | – |
| 0x3742 | `B3-qr4i5ce6ce7e/S23-a4ejrz5-k6ae` | EXPLOSIVE | 99 | 0.98 (0.92–1.00) | 1920 | 0.05 | 0/64 | 1/64 | 130 |
| 0xda28 | `B3-qr5cejr6c/S23-a4ejr56ak7c` | EXPLOSIVE | 99 | 0.98 (0.92–1.00) | 2144 | 0.05 | 0/64 | 1/64 | 243 |
| 0x0e52 | `B3-qr4i5cek6ce/S23-a4ejrz5` | EXPLOSIVE | 99 | 0.98 (0.92–1.00) | 2176 | 0.05 | 0/64 | 1/64 | 183 |
| 0x2c45 | `B3-q4k5ce6ce/S23-a4erz56e` | EXPLOSIVE | 97 | 0.97 (0.89–0.99) | 1696 | 0.06 | 0/64 | 2/64 | 1075 |
| 0xc985 | `B3-q4k5ce6c7/S23-a4er56k7c` | EXPLOSIVE | 97 | 0.97 (0.89–0.99) | 1984 | 0.05 | 0/64 | 2/64 | 160 |
| 0x6651 | `B3-q5cek6ce/S23-a4ejrz5-k6ek` | EXPLOSIVE | 97 | 0.97 (0.89–0.99) | 2064 | 0.05 | 0/64 | 2/64 | 291 |
| 0xc25c | `B3-qr4k5cejk6ce/S23-a4ejr5-k6k7c` | EXPLOSIVE | 96 | 0.95 (0.87–0.98) | 1632 | 0.06 | 0/64 | 3/64 | 37 |
| 0xb131 | `B3-q5cekr6c7e/S23-a4er5-k6ae7c` | EXPLOSIVE | 96 | 0.95 (0.87–0.98) | 2208 | 0.04 | 0/64 | 3/64 | 214 |
| 0xcf14 | `B3-qr4k5cek6c7e/S23-a4ejrz56k7c` | EXPLOSIVE | 96 | 0.95 (0.87–0.98) | 2272 | 0.04 | 0/64 | 3/64 | 156 |
| 0xeb51 | `B3-q5cek6ce7e/S23-a4ejr56ek7c` | EXPLOSIVE | 96 | 0.95 (0.87–0.98) | 2464 | 0.04 | 0/64 | 3/64 | 24 |
| 0xb009 | `B3-q5cej6c/S23-a4er5-k6ae7c` | EXPLOSIVE | 96 | 0.95 (0.87–0.98) | 2656 | 0.04 | 0/64 | 3/64 | 176 |
| 0xa286 | `B3-qr4ik5ce6c7c/S23-a4ejr5-k6e7c` | EXPLOSIVE | 93 | 0.92 (0.83–0.97) | 2528 | 0.04 | 0/64 | 5/64 | 149 |
| 0xea92 | `B3-qr4i5cek6c7c/S23-a4ejr56ek7c` | MOSTLY EXPLOSIVE | 90 | 0.89 (0.79–0.95) | 2464 | 0.04 | 0/64 | 7/64 | 196 |
| 0xa98a | `B3-qr4i5cej6c7/S23-a4er56e7c` | MOSTLY EXPLOSIVE | 90 | 0.89 (0.79–0.95) | 2880 | 0.04 | 0/64 | 7/64 | 105 |
| 0xf5e0 | `B3-qr5cer6ce7/S23-a4erz5-k6aek7c` | MOSTLY EXPLOSIVE | 87 | 0.86 (0.75–0.92) | 2560 | 0.04 | 0/64 | 9/64 | 95 |
| 0x7c48 | `B3-qr5cej6ce/S23-a4erz56aek` | MOSTLY EXPLOSIVE | 78 | 0.75 (0.63–0.84) | 4560 | 0.03 | 0/64 | 16/64 | 146 |
| 0x885c | `B3-qr4k5cejk6ce/S23-a4er57c` | MOSTLY EXPLOSIVE | 76 | 0.73 (0.62–0.83) | 3968 | 0.03 | 5/64 | 12/64 | 73 |
| 0x7181 | `B3-q5ce6c7/S23-a4er5-k6aek` | MOSTLY EXPLOSIVE | 75 | 0.72 (0.60–0.81) | 4240 | 0.03 | 7/64 | 11/64 | 0 |
| 0x422c | `B3-qr4k5cejr6c/S23-a4ejr5-k6k` | MOSTLY EXPLOSIVE | 59 | 0.55 (0.43–0.66) | 3776 | 0.03 | 10/64 | 19/64 | 92 |
| 0xf9c8 | `B3-qr5cej6ce7/S23-a4er56aek7c` | SOMETIMES EXPLOSIVE | 47 | 0.41 (0.29–0.53) | 4672 | 0.03 | 21/64 | 17/64 | 0 |
| 0xe032 | `B3-qr4i5cekr6c/S23-a4er5-k6ek7c` | SOMETIMES EXPLOSIVE | 42 | 0.36 (0.25–0.48) | 4512 | 0.03 | 25/64 | 16/64 | 0 |
| 0xd420 | `B3-qr5cer6c/S23-a4erz5-k6ak7c` | RARELY EXPLOSIVE | 16 | 0.06 (0.02–0.15) | 3248 | 0.02 | 13/64 | 40/64 | 90 |
| 0x7520 | `B3-qr5cer6c7e/S23-a4erz5-k6aek` | GROWING | 2 | 0.00 (0.00–0.06) | – | – | 8/64 | 50/64 | 304 |
| 0xc958 | `B3-qr5cejk6ce7e/S23-a4er56k7c` | GROWING | 1 | 0.00 (0.00–0.06) | – | – | 7/64 | 52/64 | 401 |
| 0x5894 | `B3-qr4k5cek6c7c/S23-a4er56ak` | STABLE | 1 | 0.00 (0.00–0.06) | – | – | 2/64 | 54/64 | 374 |
| 0x02ec | `B3-qr4k5cejr6ce7c/S23-a4ejr5-k` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 246 |
| 0x064c | `B3-qr4k5cej6ce/S23-a4ejrz5-k` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 184 |
| 0x53d0 | `B3-qr5cek6ce7/S23-a4ejr5-k6ak` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 184 |
| 0x0f24 | `B3-qr4k5cer6c7e/S23-a4ejrz5` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 180 |
| 0xf848 | `B3-qr5cej6ce/S23-a4er56aek7c` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 122 |
| 0xb4e0 | `B3-qr5cer6ce7c/S23-a4erz5-k6ae7c` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 120 |
| 0x97c0 | `B3-qr5ce6ce7/S23-a4ejrz5-k6a7c` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 104 |
| 0x0114 | `B3-qr4k5cek6c7e/S23-a4er5-k` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 72 |
| 0x8298 | `B3-qr5cejk6c7c/S23-a4ejr5-k7c` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 66 |
| 0x7004 | `B3-qr4k5ce6c/S23-a4er5-k6aek` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 28 |
| 0x0b00 | `B3-qr5ce6c7e/S23-a4ejr5` | STABLE | 0 | 0.00 (0.00–0.06) | – | – | 0/64 | 64/64 | 12 |
