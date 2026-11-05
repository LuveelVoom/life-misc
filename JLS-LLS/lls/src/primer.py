from typing import Callable
from src.literal_manipulation import negate
from collections.abc import Generator, Iterable

literal = str  # JLS-LLS: `type` statements need Python 3.12; macOS has 3.9
clause = list[literal]

def is_power(n):
    return n != 0 and (n & (n - 1)) == 0

def make_prime(n, fn):
    t = []
    for i in range(1 << n):
        if not fn([(i >> j & 1) == 1 for j in range(n)]):
            t.append([i, 0, 0])
    for i, c in enumerate(t):
        for j, d in enumerate(t):
            diff = c[0] ^ d[0]
            if (is_power(diff)): c[2] |= diff
    for v in range(n):
        s = []
        j = 0
        for i, c in enumerate(t):
            if not c[2]: yield c[0], c[1]
            else:
                if c[2] >> (v+1):
                    s.append(c)
                if not c[0] >> v & 1 and c[2] >> v & 1:
                    while t[j][0] != c[0] | 1 << v or t[j][1] != c[1]:
                        j += 1
                    s.append([c[0], c[1] | 1 << v, c[2] & t[j][2] & ~(1 << v)])
        t = s
    for c in t:
        if not c[2]: yield c[0], c[1]

def instantiate_prime(var, enc):
    for a, b in enc:
        clause = [negate(x, a >> i & 1) for i, x in enumerate(var) if not b >> i & 1]
        yield clause
