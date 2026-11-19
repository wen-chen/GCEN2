"""R's default random number generator, for the few places where GCEN2
reproduces an R function that draws random numbers.

``RRandom(seed)`` is R's state after ``set.seed(seed)`` with the default
kinds (``RNGkind("Mersenne-Twister", "Inversion", "Rejection")``, R >= 3.6):
:meth:`unif_rand` is ``runif(1)`` and :meth:`sample` is
``sample.int(n, size)`` without replacement (both of R's algorithms: the
partial shuffle, and the hash version R uses for n > 1e7 and size <= n / 2).
"""

from __future__ import annotations

import math

import numpy as np

_N, _M = 624, 397
# numpy scalars: with numpy 1.x, uint64 op Python int is promoted to float64
_UPPER, _LOWER = np.uint64(0x80000000), np.uint64(0x7FFFFFFF)
_MATRIX_A = np.uint64(0x9908B0DF)
_I2_32M1 = 2.328306437080797e-10  # 1 / (2^32 - 1)


class RRandom:
    """R's Mersenne-Twister after ``set.seed(seed)``."""

    def __init__(self, seed: int):
        s = int(seed) & 0xFFFFFFFF
        for _ in range(50):  # initial scrambling
            s = (69069 * s + 1) & 0xFFFFFFFF
        seeds = []
        for _ in range(_N + 1):
            s = (69069 * s + 1) & 0xFFFFFFFF
            seeds.append(s)
        # dummy[0] is mti, set to N by FixupSeeds; mt = dummy[1:]
        self._mt = np.array(seeds[1:], dtype=np.uint32)
        self._mti = _N
        self._buf = np.empty(0)
        self._pos = 0

    def _generate(self) -> None:
        mt = self._mt.astype(np.uint64)
        mag = np.array([0, _MATRIX_A], dtype=np.uint64)
        for kk in range(_N - _M):
            y = (mt[kk] & _UPPER) | (mt[kk + 1] & _LOWER)
            mt[kk] = mt[kk + _M] ^ (y >> np.uint64(1)) ^ mag[int(y & np.uint64(1))]
        for kk in range(_N - _M, _N - 1):
            y = (mt[kk] & _UPPER) | (mt[kk + 1] & _LOWER)
            mt[kk] = mt[kk + (_M - _N)] ^ (y >> np.uint64(1)) ^ mag[int(y & np.uint64(1))]
        y = (mt[_N - 1] & _UPPER) | (mt[0] & _LOWER)
        mt[_N - 1] = mt[_M - 1] ^ (y >> np.uint64(1)) ^ mag[int(y & np.uint64(1))]
        self._mt = mt.astype(np.uint32)
        y = mt.copy()
        y ^= y >> np.uint64(11)
        y ^= (y << np.uint64(7)) & np.uint64(0x9D2C5680)
        y ^= (y << np.uint64(15)) & np.uint64(0xEFC60000)
        y ^= y >> np.uint64(18)
        y &= np.uint64(0xFFFFFFFF)
        u = y.astype(np.float64) * 2.3283064365386963e-10
        # fixup: inside (0, 1)
        u[u <= 0.0] = 0.5 * _I2_32M1
        u[1.0 - u <= 0.0] = 1.0 - 0.5 * _I2_32M1
        self._buf = u
        self._pos = 0

    def unif_rand(self) -> float:
        """One draw of ``runif(1)``."""
        if self._pos >= self._buf.size:
            self._generate()
        v = self._buf[self._pos]
        self._pos += 1
        return float(v)

    def _rbits(self, bits: int) -> int:
        v = 0
        for _ in range(0, bits + 1, 16):
            v = 65536 * v + int(math.floor(self.unif_rand() * 65536))
        return v & ((1 << bits) - 1)

    def unif_index(self, dn: float) -> int:
        """R_unif_index(dn): uniform on 0 .. dn - 1 (rejection sampling)."""
        if dn <= 0:
            return 0
        bits = int(math.ceil(math.log2(dn)))
        while True:
            dv = self._rbits(bits)
            if dn > dv:
                return dv

    def sample(self, n: int, size: int) -> np.ndarray:
        """``sample.int(n, size)`` (without replacement): 1-based values."""
        n, size = int(n), int(size)
        if size > n:
            raise ValueError("cannot take a sample larger than the population")
        out = np.empty(size, dtype=np.int64)
        if n > 1e7 and size <= n / 2:  # R's hash version (sample2)
            seen = set()
            i = 0
            while i < size:
                v = self.unif_index(n) + 1
                if v not in seen:
                    seen.add(v)
                    out[i] = v
                    i += 1
            return out
        x = np.arange(n, dtype=np.int64)
        m = n
        for i in range(size):
            j = self.unif_index(m)
            out[i] = x[j] + 1
            m -= 1
            x[j] = x[m]
        return out
