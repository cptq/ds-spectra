"""Faster search utilities for Perfect-Mirsky exceptions (pairs only).

This keeps the older permhull.py intact and focuses on a faster
parallel_search_exception implementation.
"""
import itertools
import multiprocessing as mp

import numpy as np

_T_GRID = None
_IN_RAD = None
_X_LO = None
_X_HI = None
_M = None
_B = None
_IDENTITY = None


def _poly_boundary_arrays(k):
    pts = np.arange(0, k // 2, dtype=float)
    angles1 = 2 * np.pi * pts / k
    angles2 = 2 * np.pi * (pts + 1) / k
    x1 = np.cos(angles1)
    y1 = np.sin(angles1)
    x2 = np.cos(angles2)
    y2 = np.sin(angles2)
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    return x2, x1, m, b


def _pm_boundary_arrays(n):
    assert n >= 3
    x_lo = []
    x_hi = []
    m = []
    b = []
    for k in range(3, n + 1):
        poly_x_lo, poly_x_hi, poly_m, poly_b = _poly_boundary_arrays(k)
        x_lo.append(poly_x_lo)
        x_hi.append(poly_x_hi)
        m.append(poly_m)
        b.append(poly_b)
    return (
        np.concatenate(x_lo),
        np.concatenate(x_hi),
        np.concatenate(m),
        np.concatenate(b),
    )


def _in_region_mask(vals, x_lo, x_hi, m, b, eps=1e-14):
    if vals.size == 0:
        return np.zeros(0, dtype=bool)
    x = vals.real[:, None]
    y = np.abs(vals.imag)[:, None]
    in_x = (x >= (x_lo - eps)) & (x <= (x_hi + eps))
    below = y <= (m * x + b + eps)
    return np.any(in_x & below, axis=1)


def _perm_to_mat(perm):
    perm = np.asarray(perm, dtype=int)
    global _IDENTITY
    if _IDENTITY is None or _IDENTITY.shape[0] != len(perm):
        _IDENTITY = np.eye(len(perm), dtype=float)
    return _IDENTITY[perm - 1]


def symmetric_group(n):
    for perm in itertools.permutations(range(1, n + 1)):
        yield _perm_to_mat(perm)


def accel_asc(n):
    a = [0 for _ in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]


def cycle_types(n):
    partitions = accel_asc(n)
    types = []
    for parts in partitions:
        cyc_type = []
        parts = reversed(parts)
        i = 0
        for cyc_len in parts:
            cyc = [i + k + 1 for k in range(1, cyc_len)] + [i + 1]
            cyc_type.extend(cyc)
            i += cyc_len
        types.append(cyc_type)
    for perm in types:
        yield _perm_to_mat(perm)


def _pair_exception_eigval(C, P):
    diff = C - P
    mats = P[None, :, :] + _T_GRID[:, None, None] * diff[None, :, :]
    try:
        vals = np.linalg.eigvals(mats).reshape(-1)
        mask = (vals.imag > 0) & (vals.real != 0) & (np.abs(vals) > _IN_RAD)
        if not np.any(mask):
            return None
        cand = vals[mask]
        inside = _in_region_mask(cand, _X_LO, _X_HI, _M, _B)
        if np.all(inside):
            return None
        idx = np.where(~inside)[0][0]
        return cand[idx]
    except Exception:
        for t in _T_GRID:
            vals = np.linalg.eigvals(P + t * diff)
            mask = (vals.imag > 0) & (vals.real != 0) & (np.abs(vals) > _IN_RAD)
            if not np.any(mask):
                continue
            cand = vals[mask]
            inside = _in_region_mask(cand, _X_LO, _X_HI, _M, _B)
            if np.all(inside):
                continue
            idx = np.where(~inside)[0][0]
            return cand[idx]
    return None


def _worker_cycle_type(C):
    for P in symmetric_group(_IDENTITY.shape[0]):
        eigval = _pair_exception_eigval(C, P)
        if eigval is not None:
            return eigval, C, P
    return None


def _init_worker(n, num_incr):
    global _T_GRID, _IN_RAD, _X_LO, _X_HI, _M, _B, _IDENTITY
    _T_GRID = np.linspace(0.0, 1.0, num=num_incr, dtype=float)
    _IN_RAD = np.cos(np.pi / n)
    _X_LO, _X_HI, _M, _B = _pm_boundary_arrays(n)
    _IDENTITY = np.eye(n, dtype=float)


def parallel_search_exception(n, num_incr=10, processes=None, chunksize=1):
    """Search for exceptions using pairs, parallelized by cycle type."""
    if processes is None:
        processes = max(mp.cpu_count() - 1, 1)
    global _IDENTITY
    _IDENTITY = np.eye(n, dtype=float)
    cycle_type_mats = list(cycle_types(n))
    pool = mp.Pool(
        processes=processes,
        initializer=_init_worker,
        initargs=(n, num_incr),
    )
    try:
        for result in pool.imap_unordered(_worker_cycle_type, cycle_type_mats, chunksize=chunksize):
            if result is None:
                continue
            eigval, C, P = result
            print(f"EXCEPTION FOUND FOR n = {n}", "\n")
            print("eigenvalue:", eigval, "\n")
            print(C, "\n")
            print(P, "\n")
            pool.terminate()
            return result
    finally:
        pool.close()
        pool.join()
    print(f"No exception: n = {n}")
    return None
