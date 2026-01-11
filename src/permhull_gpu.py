"""GPU-backed search utilities for Perfect-Mirsky exceptions (pairs only).

This is a basic GPU implementation meant for testing with Modal.
"""
import itertools
import time

import numpy as np

try:
    import torch
except Exception as exc:  # pragma: no cover - only raised when GPU deps missing
    raise ImportError(
        "torch is required for GPU search. Install torch with CUDA or run via Modal."
    ) from exc


def _poly_boundary_arrays(k):
    """Return boundary arrays (x_lo, x_hi, m, b) for a k-gon."""
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
    """Return concatenated boundary arrays for the Perfect-Mirsky region."""
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
    """Vectorized region membership for complex vals (True if inside)."""
    if vals.size == 0:
        return np.zeros(0, dtype=bool)
    x = vals.real[:, None]
    y = np.abs(vals.imag)[:, None]
    in_x = (x >= (x_lo - eps)) & (x <= (x_hi + eps))
    below = y <= (m * x + b + eps)
    return np.any(in_x & below, axis=1)


def _perm_to_mat(perm):
    """Convert a 1-based permutation into its permutation matrix."""
    perm = np.asarray(perm, dtype=int)
    n = len(perm)
    return np.eye(n, dtype=float)[perm - 1]


def symmetric_group(n):
    """Yield permutation matrices for S_n in lexicographic order."""
    for perm in itertools.permutations(range(1, n + 1)):
        yield _perm_to_mat(perm)


def accel_asc(n):
    """Yield integer partitions of n (Kelleher's accel_asc algorithm)."""
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
    """Yield one permutation matrix representative per cycle type in S_n."""
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


def _eigvals_batch(P, diff, t_grid):
    mats = P.unsqueeze(0) + t_grid[:, None, None] * diff
    return torch.linalg.eigvals(mats).reshape(-1)


def gpu_search_exception(n, num_incr=10, max_pairs=None):
    """Search for exceptions using GPU eigenvalues over pair convex combos.

    Args:
        n: size of the problem.
        num_incr: number of t samples along each pair segment.
        max_pairs: optional cap for the number of pairs to check.

    Returns:
        dict with search stats and, if found, eigenvalue (real/imag tuple) and pair data.
    """
    start = time.perf_counter()
    in_rad = np.cos(np.pi / n)
    x_lo, x_hi, m, b = _pm_boundary_arrays(n)
    device = torch.device("cuda")
    t_grid = torch.linspace(0.0, 1.0, steps=num_incr, device=device, dtype=torch.float64)
    in_rad_t = torch.tensor(in_rad, device=device, dtype=torch.float64)

    checked = 0
    for C in cycle_types(n):
        C_gpu = torch.tensor(C, device=device, dtype=torch.float64)
        for P in symmetric_group(n):
            if max_pairs is not None and checked >= max_pairs:
                elapsed = time.perf_counter() - start
                return {
                    "n": int(n),
                    "pairs_checked": int(checked),
                    "pairs_checked_str": f"{checked:,}",
                    "elapsed_s": float(elapsed),
                    "found": False,
                }
            checked += 1
            P_gpu = torch.tensor(P, device=device, dtype=torch.float64)
            diff = C_gpu - P_gpu
            vals = _eigvals_batch(P_gpu, diff, t_grid)
            mask = (vals.imag > 0) & (vals.real != 0) & (vals.abs() > in_rad_t)
            if torch.any(mask):
                cand = vals[mask].cpu().numpy()
                inside = _in_region_mask(cand, x_lo, x_hi, m, b)
                if not np.all(inside):
                    idx = np.where(~inside)[0][0]
                    val = cand[idx]
                    elapsed = time.perf_counter() - start
                    return {
                        "n": int(n),
                        "eigenvalue": (float(val.real), float(val.imag)),
                        "C": C.tolist() if hasattr(C, "tolist") else C,
                        "P": P.tolist() if hasattr(P, "tolist") else P,
                        "pairs_checked": int(checked),
                        "pairs_checked_str": f"{checked:,}",
                        "elapsed_s": float(elapsed),
                        "found": True,
                    }
    elapsed = time.perf_counter() - start
    return {
        "n": int(n),
        "pairs_checked": int(checked),
        "pairs_checked_str": f"{checked:,}",
        "elapsed_s": float(elapsed),
        "found": False,
    }
