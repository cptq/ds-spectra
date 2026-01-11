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


def _perm_to_mat_zero_based(perm, n):
    """Convert a 0-based permutation into its permutation matrix."""
    perm = np.asarray(perm, dtype=int)
    return np.eye(n, dtype=float)[perm]


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


def _cycle_type_perms(n):
    """Return 0-based permutation arrays for each cycle type."""
    partitions = accel_asc(n)
    perms = []
    for parts in partitions:
        cyc_type = []
        parts = reversed(parts)
        i = 0
        for cyc_len in parts:
            cyc = [i + k + 1 for k in range(1, cyc_len)] + [i + 1]
            cyc_type.extend(cyc)
            i += cyc_len
        perms.append(np.asarray(cyc_type, dtype=np.int64) - 1)
    return perms


def _in_region_mask_torch(vals, x_lo, x_hi, m, b, eps=1e-14):
    """Torch region membership for complex vals (True if inside)."""
    if vals.numel() == 0:
        return torch.zeros(0, dtype=torch.bool, device=vals.device)
    x = vals.real[:, None]
    y = vals.imag.abs()[:, None]
    in_x = (x >= (x_lo - eps)) & (x <= (x_hi + eps))
    below = y <= (m * x + b + eps)
    return (in_x & below).any(dim=1)


def _eigvals_batch(P_batch, diff_batch, t_grid, fallback_counts=None):
    mats = P_batch[:, None, :, :] + t_grid[None, :, None, None] * diff_batch[:, None, :, :]
    n = P_batch.shape[-1]
    mats = mats.reshape(-1, n, n)
    try:
        vals = torch.linalg.eigvals(mats)
    except RuntimeError:
        # Some batches can be ill-conditioned or nearly defective; GPU eigvals may fail.
        # Fall back to CPU and, if necessary, NumPy to keep the search running.
        if fallback_counts is not None:
            fallback_counts["gpu_to_cpu"] += 1
        try:
            vals = torch.linalg.eigvals(mats.cpu()).to(mats.device)
        except RuntimeError:
            if fallback_counts is not None:
                fallback_counts["cpu_to_numpy"] += 1
            mats_np = mats.cpu().numpy()
            vals_np = np.stack([np.linalg.eigvals(m) for m in mats_np], axis=0)
            vals = torch.from_numpy(vals_np).to(mats.device)
    return vals.reshape(P_batch.shape[0], t_grid.shape[0], n)


def gpu_search_exception(n, num_incr=10, max_pairs=None, batch_size=64, device=None):
    """Search for exceptions using GPU eigenvalues over pair convex combos.

    Args:
        n: size of the problem.
        num_incr: number of t samples along each pair segment.
        max_pairs: optional cap for the number of pairs to check.
        batch_size: number of permutation pairs to batch per eigensolve.
        device: torch device string or object (defaults to CUDA if available).

    Returns:
        dict with search stats and, if found, eigenvalue (real/imag tuple) and pair data.
    """
    batch_size = max(int(batch_size), 1)
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    device = torch.device(device)
    start = time.perf_counter()
    in_rad = np.cos(np.pi / n)
    x_lo, x_hi, m, b = _pm_boundary_arrays(n)
    dtype = torch.float64
    t_grid = torch.linspace(0.0, 1.0, steps=num_incr, device=device, dtype=dtype)
    in_rad_t = torch.tensor(in_rad, device=device, dtype=dtype)
    x_lo_t = torch.tensor(x_lo, device=device, dtype=dtype)
    x_hi_t = torch.tensor(x_hi, device=device, dtype=dtype)
    m_t = torch.tensor(m, device=device, dtype=dtype)
    b_t = torch.tensor(b, device=device, dtype=dtype)

    perms_cpu = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
    perms_gpu = torch.from_numpy(perms_cpu).to(device)
    eye = torch.eye(n, device=device, dtype=dtype)
    cycle_perms = _cycle_type_perms(n)

    checked = 0
    fallback_counts = {"gpu_to_cpu": 0, "cpu_to_numpy": 0}
    with torch.no_grad():
        for C_perm in cycle_perms:
            C_gpu = eye[torch.tensor(C_perm, device=device)]
            total_pairs = perms_cpu.shape[0]
            step = batch_size
            for batch_start in range(0, total_pairs, step):
                remaining = None if max_pairs is None else max_pairs - checked
                if remaining is not None and remaining <= 0:
                    elapsed = time.perf_counter() - start
                    return {
                        "n": int(n),
                        "pairs_checked": int(checked),
                        "pairs_checked_str": f"{checked:,}",
                        "elapsed_s": float(elapsed),
                        "found": False,
                    }
                current_batch = min(step, total_pairs - batch_start)
                if remaining is not None:
                    current_batch = min(current_batch, remaining)
                if current_batch <= 0:
                    elapsed = time.perf_counter() - start
                    return {
                        "n": int(n),
                        "pairs_checked": int(checked),
                        "pairs_checked_str": f"{checked:,}",
                        "elapsed_s": float(elapsed),
                        "found": False,
                    }
                batch_end = batch_start + current_batch
                P_batch = eye[perms_gpu[batch_start:batch_end]]
                diff = C_gpu.unsqueeze(0) - P_batch
                vals = _eigvals_batch(P_batch, diff, t_grid, fallback_counts).reshape(-1)
                mask = (vals.imag > 0) & (vals.real != 0) & (vals.abs() > in_rad_t)
                if torch.any(mask):
                    cand_idx = torch.nonzero(mask, as_tuple=False).flatten()
                    cand = vals[cand_idx]
                    inside = _in_region_mask_torch(cand, x_lo_t, x_hi_t, m_t, b_t)
                    if not torch.all(inside):
                        bad_rel = torch.nonzero(~inside, as_tuple=False)[0].item()
                        bad_idx = cand_idx[bad_rel].item()
                        pair_idx = bad_idx // (num_incr * n)
                        val = vals[bad_idx].cpu().numpy()
                        checked = checked + int(pair_idx) + 1
                        P_perm = perms_cpu[batch_start + int(pair_idx)]
                        C_mat = _perm_to_mat_zero_based(C_perm, n)
                        P_mat = _perm_to_mat_zero_based(P_perm, n)
                        elapsed = time.perf_counter() - start
                        return {
                            "n": int(n),
                            "eigenvalue": (float(val.real), float(val.imag)),
                            "C": C_mat.tolist(),
                            "P": P_mat.tolist(),
                            "pairs_checked": int(checked),
                            "pairs_checked_str": f"{checked:,}",
                            "elapsed_s": float(elapsed),
                            "fallback_counts": fallback_counts,
                            "found": True,
                        }
                checked += current_batch
    elapsed = time.perf_counter() - start
    return {
        "n": int(n),
        "pairs_checked": int(checked),
        "pairs_checked_str": f"{checked:,}",
        "elapsed_s": float(elapsed),
        "fallback_counts": fallback_counts,
        "found": False,
    }
