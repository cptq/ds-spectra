"""Sinkhorn-based continuous search for Perfect-Mirsky violations."""
import itertools
import math
import time

import numpy as np

try:
    import torch
except Exception as exc:  # pragma: no cover - only raised when GPU deps missing
    raise ImportError(
        "torch is required for Sinkhorn search. Install torch with CUDA or run via Modal."
    ) from exc

_PERM_CACHE = {}


def sinkhorn(logits, iters=10, eps=1e-8, temperature=1.0):
    scale = max(float(temperature), 1e-6)
    x = torch.exp(logits / scale)
    for _ in range(iters):
        x = x / (x.sum(dim=-1, keepdim=True) + eps)
        x = x / (x.sum(dim=-2, keepdim=True) + eps)
    return x


def pm_radius_table(n, bins=4096, device=None, dtype=torch.float32):
    if n < 3:
        raise ValueError("n must be >= 3")
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    device = torch.device(device)
    theta = torch.linspace(0.0, 2 * math.pi, bins + 1, device=device, dtype=dtype)[:-1]
    ks = torch.arange(3, n + 1, device=device, dtype=dtype)
    seg = 2 * math.pi / ks
    angle = torch.remainder(theta[:, None], seg[None, :])
    r_k = torch.cos(math.pi / ks) / torch.cos(angle - math.pi / ks)
    return r_k.max(dim=1).values


def pm_radius_from_table(theta, table):
    bins = table.shape[0]
    theta = torch.remainder(theta, 2 * math.pi)
    scaled = theta * (bins / (2 * math.pi))
    idx0 = torch.floor(scaled).long() % bins
    idx1 = (idx0 + 1) % bins
    t = scaled - idx0.to(scaled.dtype)
    return table[idx0] * (1 - t) + table[idx1] * t


def eigvals_batch(mats):
    try:
        return torch.linalg.eigvals(mats)
    except RuntimeError:
        try:
            return torch.linalg.eigvals(mats.cpu()).to(mats.device)
        except RuntimeError:
            mats_np = mats.detach().cpu().numpy()
            vals_np = np.stack([np.linalg.eigvals(m) for m in mats_np], axis=0)
    return torch.from_numpy(vals_np).to(mats.device)


def score_eigvals(eigvals, table, temp=0.02, exclude_tol=1e-6):
    theta = torch.atan2(eigvals.imag, eigvals.real)
    r = pm_radius_from_table(theta, table)
    excess = eigvals.abs() - r
    mask = torch.abs(eigvals - 1) < exclude_tol
    excess = excess.masked_fill(mask, -1e9)
    if temp is None or temp <= 0:
        score = excess.max(dim=-1).values
    else:
        score = temp * torch.logsumexp(excess / temp, dim=-1)
    return score, excess


def pm_radius_exact(theta, n):
    theta = np.remainder(theta, 2 * math.pi)
    r_max = np.zeros_like(theta, dtype=float)
    for k in range(3, n + 1):
        seg = 2 * math.pi / k
        angle = np.remainder(theta, seg)
        denom = np.cos(angle - (math.pi / k))
        r_k = (math.cos(math.pi / k) / denom)
        r_max = np.maximum(r_max, r_k)
    return r_max


def max_excess_exact(eigvals, n):
    theta = np.angle(eigvals)
    r = pm_radius_exact(theta, n)
    excess = np.abs(eigvals) - r
    idx = int(np.argmax(excess))
    return float(excess[idx]), eigvals[idx]


def perm_mats(n, device, dtype, max_perm_cache=200000):
    total = math.factorial(n)
    if total > max_perm_cache:
        return None
    key = (n, str(device), dtype)
    cached = _PERM_CACHE.get(key)
    if cached is not None:
        return cached
    perms = list(itertools.permutations(range(n)))
    perm_idx = torch.tensor(perms, device=device, dtype=torch.long)
    mats = torch.eye(n, device=device, dtype=dtype)[perm_idx]
    _PERM_CACHE[key] = mats
    return mats


def perm_mix_batch(perms, batch_size, mix_k=2, mix_alpha=0.3, generator=None):
    num_perms = perms.shape[0]
    device = perms.device
    idx = torch.randint(num_perms, (batch_size, mix_k), device=device, generator=generator)
    sel = perms[idx]
    weights = torch.rand((batch_size, mix_k), device=device, generator=generator)
    if mix_alpha is not None and mix_alpha > 0:
        weights = weights ** (1.0 / float(mix_alpha))
    weights = weights / weights.sum(dim=1, keepdim=True)
    mats = (sel * weights[:, :, None, None]).sum(dim=1)
    return mats


def summarize_candidates(logits, n, table, sinkhorn_iters=10, temperature=1.0, score_temp=0.02):
    mats = sinkhorn(logits, iters=sinkhorn_iters, temperature=temperature)
    eigvals = eigvals_batch(mats)
    score, excess = score_eigvals(eigvals, table, temp=score_temp)
    max_excess = excess.max(dim=-1).values
    best_idx = excess.argmax(dim=-1)
    rows = torch.arange(eigvals.shape[0], device=eigvals.device)
    best_vals = eigvals[rows, best_idx]
    return (
        score.detach().cpu(),
        max_excess.detach().cpu(),
        best_vals.detach().cpu(),
        mats.detach().cpu(),
    )


def stage1_random_search(
    n,
    batch_size=256,
    num_batches=50,
    top_k=8,
    sinkhorn_iters=10,
    temperature=1.0,
    score_temp=0.02,
    table_bins=4096,
    init_mode="gaussian",
    mix_k=2,
    mix_alpha=0.3,
    max_perm_cache=200000,
    device=None,
    dtype=torch.float32,
    seed=None,
):
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    device = torch.device(device)
    table = pm_radius_table(n, bins=table_bins, device=device, dtype=dtype)
    gen = None
    if seed is not None:
        gen = torch.Generator(device=device)
        gen.manual_seed(int(seed))

    best_scores = None
    best_logits = None
    perms = None
    use_perm_mix = init_mode in {"perm_mix", "hybrid"}
    if use_perm_mix:
        perms = perm_mats(n, device, dtype, max_perm_cache=max_perm_cache)
        if perms is None:
            use_perm_mix = False
    with torch.no_grad():
        for _ in range(num_batches):
            if use_perm_mix and (init_mode == "perm_mix" or (init_mode == "hybrid" and torch.rand(()) < 0.5)):
                mats = perm_mix_batch(
                    perms,
                    batch_size,
                    mix_k=mix_k,
                    mix_alpha=mix_alpha,
                    generator=gen,
                )
                eigvals = eigvals_batch(mats)
                scores, excess = score_eigvals(eigvals, table, temp=score_temp)
                rank_scores = excess.max(dim=-1).values
                logits = (mats + 1e-8).log()
            else:
                logits = torch.randn(
                    (batch_size, n, n),
                    device=device,
                    dtype=dtype,
                    generator=gen,
                )
                scores, rank_scores, _, _ = summarize_candidates(
                    logits,
                    n,
                    table,
                    sinkhorn_iters=sinkhorn_iters,
                    temperature=temperature,
                    score_temp=score_temp,
                )
            top = min(top_k, scores.shape[0])
            top_vals, top_idx = torch.topk(rank_scores, k=top)
            top_vals_cpu = top_vals.detach().cpu()
            logits_cpu = logits[top_idx].detach().cpu()
            if best_scores is None:
                best_scores = top_vals_cpu.clone()
                best_logits = logits_cpu.clone()
            else:
                combo_scores = torch.cat([best_scores, top_vals_cpu], dim=0)
                combo_logits = torch.cat([best_logits, logits_cpu], dim=0)
                keep = min(top_k, combo_scores.shape[0])
                new_scores, new_idx = torch.topk(combo_scores, k=keep)
                best_scores = new_scores
                best_logits = combo_logits[new_idx]
    return {
        "table": table,
        "best_scores": best_scores if best_scores is not None else torch.empty(0),
        "best_logits": best_logits if best_logits is not None else torch.empty(0),
    }


def refine_candidates(
    logits,
    n,
    table,
    steps=50,
    lr=0.05,
    sinkhorn_iters=10,
    temperature=1.0,
    temp_end=None,
    temp_anneal="linear",
    score_temp=0.02,
    entropy_weight=0.0,
):
    device = logits.device
    logits = logits.clone().detach().requires_grad_(True)
    opt = torch.optim.Adam([logits], lr=lr)
    eps = 1e-8
    if temp_end is None:
        temp_end = temperature
    if temp_anneal not in {"linear", "exp"}:
        raise ValueError("temp_anneal must be 'linear' or 'exp'")
    for step in range(steps):
        if steps <= 1 or temperature == temp_end:
            temp = temperature
        else:
            frac = step / (steps - 1)
            if temp_anneal == "exp":
                temp = temperature * ((temp_end / max(temperature, 1e-6)) ** frac)
            else:
                temp = temperature + (temp_end - temperature) * frac
        opt.zero_grad(set_to_none=True)
        mats = sinkhorn(logits, iters=sinkhorn_iters, temperature=temp)
        eigvals = eigvals_batch(mats)
        score, _ = score_eigvals(eigvals, table, temp=score_temp)
        loss = -score.mean()
        if entropy_weight:
            entropy = (mats * (mats + eps).log()).sum(dim=(-1, -2)).mean()
            loss -= entropy_weight * entropy
        loss.backward()
        opt.step()
    return logits.detach()


def sinkhorn_pipeline(
    n,
    batch_size=256,
    num_batches=50,
    top_k=8,
    sinkhorn_iters=10,
    temperature=1.0,
    temp_end=None,
    temp_anneal="linear",
    score_temp=0.02,
    table_bins=4096,
    init_mode="gaussian",
    mix_k=2,
    mix_alpha=0.3,
    max_perm_cache=200000,
    opt_steps=0,
    opt_lr=0.05,
    entropy_weight=0.0,
    device=None,
    dtype=torch.float32,
    seed=None,
):
    start = time.perf_counter()
    stage1 = stage1_random_search(
        n,
        batch_size=batch_size,
        num_batches=num_batches,
        top_k=top_k,
        sinkhorn_iters=sinkhorn_iters,
        temperature=temperature,
        score_temp=score_temp,
        table_bins=table_bins,
        init_mode=init_mode,
        mix_k=mix_k,
        mix_alpha=mix_alpha,
        max_perm_cache=max_perm_cache,
        device=device,
        dtype=dtype,
        seed=seed,
    )
    table = stage1["table"]
    best_logits = stage1["best_logits"]
    best_scores = stage1["best_scores"]
    if best_logits.numel() == 0:
        return {
            "elapsed_s": time.perf_counter() - start,
            "top": [],
        }

    device = table.device
    logits = best_logits.to(device)
    if opt_steps > 0:
        logits = refine_candidates(
            logits,
            n,
            table,
            steps=opt_steps,
            lr=opt_lr,
            sinkhorn_iters=sinkhorn_iters,
            temperature=temperature,
            temp_end=temp_end,
            temp_anneal=temp_anneal,
            score_temp=score_temp,
            entropy_weight=entropy_weight,
        )

    scores, max_excess, best_vals, mats = summarize_candidates(
        logits,
        n,
        table,
        sinkhorn_iters=sinkhorn_iters,
        temperature=temperature,
        score_temp=score_temp,
    )
    order = torch.argsort(max_excess, descending=True)
    top = []
    exact_scores = []
    exact_eigs = []
    mats_np = mats.detach().cpu().numpy()
    order_top = order[:top_k]
    for idx in order_top:
        exact_excess, exact_eigval = max_excess_exact(np.linalg.eigvals(mats_np[idx]), n)
        exact_scores.append(exact_excess)
        exact_eigs.append(exact_eigval)
    for i, idx in enumerate(order_top):
        eigval = best_vals[idx]
        top.append(
            {
                "score": float(scores[idx]),
                "max_excess": float(max_excess[idx]),
                "exact_excess": float(exact_scores[i]),
                "eigenvalue": (float(eigval.real), float(eigval.imag)),
            }
        )

    best_exact_idx = int(np.argmax(exact_scores)) if exact_scores else 0
    best_mat = mats_np[order_top[best_exact_idx]]
    exact_excess = exact_scores[best_exact_idx] if exact_scores else 0.0
    exact_eigval = exact_eigs[best_exact_idx] if exact_eigs else 0.0 + 0.0j
    elapsed = time.perf_counter() - start
    return {
        "elapsed_s": float(elapsed),
        "top": top,
        "best_matrix": best_mat.tolist(),
        "best_exact_excess": exact_excess,
        "best_exact_eigenvalue": (float(exact_eigval.real), float(exact_eigval.imag)),
        "stage1_top_scores": best_scores.tolist(),
    }
