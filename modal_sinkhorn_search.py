"""Modal entrypoint for Sinkhorn-based continuous search.

Run:
  modal run modal_sinkhorn_search.py
"""
import modal

app = modal.App("ds-spectra-sinkhorn-search")

image = (
    modal.Image.from_registry("nvidia/cuda:12.1.1-devel-ubuntu22.04", add_python="3.9")
    .pip_install(
        "numpy",
        "torch",
        index_url="https://download.pytorch.org/whl/cu121",
    )
    .add_local_dir("src", remote_path="/root/src", copy=True)
)


@app.function(
    gpu="A10G",
    image=image,
    timeout=60 * 60,
)
def sinkhorn_search(
    n=12,
    batch_size=256,
    num_batches=50,
    top_k=8,
    sinkhorn_iters=10,
    score_temp=0.02,
    table_bins=4096,
    init_mode="gaussian",
    mix_k=2,
    mix_alpha=0.3,
    max_perm_cache=200000,
    opt_steps=0,
    opt_lr=0.05,
    entropy_weight=0.0,
    seed=-1,
):
    import sys

    sys.path.append("/root/src")
    from sinkhorn_search import sinkhorn_pipeline

    seed = None if seed is None or seed < 0 else seed
    return sinkhorn_pipeline(
        n=n,
        batch_size=batch_size,
        num_batches=num_batches,
        top_k=top_k,
        sinkhorn_iters=sinkhorn_iters,
        score_temp=score_temp,
        table_bins=table_bins,
        init_mode=init_mode,
        mix_k=mix_k,
        mix_alpha=mix_alpha,
        max_perm_cache=max_perm_cache,
        opt_steps=opt_steps,
        opt_lr=opt_lr,
        entropy_weight=entropy_weight,
        seed=seed,
    )


@app.local_entrypoint()
def main(
    n: int = 12,
    batch_size: int = 256,
    num_batches: int = 50,
    top_k: int = 8,
    sinkhorn_iters: int = 10,
    score_temp: float = 0.02,
    table_bins: int = 4096,
    init_mode: str = "gaussian",
    mix_k: int = 2,
    mix_alpha: float = 0.3,
    max_perm_cache: int = 200000,
    opt_steps: int = 0,
    opt_lr: float = 0.05,
    entropy_weight: float = 0.0,
    seed: int = -1,
    score_threshold: float = 1e-4,
):
    result = sinkhorn_search.remote(
        n=n,
        batch_size=batch_size,
        num_batches=num_batches,
        top_k=top_k,
        sinkhorn_iters=sinkhorn_iters,
        score_temp=score_temp,
        table_bins=table_bins,
        init_mode=init_mode,
        mix_k=mix_k,
        mix_alpha=mix_alpha,
        max_perm_cache=max_perm_cache,
        opt_steps=opt_steps,
        opt_lr=opt_lr,
        entropy_weight=entropy_weight,
        seed=seed,
    )
    top = result.get("top", [])
    best = top[0] if top else None
    best_score = best.get("score") if best else None
    is_counterexample = best_score is not None and best_score > score_threshold
    result["best_score"] = best_score
    result["counterexample_found"] = bool(is_counterexample)
    if is_counterexample:
        print("!!! COUNTEREXAMPLE FOUND: YES !!!")
        print(f"best_score: {best_score}")
        print(f"eigenvalue: {best.get('eigenvalue')}")
        print("!!! COUNTEREXAMPLE FOUND: YES !!!")
    else:
        print("COUNTEREXAMPLE FOUND: NO")
        print(f"best_score: {best_score}")
        print(f"threshold: {score_threshold}")
    print(result)
