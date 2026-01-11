"""Modal entrypoint for GPU search testing.

Run:
  modal run modal_gpu_search.py
"""
import modal

app = modal.App("ds-spectra-gpu-search")

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
def gpu_search(n=4, num_incr=100, max_pairs=200, batch_size=64):
    import sys

    sys.path.append("/root/src")
    from permhull_gpu import gpu_search_exception

    if max_pairs is not None and max_pairs <= 0:
        max_pairs = None
    if batch_size is not None and batch_size <= 0:
        batch_size = 64
    return gpu_search_exception(
        n,
        num_incr=num_incr,
        max_pairs=max_pairs,
        batch_size=batch_size,
    )


@app.local_entrypoint()
def main(n: int = 4, num_incr: int = 100, max_pairs: int = 200, batch_size: int = 64):
    if max_pairs is not None and max_pairs <= 0:
        max_pairs = None
    if batch_size is not None and batch_size <= 0:
        batch_size = 64
    result = gpu_search.remote(
        n=n,
        num_incr=num_incr,
        max_pairs=max_pairs,
        batch_size=batch_size,
    )
    if result is None:
        print(result)
        return
    pairs_str = result.get("pairs_checked_str", f"{result.get('pairs_checked', 0):,}")
    elapsed = result.get("elapsed_s")
    if elapsed is not None:
        print(f"pairs_checked: {pairs_str}")
        print(f"elapsed_s: {elapsed:.3f}")
    print(result)
