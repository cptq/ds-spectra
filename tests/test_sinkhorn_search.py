import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src/")))


def _skip(msg):
    try:
        import pytest  # type: ignore
    except Exception:
        print(msg)
        return
    pytest.skip(msg)


def _torch_available():
    try:
        import torch  # type: ignore
    except Exception:
        return False
    return True


def test_sinkhorn_n4_n5():
    if not _torch_available():
        _skip("torch not available; skipping Sinkhorn search test.")
        return

    from sinkhorn_search import sinkhorn_pipeline

    common_kwargs = dict(
        batch_size=256,
        num_batches=20,
        top_k=8,
        init_mode="perm_mix",
        mix_k=2,
        mix_alpha=0.2,
        seed=0,
        device="cpu",
        table_bins=4096,
        opt_steps=0,
    )

    res4 = sinkhorn_pipeline(4, **common_kwargs)
    assert res4["best_exact_excess"] <= 1e-6

    res5 = sinkhorn_pipeline(5, **common_kwargs)
    assert res5["best_exact_excess"] > 1e-4


if __name__ == "__main__":
    test_sinkhorn_n4_n5()
    print("Sinkhorn search tests pass")
