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


def test_gpu_search_n4_n5():
    if not _torch_available():
        _skip("torch not available; skipping GPU search test.")
        return

    from permhull_gpu import gpu_search_exception

    res4 = gpu_search_exception(
        4,
        num_incr=400,
        max_pairs=None,
        pair_batch=32,
        device="cpu",
    )
    assert res4["found"] is False

    res5 = gpu_search_exception(
        5,
        num_incr=400,
        max_pairs=None,
        pair_batch=32,
        device="cpu",
    )
    assert res5["found"] is True


if __name__ == "__main__":
    test_gpu_search_n4_n5()
    print("GPU search tests pass")
