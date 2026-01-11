"""Compare old vs new parallel search runtime."""
import argparse
from time import perf_counter

from permhull import parallel_search_exception as old_search
from permhull_fast import parallel_search_exception as new_search


def _run(label, func, n, num_incr):
    print(f"{label} (n={n}, num_incr={num_incr})")
    start = perf_counter()
    func(n, num_incr=num_incr)
    elapsed = perf_counter() - start
    print(f"{label} elapsed: {elapsed:.3f}s\n")


def main():
    parser = argparse.ArgumentParser(description="Compare permhull search runtimes.")
    parser.add_argument("-n", type=int, default=4, help="problem size (default: 4)")
    parser.add_argument("--num-incr", type=int, default=10, help="number of t samples (default: 10)")
    args = parser.parse_args()

    _run("old", old_search, args.n, args.num_incr)
    _run("new", new_search, args.n, args.num_incr)


if __name__ == "__main__":
    main()
