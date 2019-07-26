import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))
from permhull import *
from evenpermhull import *

def test_perms():
    from math import factorial
    # tests that generated groups are correct size
    for n in range(2,7):
        assert len(list(symmetric_group(n))) == factorial(n)
        assert len(list(alternating_group(n))) == factorial(n)//2

def test_cycle_types():
    assert len(list(even_cycle_types(6))) == 6
    assert len(list(even_cycle_types(7))) == 8
    num_partitions = [0,1,2,3,5,7,11,15,22,30,42,
                      56,77,101,135,176]
    num_even_cycle_types = [0,1,1,2,3,4,6,8]
    for n in range(1,15):
        cycles = cycle_types(n)
        assert len(list(cycles)) == num_partitions[n]
        assert all(len(lst)==n for lst in cycles)
        for i in range(len(list(cycles))):
            for j in range(i+1, len(list(cycles))):
                assert cycles[i] != cycles[j]
        even_cycles = even_cycle_types(n)        
        for i in range(len(list(even_cycles))):
            for j in range(i+1, len(list(even_cycles))):
                assert not np.all(even_cycles[i] != even_cycles[j])
        assert all(len(lst)==n for lst in even_cycle_types(n))

def test_all_comb_coeff():
    from scipy.special import binom
    # checks that all tuples of convex coefficients sum to 1
    eps = 1e-14
    assert all(abs(sum(c)-1)< eps for c in all_comb_coeff(2,1,.049))         
    assert all(abs(sum(c)-1)< eps for c in all_comb_coeff(2,1,.1))         
    assert all(abs(sum(c)-1)< eps for c in all_comb_coeff(2,1,.12))         
    assert all(abs(sum(c)-1)< eps for c in all_comb_coeff(3,1,.1))         
    assert all(abs(sum(c)-1)< eps for c in all_comb_coeff(3,1,.12))         
    assert all(abs(sum(c)-1)< eps for c in all_comb_coeff(4,1,.1))         
    assert all(abs(sum(c)-1)< eps for c in all_comb_coeff(4,1,.12))         

def test_in_region():

    # PM_3
    x_ranges, lines = pm_boundary(3)
    assert in_region((.49 + .0j), x_ranges, lines)
    assert in_region((.49 + .2j), x_ranges, lines)
    assert not in_region((.49 - .5j), x_ranges, lines)
    assert not in_region((-.51 + .49j), x_ranges, lines)

    # PM_4
    x_ranges, lines = pm_boundary(4)
    assert in_region((.49 + .0j), x_ranges, lines)
    assert in_region((.49 + .2j), x_ranges, lines)
    assert not in_region((.49 - .52j), x_ranges, lines)
    assert in_region((-.51 + .49j), x_ranges, lines)


if __name__ == "__main__":
    test_perms()
    test_cycle_types()
    test_all_comb_coeff()
    test_in_region()
    print("All tests pass")


