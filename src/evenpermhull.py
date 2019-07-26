""" Main scripts to compute EDS_n, the eigenvalues achievable by matrices in the convex hull of the even permutation matrices.
"""
import itertools
import csv
from math import factorial

import numpy as np
import multiprocessing as mp

from permhull import symmetric_group, cycle_types, all_comb_coeff

def alternating_group(n):
    """ generator for alternating group"""
    for p in filter(lambda p: np.linalg.det(p)>0, symmetric_group(n)):
        yield p

def even_cycle_types(n):
    """ list of one representative matrix for each cycle type
        for permutations of size n
    """
    return [p for p in cycle_types(n) if 
                 np.linalg.det(p)>0]

def compute_all_combs_EDS_n(n, num_incr=30, 
                save_file="all_comb_EDS", save_chunk=10000):
    """ n is size of problem
        num_incr is mesh size
        all convex combinations!
        very expensive
    """
    eigvals = []
    mesh = np.linspace(0,1,num=num_incr)
    assert n <= 5, "Very expensive to compute all combinations!!" 
    perms = list(alternating_group(n))
    count=0
    with open(save_file+f"_{n}", "w") as f:
        for coeffs in all_comb_coeff(factorial(n)//2, 1, 1/num_incr):
            mat = np.zeros((n,n))
            for i in range(factorial(n)//2):
                mat += coeffs[i]*perms[i]
            vals = np.linalg.eigvals(mat)
            eigvals.extend(vals)
            count+=len(vals)
            if len(eigvals) >= save_chunk:
                    f.writelines(str(val)+"\n" for val in eigvals)
                    count = 0
                    del eigvals
                    eigvals=[]
        f.writelines(str(val)+"\n" for val in eigvals)

def compute_eff_EDS_n(n, num_incr=100, save_file="eff_EDS",
                     save_chunk=10000):
    """ n is size of problem
        num_incr is mesh size
    """
    eigvals = []
    in_rad = np.cos(np.pi/(n-1)) if n%2 == 0 else np.cos(np.pi/n)
    count = 0
    with open(save_file+f"_{n}", "w") as f:
        for C in even_cycle_types(n):
            for P in alternating_group(n):
                vals = [val for t in np.linspace(0,1,num=num_incr)
                        for val in np.linalg.eigvals(t*C + (1-t)*P)]
                eigvals.extend(filter(
                        lambda val: val.imag > 0 and val.real != 0 and
                                abs(val)>in_rad, vals))
                count += num_incr*n
                if count >= save_chunk:
                    f.writelines(str(val)+"\n" for val in eigvals)
                    count = 0
                    del eigvals
                    eigvals=[]
        f.writelines(str(val)+"\n" for val in eigvals)

def compute_pairs_EDS_n(n, num_incr=100, save_file="all_EDS",
                     save_chunk=10000):
    """ n is size of problem
        num_incr is mesh size
    """
    eigvals = []
    for C in even_cycle_types(n):
        for P in alternating_group(n):
            vals = [val for t in np.linspace(0,1,num=num_incr)
                    for val in np.linalg.eigvals(t*C + (1-t)*P)]
            eigvals.extend(vals)
    return eigvals
