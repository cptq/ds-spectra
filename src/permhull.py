""" Main scripts to compute DS_n and search for exceptions to Perfect-Mirsky.
"""
import itertools
import csv
from math import factorial

import numpy as np
import multiprocessing as mp

def poly_boundary(k):
    """ represents boundary of kth roots of unity polygon
        with x_ranges between vertices of k-gon
        and equations of lines for sides of k-gon
    """
    x_ranges = []
    lines = []
    for pt in range(0, k//2):
        x1, x2 = np.cos(pt*2*np.pi / k), np.cos((pt+1)*2*np.pi/k)
        y1, y2 = np.sin(pt*2*np.pi / k), np.sin((pt+1)*2*np.pi/k)
        m = (y2 - y1)/(x2-x1)
        b = y1-m*x1
        x_ranges.append((x2,x1))
        lines.append((m,b))
    return x_ranges, lines

def pm_boundary(n):
    """ represents boundary of perfect mirsky region
        of order n by boundaries of k-gons for k <= n
        only use for n >= 3
    """
    assert n >= 3
    x_ranges = []
    lines = []
    for k in range(3, n+1):
        poly_x, poly_lines = poly_boundary(k)
        x_ranges.extend(poly_x)
        lines.extend(poly_lines)
    return x_ranges, lines

def in_region(val, x_ranges, lines, eps=1e-14) -> bool:
    """ checks if val is in region with
        boundary determined by x_ranges, lines
        x_ranges and lines are returned from
        pm_boundary(n)
    """
    for i in range(len(x_ranges)):
        x1, x2 = x_ranges[i]
        if val.real >= x1-eps and val.real <= x2+eps:
            m, b = lines[i]
            if abs(val.imag) < m*val.real + b + eps:
                return True
    return False

def get_out_region(vals, n):
    """ returns points in vals that are not
        in perfect mirksy region of order n
    """
    x_ranges, lines = pm_boundary(n)
    return list(filter(lambda v: not in_region(v, x_ranges, lines),
                       vals))


def symmetric_group(n):
    """ symmetric group of size n
        e.g. the 3-cycle (123) is represented as [2,3,1]
    """
    for perm in itertools.permutations([i for i in range(1,n+1)]):
        yield perm_to_mat(perm)

def rev_symmetric_group(n):
    """ symmetric group of size n
        reversed order
    """
    for perm in itertools.permutations(reversed([i for i in range(1,n+1)])):
        yield perm_to_mat(perm)
        
def perm_to_mat(perm):
    """ converts permutation to matrix
        perm is permutation of [1,...,n]
    """
    n = len(perm)
    P = np.zeros((n,n))
    for i in range(n):
        P[i, perm[i]-1] = 1
    return P

def cycle_types(n):
    """ list of one representative matrix for each cycle type
        for permutations of size n
    """
    partitions = accel_asc(n)
    types = []
    for parts in partitions:
        cyc_type = []
        parts = reversed(parts)
        i = 0
        for cyc_len in parts:
            cyc = [i+k+1 for k in range(1,cyc_len)] + [i+1]
            cyc_type.extend(cyc)
            i+=cyc_len
        types.append(cyc_type)
    return map(perm_to_mat, types)

def accel_asc(n):
    """ compute partitions of size n
        algorithm by jerome kelleher
    """
    a = [0 for i in range(n + 1)]
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

def conv_comb_pairs(As, incr=0.01):
    """ eigenvalues achieved by convex
        combinations of pairs of matrices in As
    """
    n = len(As)
    eigvals = []
    for p in itertools.combinations(As, 2):
        for t in np.linspace(0,1,int(1/incr)):
            vals = np.linalg.eigvals(t*p[0] + (1-t)*p[1])
            eigvals.extend(vals)
    return eigvals

def all_comb_coeff(k, target, incr):
    """ returns list of all combinations
        of k nonnegative coefficients that
        add up to target and scale by incr
    """
    if k == 1:
        return [target]
    if target <= incr:
        lst = [[0]*k for _ in range(k)]
        for i in range(k):
            lst[i][i] = target
        return lst
    return [np.append([t], lst) for 
            t in np.append(np.arange(0,target,incr),target) for 
            lst in all_comb_coeff(k-1, (target-t), incr)]

def conv_comb_eigv_k(As, num_incr=50):
    """ computes eigenvalues of convex combinations of
        k-tuples of matrices in As,
        where As contains k matrices
    """
    k = len(As)
    eigvals = []
    coeffs = all_comb_coeff(k, 1, 1/num_incr)
    for c in coeffs:
        mat = sum(c[i]*As[i] for i in range(len(c)))
        vals = np.linalg.eigvals(mat)
        eigvals.extend(vals)
    return eigvals

def compute_all_combs_DS_n(n, num_incr=30):
    """ n is size of problem
        num_incr is mesh size
        all convex combinations!
        very expensive
    """
    eigvals = []
    mesh = np.linspace(0,1,num=num_incr)
    assert n <= 5, "Very expensive to compute all combinations!!"
    perms = list(symmetric_group(n))
    for coeffs in all_comb_coeff(factorial(n), 1, 1/num_incr):
        mat = np.zeros((n,n))
        for i in range(factorial(n)):
            mat += coeffs[i]*perms[i]
        vals = np.linalg.eigvals(mat)
        eigvals.extend(vals)
    return eigvals

def compute_pairs_DS_n(n, num_incr=100, save_file="pairs_DS",
                    save_chunk=10000):
    """ 
        omits inradius, real eigenvalues, and lower half plane
        n is size of problem
        num_incr is mesh size
    """
    eigvals = []
    in_rad = np.cos(np.pi/n)
    count=0
    with open(save_file+f"_{n}", "w") as f:
        for C in cycle_types(n):
            for P in symmetric_group(n):
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

def compute_trips_DS_n(n, num_incr=20,
            save_file="trips_DS", save_chunk=100000):
 
    """ 
        omits inradius, real eigenvalues, and lower half plane
        n is size of problem
        num_incr is mesh size
    """
    eigvals = []
    coeffs = all_comb_coeff(3, 1, 1/num_incr)
    in_rad = np.cos(np.pi/n)
    count = 0
    with open(save_file+f"_{n}", "w") as f:
        for C in cycle_types(n):
            for P1 in symmetric_group(n):
                for P2 in symmetric_group(n):
                    mat = np.zeros((n,n))
                    for c in coeffs:
                        vals = np.linalg.eigvals(
                            c[0]*C + c[1]*P1 + c[2]*P2)
                        eigvals.extend(filter(
                        lambda val: val.imag > 0 and val.real != 0 and
                                abs(val)>in_rad, vals))

                        count += len(vals)
                    if count >= save_chunk:
                        f.writelines(str(val)+"\n" for val in eigvals)
                        count = 0
                        del eigvals
                        eigvals=[]
        f.writelines(str(val)+"\n" for val in eigvals)

def search_exception_mids(n, num_incr=1000):
    """ search for exception to perfect mirsky for
        given n, checking midpoints of pairs
    """
    eigvals = []
    in_rad = np.cos(np.pi/n)
    x_ranges, lines = pm_boundary(n)
    for P in reversed(list(cycle_types(n))):
        for Q in rev_symmetric_group(n):
            vals = np.linalg.eigvals(1/2*P + 1/2*Q)
            eigvals = list(filter(
            lambda val: val.imag > 0 and val.real != 0 and
              abs(val)>in_rad, vals))
            if len(eigvals) > 0:
                exc_lst = list(filter(lambda val: not in_region(val, x_ranges,lines),
                             eigvals))
                if len(exc_lst) > 0:
                    print(f"EXCEPTION FOUND FOR n = {n}", "\n")
                    print("eigenvalue:", exc_lst[0], "\n")
                    print(P, "\n")
                    print(Q, "\n")
                    return
    print(f"No exception: n = {n}")

def make_inequiv_pairs(n):
    """ Convert from gap permutation pairs output to txt format
        gap format is assumed to be stored in Pairs-n-1-1.txt
    """
    path = "data10/"
    fname = f"Pairs-{n}-1-1.txt"
    with open(f"{path}{fname}", "r") as f:
        lst = f.readlines()
    inequiv_parser_pairs(lst, f"{path}New_{fname}", n)

def inequiv_parser_pairs(lst, out_file, n, save_chunk=100000):
    """ parse gap txt output in lst
        make new format
    """
    pairs_lst = []
    count = 0
    with open(out_file, "w") as f:
        i = 0
        while i < len(lst):
            r = lst[i]
            row_count = 0
            # cycles of first perm in pair
            middle = r[1:-1]
            left_paren = middle.find("[")
            right_paren = middle.find("],")
            first_perm = middle[left_paren:right_paren+1]
            if first_perm.find(".") < 0:
                first_perm = first_perm.split()
                first_lst = [int(i.strip(",")) for i in first_perm if i.strip(",").isdigit()]
            else:
                first_lst = [i+1 for i in range(n)] # is identity

            # cycles of second perm in pair
            sec_half = middle[right_paren+1:]
            sec_left_paren = sec_half.find("[")
            if sec_left_paren < 0: # on a half line
                i += 1
                r = lst[i]
                middle = r[1:-1]
                sec_left_paren = middle.find("[")
                sec_half = middle
            sec_right_paren = sec_half.rfind("]")
            sec_perm = sec_half[sec_left_paren:sec_right_paren+1]
            if sec_perm.find(".") < 0:
                sec_perm = sec_perm.split()
                sec_lst = [int(i.strip(",")) for i in sec_perm if i.strip(",").isdigit()]
            else:
                sec_lst = [i+1 for i in range(n)] # is identity

            pairs_lst.append((first_lst, sec_lst))
            count += 1
            if count >= save_chunk:
                # save memory by writing in chunks
                f.writelines((str(pair)[1:-1]+"\n" for pair in pairs_lst))
                del pairs_lst
                pairs_lst = []
                count = 0
            i+=1

        f.writelines((str(pair)[1:-1]+"\n" for pair in pairs_lst))
    del pairs_lst
    print("done saving")

def read_inequiv_pairs(n):
    """ read New_Pairs-n-1-1.txt and return pairs of permutations
    """
    with open(f"data10/New_Pairs-{n}-1-1.txt", "r") as f:
        lst = f.readlines()
    perms = []
    for r in lst:
        middle = r[0:-1]
        left_paren = middle.find("[")
        right_paren = middle.find("],")
        first_perm = middle[left_paren+1:right_paren]
        first_perm = first_perm.replace(",", " ").split()
        first_lst = [int(i) for i in first_perm]

        sec_half = middle[right_paren+1:]
        sec_left_paren = sec_half.find("[")
        sec_right_paren = sec_half.rfind("]")
        sec_perm = sec_half[sec_left_paren+1:sec_right_paren]
        sec_perm = sec_perm.replace(",", " ").split()
        sec_lst = [int(i) for i in sec_perm]
        perms.append((first_lst, sec_lst))
    return perms

def compute_inequiv_pairs_DS_n(n, num_incr=100, save_file="new_pairs_DS",
                    save_chunk=10000):
    """ Computes eigenvalues of inequivalent pairs
        n is size of problem
        num_incr is mesh size
    """
    eigvals = []
    in_rad = np.cos(np.pi/n)
    in_rad = -1 # save all
    count=0
    with open(save_file+f"_{n}", "w") as f:
        writer = csv.writer(f)
        pairs = map(lambda t: (perm_to_mat(t[0]),
                               perm_to_mat(t[1])),
                               read_inequiv_pairs(n))
        for P, Q in pairs:
            vals = [val for t in np.linspace(0,1,num=num_incr)
                    for val in np.linalg.eigvals(t*P + (1-t)*Q)]
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

def inequiv_search_exception(n, num_incr=1000):
    """ search for exception to perfect mirsky for
        given n, checking pairs
        uses inequivalent pairs 
    """
    eigvals = []
    in_rad = np.cos(np.pi/n)
    count=0
    x_ranges, lines = pm_boundary(n)
    pairs = map(lambda t: (perm_to_mat(t[0]),
                               perm_to_mat(t[1])),
                               read_inequiv_pairs(n))
    for P, Q in pairs:
        vals = [val for t in np.linspace(0,1,num=num_incr)
          for val in np.linalg.eigvals(t*P + (1-t)*Q)]
        eigvals = list(filter(
            lambda val: val.imag > 0 and val.real != 0 and
            abs(val)>in_rad, vals))
        if len(eigvals) > 0:
            exc_lst = list(filter(lambda val: not in_region(val, x_ranges,lines),
                      eigvals))
            if len(exc_lst) > 0:
                print(f"EXCEPTION FOUND FOR n = {n}", "\n")
                print("eigenvalue:", exc_lst[0], "\n")
                print(P, "\n")
                print(Q, "\n")
                return
    print(f"No exception: n = {n}")

def search_exception(n, num_incr=1000):
    """ search for exception to perfect mirsky for
        given n, checking pairs
        checks in reversed order to check pairs with permutations
        of less fixed points first
    """
    eigvals = []
    in_rad = np.cos(np.pi/n)
    x_ranges, lines = pm_boundary(n)
    for P in reversed(list(cycle_types(n))):
        for Q in rev_symmetric_group(n):
            vals = [val for t in np.linspace(0,1,num=num_incr)
              for val in np.linalg.eigvals(t*P + (1-t)*Q)]
            eigvals = list(filter(
            lambda val: val.imag > 0 and val.real != 0 and
              abs(val)>in_rad, vals))
            if len(eigvals) > 0:
                exc_lst = list(filter(lambda val: not in_region(val, x_ranges,lines),
                             eigvals))
                if len(exc_lst) > 0:
                    print(f"EXCEPTION FOUND FOR n = {n}", "\n")
                    print("eigenvalue:", exc_lst[0], "\n")
                    print(P, "\n")
                    print(Q, "\n")
                    return
    print(f"No exception: n = {n}")

def search_trips_exception(n, num_incr=100):
    """ search for exception to perfect mirsky for
        given n, checking triples
    """
    eigvals = []
    in_rad = np.cos(np.pi/n)
    count=0
    coeffs = all_comb_coeff(3, 1, 1/num_incr)
    x_ranges, lines = pm_boundary(n)
    for P in reversed(list(cycle_types(n))):
        for Q in rev_symmetric_group(n):
            for R in rev_symmetric_group(n):
                    for c in coeffs:
                        vals = np.linalg.eigvals(
                            c[0]*P + c[1]*Q + c[2]*R)
                        eigvals.extend(filter(
                        lambda val: val.imag > 0 and val.real != 0 and
                                abs(val)>in_rad, vals))
                    if len(eigvals) > 0:
                        exc_lst = list(filter(lambda val: not in_region(val, x_ranges,lines), eigvals))
                        if len(exc_lst) > 0:
                            print(f"EXCEPTION FOUND FOR n = {n}", "\n")
                            print("eigenvalue:", exc_lst[0], "\n")
                            print(P, "\n")
                            print(Q, "\n")
                            print(R, "\n")
                            return
    print(f"No exception: n = {n}")

def worker_pairs(As, num_incr=10):
    """ returns eigenvalues of convex combinations of pairs for As"""
    return conv_comb_pairs(As, 1/num_incr)

def parallel_search_exception(n):
    """ searches for exception to the perfect mirsky conjecture
        of size n, using multiprocessing for parallelization
        uses one less than total cpus available
        n is size of problem
        num_incr is controlled by worker_pairs !!
    """
    eigvals = []
    in_rad = np.cos(np.pi/n)
    pool = mp.Pool(max(mp.cpu_count()-1, 1)) 
    x_ranges, lines = pm_boundary(n)
    for vals in pool.imap(worker_pairs, itertools.product(cycle_types(n), symmetric_group(n))):
            #((C,P) for C in cycle_types(n) for P in symmetric_group(n))):
        eigvals.extend(filter(lambda val: val.imag > 0 and val.real != 0 and
            abs(val)>in_rad, vals))

        if len(eigvals) > 0:
                exc_lst = list(filter(lambda val: not in_region(val, x_ranges,lines), eigvals))
                if len(exc_lst) > 0:
                    print(f"EXCEPTION FOUND FOR n = {n}", "\n")
                    print("eigenvalue:", exc_lst[0], "\n")
                    return
    print(f"No exception: n = {n}")

