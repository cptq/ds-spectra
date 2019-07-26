""" Helper functions for plotting figures
"""
import numpy as np
import matplotlib.pyplot as plt

from permhull import *

# set matplotlib parameters
params = {
       'axes.labelsize': 8,
       'font.size': 8,
       'legend.fontsize': 8,
       'xtick.labelsize': 8,
       'ytick.labelsize': 8,
       'text.usetex': False,
       'figure.figsize': [4, 4],
        "font.family": "serif",
    }
plt.rcParams.update(params)

def plt_half_pm(n, linewidth=1.5):
    """ plots subset of perfect-mirsky region of size n
        in upper half plane
    """
    for k in range(2,n+1):
        rts = [np.exp(2*l * np.pi * 1j/k) for l in range(k)]
        rts.append(rts[0])
        plt.plot(np.real(rts), np.imag(rts), linewidth=linewidth, 
                  color="black")
        plt.ylim(-0.001,1)
        
def plt_pm(n, linewidth=1.5):
    """ plots perfect-mirsky region of size n
    """
    for k in range(2,n+1):
        rts = [np.exp(2*l * np.pi * 1j/k) for l in range(k)]
        rts.append(rts[0])
        plt.plot(np.real(rts), np.imag(rts), linewidth=linewidth,
                 color="black") 

def plt_K_n(n, num_incr=250):
    """ Plots Karpelevich region of size n
    """
    if n == 4: return plt_K_4(num_incr=num_incr)
    elif n == 5: return plt_K_5(num_incr=num_incr)
    else: raise ValueError("Karpelevich region Unavailable for this n")

def plt_K_4(num_incr):
    """ Plots Karpelevich region of size 4
    """
    ts = np.linspace(0, 1, num_incr)
    pts_lst = []
    x_ranges, lines = pm_boundary(4)
    for t in ts:
        s = 1-t
        pts_lst.extend([t*np.exp(1j*np.pi/2) + s])
        pts_lst.extend([t*np.exp(-1j*np.pi/2) + s])
        
        pts = np.roots([1,0,0,-s,-t])
        pts = [pt for pt in pts if not in_region(pt, x_ranges, lines)]
        pts_lst.extend(pts)
        
        pts = np.roots([1,0,-2*t,-s**2,+t**2])
        pts = pts[pts != 1]
        pts = [pt for pt in pts if not in_region(pt, x_ranges, lines)]
        pts_lst.extend(pts)
    pts_lst = list(filter(lambda v: v.real < .98, pts_lst)) # numerical error
    return pts_lst

def plt_K_5(num_incr):
    """ Plots Karpelevich region of size 5
    """
    ts = np.linspace(0, 1, num_incr)
    pts_lst = []
    x_ranges, lines = pm_boundary(5)
    for t in ts:
        s = 1-t
        pts_lst.extend([t*np.exp(2j*np.pi/5) + s])
        pts_lst.extend([t*np.exp(-2j*np.pi/5) + s])
        pts = np.roots([1,0,0,0,-s,-t])
        pts = [pt for pt in pts if not in_region(pt, x_ranges, lines)]
        pts_lst.extend(pts)
        
        pts = np.roots([1,0,0,-s,-t])
        pts = pts[pts != 1]
        pts = [pt for pt in pts if not in_region(pt, x_ranges, lines)]
        pts_lst.extend(pts)
        
        pts = np.roots([1,0,0,-s,0,-t])
        pts = pts[pts != 1]
        pts = [pt for pt in pts if not in_region(pt, x_ranges, lines)]
        pts_lst.extend(pts)
        
        pts = np.roots([1,0,-2*t,0,t**2,-s**2])
        pts = pts[pts != 1]
        pts = [pt for pt in pts if not in_region(pt, x_ranges, lines)]
        pts_lst.extend(pts)
    pts_lst = list(filter(lambda v: v.real < .98, pts_lst)) # numerical error
    return pts_lst

