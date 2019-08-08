""" Generates plots as seen in the paper
"""

import numpy as np
import matplotlib.pyplot as plt

from helperplot import *
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


def plot1(n, save_file="temp"):
    """ Plots pairs saved in data/full_pairs_DS_n
    """
    with open(f"data/full_pairs_DS_{n}", "r") as f:
        eigvals = f.readlines()
    eigvals = [complex(s[:len(s)-1]) for s in eigvals]
    plt.figure(figsize=(7,3.5))
    plt.plot(np.real(eigvals), np.imag(eigvals), "o", markersize=.7,
               color="crimson", markeredgewidth=.03, markeredgecolor="black")
    plt_half_pm(n, linewidth=.5)
    plt.axis("off")
    plt.savefig(save_file, format="jpg", dpi=600, bbox_inches="tight")
    plt.close()

def plot2(save_file="temp"):
    """ Plots zoomed in figure of exception curve in DS_5
        as well as other pair that is near the boundary
        Exceptional pair: (123)(45) | (1523)
        Other pair: (1523) | (12453)
    """
    K1 = np.array([[0,1,0,0,0],[0,0,1,0,0],[1,0,0,0,0],[0,0,0,0,1],[0,0,0,1,0]]) # (123)(45)
    K2 = np.array([[0,0,0,0,1],[0,0,1,0,0],[1,0,0,0,0],[0,0,0,1,0],[0,1,0,0,0]]) # (1523)
    K3 = np.array([[0,1,0,0,0],[0,0,0,1,0],[1,0,0,0,0],[0,0,0,0,1],[0,0,1,0,0]]) # (12453)

    plt.figure(figsize=(3,2.5))

    eigvals = conv_comb_pairs([K1, K2], 0.0001)
    plt.plot(np.real(eigvals), np.imag(eigvals), "o", markersize=1.0,
               color="orange", markeredgewidth=.01, markeredgecolor="black")

    eigvals = conv_comb_pairs([K2, K3], 0.0001)
    plt.plot(np.real(eigvals), np.imag(eigvals), "o", markersize=1.0,
               color="orange", markeredgewidth=.01, markeredgecolor="black")

    plt_half_pm(5, linewidth=1)
    x_low, x_high = -.34, -.20
    y_low, y_high = .72, .80
    plt.xlim(x_low,x_high)
    plt.ylim(y_low,y_high)
    plt.locator_params(axis='y', nbins=4)
    plt.locator_params(axis='x', nbins=4)
    plt.savefig(save_file, format="jpg", dpi=600, bbox_inches="tight")
    plt.close()

def plot3(n, save_file="temp"):
    """ Plots perfect mirsky region of order n and Karpelevich region of order n-1
        Currently implemented for n=5,6, since Karpelevich region difficult to
        compute for general n.
    """
    plt.figure(figsize=(4,4))
    plt_pm(n)
    x_ranges, lines = pm_boundary(n)
    pts_lst = plt_K_n(n-1)

    # fill in region outside PM:
    print("K_n pts:", len(pts_lst))
    num_incr = 200
    ts = np.linspace(.8,1,num_incr)
    extr_lst = []
    for pt in pts_lst:
        extr_lst.extend((ts*pt))
    extr_lst = list(filter(lambda v: not in_region(v, x_ranges, lines) and 
                           v.real != 0, extr_lst))
    plt.plot(np.real(extr_lst), np.imag(extr_lst), "o",
             color="steelblue", markersize=.2)
    print("filtered points:", len(extr_lst))
    plt.xticks(np.arange(-1,1.01,.5))
    plt.yticks(np.arange(-1,1.01,.5))
    plt.savefig(save_file, format="jpg", dpi=600)
    plt.close()

def plot4(save_file="temp"):
    """ Plots zoomed in exceptional curve in DS_5 and triples that are outside of 
        the perfect mirsky region of order 5
    """
    import pickle

    plt.figure(figsize=(3,2.2))
    K1 = np.array([[0,1,0,0,0],[0,0,1,0,0],[1,0,0,0,0],[0,0,0,0,1],[0,0,0,1,0]]) # (123)(45)
    K2 = np.array([[0,0,0,0,1],[0,0,1,0,0],[1,0,0,0,0],[0,0,0,1,0],[0,1,0,0,0]]) # (1523)
    x_ranges, lines = pm_boundary(5)

    with open("data/filt_prec_trips_DS_5.pkl", "rb") as f:
        filt_trips = pickle.load(f)
    print("number of triples:",len(filt_trips))
    plt_pm(5)
    exc_pair = list(filter(lambda v: not in_region(v, x_ranges, lines), 
                           conv_comb_pairs([K1,K2], 0.0001)))
    plt.plot(np.real(exc_pair), np.imag(exc_pair), "o", 
                 color="orange", markersize=.4)
    plt.plot(np.real(filt_trips), np.imag(filt_trips), "o", 
                 color="purple", markersize=.25)

    x_low, x_high = -.325, -.25
    y_low, y_high = .74, .775
    plt.xlim(x_low,x_high)
    plt.ylim(y_low,y_high)
    plt.locator_params(axis='y', nbins=4)
    plt.locator_params(axis='x', nbins=4)
    plt.savefig(save_file, format="jpg", dpi=600)

if __name__ == "__main__":
    plot1(6)
    #plot2()
    #plot3(5)
    #plot4()
