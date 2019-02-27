from cataloging_classes import *
import matplotlib.pyplot as plt
import matplotlib.style as style
import numpy as np
import scipy.stats as stats
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit

style.use('ggplot')


tableau20 = [(31, 119, 180), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)   

dd_1a = pickle.load(open(cwd + "/Pickles/dd_sne1a.p", "rb"))
dd_2 = pickle.load(open(cwd + "/Pickles/dd_sne2.p", "rb"))
dd_uk = pickle.load(open(cwd + "/Pickles/dd_sneuk.p", "rb"))



def histograph(data, xlab, ylab, title, error=None, nbins=50, filename=None, cum=False):
    plt.hist(data, bins=nbins, cumulative=cum, edgecolor='black')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def cumulative_graph(data_lists, types, xlab, ylab, title, nbins=100, filename=None, ks=False):
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    min_max = (min_list_of_lists(data_lists), max_list_of_lists(data_lists))
    colors = tableau20[:len(data_lists)]
    patches = []
    for l, c, sne in zip(data_lists, colors, types):
        n, bins, _ = plt.hist(l, nbins, min_max, cumulative=True, fill=False, linewidth=0, density=True)
        n = np.insert(n, 0, 0.)
        plt.plot(bins, n, color=c)
        patches.append(mpatches.Patch(color=c, label=sne))
    if len(patches) > 1:
        plt.legend(handles=patches)
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
    if ks:
        if len(data_lists) != 2:
            raise ValueError("To perform a 2-sample KS test there has to be exactly 2 samples!")
        print(stats.ks_2samp(data_lists[0], data_lists[1]))


def vs_graph(db, col1, col2, xlab='', ylab='', title='', filename=None, linear=True):
    """
    Takes a database and 2 columns and plots them.
    """
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    x = db[col1].values.tolist()
    print(len(x))
    y = db[col2].values.tolist()
    print(len(y))
    try:
        plt.scatter(x, y)
    except ValueError:
        x = [row[0] for row in x]
        plt.scatter(x, y)
    plt.show()





    