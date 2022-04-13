import numpy as np
import pylab as pl
import matplotlib
from cycler import cycler


def set_style():
    color_cycler = iter(['black', 'red', 'blue', 'green'])
    marker_cycler = iter(['s', 'o', '^', 'v'])

    pl.rc('text', usetex = True)
    pl.rc('font', family = 'serif', size = 15)
    pl.rc('figure', autolayout = True)
    pl.rc('lines', linestyle = '-', linewidth = 1)
    pl.rc('markers', fillstyle = 'none')
    pl.rc('axes', labelsize = 'large')
    pl.rc('xtick', top = True, bottom = True, direction = 'in')
    pl.rc('xtick.minor', visible = True)
    pl.rc('ytick', left = True, right =  True, direction = 'in')
    pl.rc('ytick.minor', visible = True)
    pl.rc('legend', fancybox = False, edgecolor = 'black')
    pl.rc('errorbar', capsize = 2 )

    return color_cycler, marker_cycler
