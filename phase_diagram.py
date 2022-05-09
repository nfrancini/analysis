import numpy as np
import pylab as pl
import sys
import os
from scipy.optimize import curve_fit
from matplotlib.pyplot import cm
from scipy.odr import odrpack
from itertools import cycle
import natsort
from tqdm import tqdm
from matplotlib.patches import Rectangle

from graph_style import *
from corr_len import *
from susc import *
from FO_intersection import *

c_cycle, m_cycle = set_style()
c = next(c_cycle)
m = next(m_cycle)

k = np.array([0.04, 0.0759, 0.295, 0.4, 0.076])
J = np.array([0.602, 1.0, 0.0, 0.2343, 1.2])

line = pl.plot(k, J, marker = m, color = c, fillstyle='full', ls ='')[0]
line.set_clip_on(False)
pl.ylim((0,1.2))

# IDENTIFICO I VARI PUNTI
shift = 0.005
pl.text(0.04+shift, 0.602, 'FO')
pl.text(0.0759+shift, 1.0, 'XY')
pl.text(0.295+shift, 0.0+2*shift, 'XY')
pl.text(0.4-10*shift, 0.2343+6*shift, 'FO+O(4)')
pl.text(0.076+shift, 1.2+shift, 'XY')

pl.show()
