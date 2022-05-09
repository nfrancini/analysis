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

aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt("./data_w_errors/discreto/J_0/L_8_J_0.000000.dat", delimiter ="\t", unpack = True)

# RIORDINO GLI ARRAY
sorted_K = aux_K[aux_K.argsort()]
sorted_J = aux_J[aux_K.argsort()]
sorted_L = aux_L[aux_K.argsort()]
sorted_corr_len = aux_corr_len[aux_K.argsort()]
sorted_err_corr_len = aux_err_corr_len[aux_K.argsort()]
sorted_susc = aux_susc[aux_K.argsort()]
sorted_err_susc = aux_err_susc[aux_K.argsort()]
sorted_U = aux_U[aux_K.argsort()]
sorted_err_U = aux_err_U[aux_K.argsort()]
sorted_C = aux_C[aux_K.argsort()]
sorted_err_C = aux_err_C[aux_K.argsort()]

c_cycle, m_cycle = set_style()
c = next(c_cycle)
m = next(m_cycle)

pl.figure(1)
pl.errorbar(sorted_K, sorted_C, sorted_err_C, ls = '', fillstyle='none', color=c, marker=m, label='L=8')
pl.xlabel(r'$\kappa$')
pl.ylabel(r'$C$')
pl.legend(loc=2)
pl.savefig('./grafici/discreto/J_0/C_vs_k.png')

pl.show()
