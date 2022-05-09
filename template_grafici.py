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

aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt('./data_w_errors/discreto/K_0.25/L_6_K_0.250000.dat', delimiter ="\t", unpack = True)

# RIORDINO GLI ARRAY
sorted_J = aux_J[aux_J.argsort()]
sorted_L = aux_L[aux_J.argsort()]
sorted_corr_len = aux_corr_len[aux_J.argsort()]
sorted_err_corr_len = aux_err_corr_len[aux_J.argsort()]
sorted_susc = aux_susc[aux_J.argsort()]
sorted_err_susc = aux_err_susc[aux_J.argsort()]
sorted_U = aux_U[aux_J.argsort()]
sorted_err_U = aux_err_U[aux_J.argsort()]
sorted_C = aux_C[aux_J.argsort()]
sorted_err_C = aux_err_C[aux_J.argsort()]
sorted_K2_g = aux_K2_g[aux_J.argsort()]
sorted_err_K2_g = aux_err_K2_g[aux_J.argsort()]
sorted_K2_sp = aux_K2_sp[aux_J.argsort()]
sorted_err_K2_sp = aux_err_K2_sp[aux_J.argsort()]

pl.figure(1)
pl.errorbar(sorted_J, sorted_K2_g, sorted_err_K2_g, ls = '', marker = m, color = c, fillstyle ='none', label = 'K2g, L =' + str(int(aux_L[0])))
c = next(c_cycle)
m = next(m_cycle)
pl.errorbar(sorted_J, sorted_K2_sp, sorted_err_K2_sp, ls = '', marker = m, color = c, fillstyle ='none', label = 'K2sp, L =' + str(int(aux_L[0])))
pl.xlabel(r'$J$')
pl.ylabel(r'$K_{2}$')
pl.legend()

aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt('./data_w_errors/discreto/K_0.25/L_8_K_0.250000.dat', delimiter ="\t", unpack = True)

# RIORDINO GLI ARRAY
sorted_J = aux_J[aux_J.argsort()]
sorted_L = aux_L[aux_J.argsort()]
sorted_corr_len = aux_corr_len[aux_J.argsort()]
sorted_err_corr_len = aux_err_corr_len[aux_J.argsort()]
sorted_susc = aux_susc[aux_J.argsort()]
sorted_err_susc = aux_err_susc[aux_J.argsort()]
sorted_U = aux_U[aux_J.argsort()]
sorted_err_U = aux_err_U[aux_J.argsort()]
sorted_C = aux_C[aux_J.argsort()]
sorted_err_C = aux_err_C[aux_J.argsort()]
sorted_K2_g = aux_K2_g[aux_J.argsort()]
sorted_err_K2_g = aux_err_K2_g[aux_J.argsort()]
sorted_K2_sp = aux_K2_sp[aux_J.argsort()]
sorted_err_K2_sp = aux_err_K2_sp[aux_J.argsort()]

pl.figure(2)
pl.errorbar(sorted_J, sorted_K2_g, sorted_err_K2_g, ls = '', marker = m, color = c, fillstyle ='none', label = 'K2g, L =' + str(int(aux_L[0])))
c = next(c_cycle)
m = next(m_cycle)
pl.errorbar(sorted_J, sorted_K2_sp, sorted_err_K2_sp, ls = '', marker = m, color = c, fillstyle ='none', label = 'K2sp, L =' + str(int(aux_L[0])))
pl.xlabel(r'$J$')
pl.ylabel(r'$K_{2}$')
pl.legend()

pl.show()
