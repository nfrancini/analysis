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

pl.xlabel(r'$R_\xi$')
pl.ylabel(r'$U$')

aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt('./data_w_errors/discreto_0.6pi/K_0.275/L_8_K_0.275000.dat', delimiter ="\t", unpack = True)

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

pl1 = pl.errorbar(sorted_corr_len/sorted_L, sorted_U, sorted_err_U, sorted_err_corr_len/sorted_L, ls='', color = c, marker = m, fillstyle = 'none')

# pl.legend(title = '$\kappa = 0.4$, discreto')

aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt('./data_w_errors/discreto_0.6pi/K_0.275/L_10_K_0.275000.dat', delimiter ="\t", unpack = True)
c = next(c_cycle)
m = next(m_cycle)

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

pl2 = pl.errorbar(sorted_corr_len/sorted_L, sorted_U, sorted_err_U, sorted_err_corr_len/sorted_L, ls='', color = c, marker = m, fillstyle = 'none')


aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt('./data_w_errors/discreto_0.6pi/K_0.275/L_12_K_0.275000.dat', delimiter ="\t", unpack = True)
c = next(c_cycle)
m = next(m_cycle)

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

pl3 = pl.errorbar(sorted_corr_len/sorted_L, sorted_U, sorted_err_U, sorted_err_corr_len/sorted_L, ls='', color = c, marker = m, fillstyle = 'none')

x, y, dx, dy = np.genfromtxt('./data_O4/O4_L16.dat', unpack = True)

pl4 = pl.errorbar(x, y, dy, dx, ls = '', color = 'green', marker = 'v', fillstyle = 'none')

x, y, dx, dy = np.genfromtxt('./data_O4/O4_L32.dat', unpack = True)

pl5 = pl.errorbar(x, y, dy, dx, ls = '', color = 'orange', marker = '<', fillstyle = 'none')

x, y, dx, dy = np.genfromtxt('./data_O4/O4_L64.dat', unpack = True)

pl6 = pl.errorbar(x, y, dy, dx, ls = '', color = 'purple', marker = '>', fillstyle = 'none')

title_proxy = Rectangle((0,0), 0, 0, color='w')
pl.legend([title_proxy, pl1, pl2, pl3, title_proxy, pl4, pl5, pl6],
           ["$\kappa=0.4$, d", "L=8", "L=16", "L=32", "$\kappa=\infty$, c", "L=16", "L=32", "L=64"])
# pl.savefig('./grafici/discreto/K_0.4/U_vs_xi_O4.png')
pl.show()
