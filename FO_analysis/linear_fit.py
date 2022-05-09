import numpy as np
import pylab as pl
import sys
import os
from scipy.optimize import curve_fit
from matplotlib.pyplot import cm
from scipy.odr import odrpack

from graph_style import *

def linear_1(x, a, b):
    return a + b*x

def linear_2(x, c, d):
    return 0.25*c*x + d

# CARICO I DATI (CHE SIANO C OD U CAMBIO SOLO IL LINK)
L, y1, dy1, y2, dy2 = np.genfromtxt('./prova_gauge.dat', delimiter = '\t', unpack = True)

# EVENTUALE TAGLIO DEI DATI
l_min = 9
l_max = 65

mask = ((L>=l_min) & (L<=l_max))
L_masked = L[mask]
y1_masked = y1[mask]
y2_masked = y2[mask]
dy1_masked = dy1[mask]
dy2_masked = dy2[mask]

# ANALISI PER JC CRITICO (USO x = 1/V, y1, dy1)
x_masked = 1.0/L_masked**3
x = 1.0/L**3

popt, pcov = curve_fit(linear_1, x_masked, y1_masked, sigma = dy1_masked, absolute_sigma = True)

# STAMPO I RISULTATI
chiq = np.sum( ((y1_masked - linear_1(x_masked, *popt)) / (dy1_masked))**2 )
ndof = len(x_masked) - len(popt)
red_chiq = chiq/ndof
perr = np.sqrt(np.diag(pcov))
print('jc = %f +- %f' % (popt[0], perr[0]))
print('m = %f +- %f' % (popt[1], perr[1]))
print('red_chiq = %f' % (red_chiq))

# PLOT
c_style, m_style = set_style()
pl.figure(1)
c = next(c_style)
m = next(m_style)

xx = np.linspace(np.min(x_masked), np.max(x_masked), 500)
y = linear_1(xx, *popt)
pl.errorbar(x_masked, y1_masked, dy1_masked, ls = '', marker = m, color = c, fillstyle = 'none', label = 'Misure')
pl.plot(xx, y, color = c, label = 'Fit')
pl.xlabel(r'$V^{-1}$')
pl.ylabel(r'$J_c^{C}$')
pl.legend()
# pl.savefig('../grafici/continuo/K_0.4/jc_vs_V.png')

# ANALISI PER L'ANDAMENTO DI C_MAX OD U
popt, pcov = curve_fit(linear_2, 1.0/x_masked, y2_masked, sigma = dy2_masked, absolute_sigma = True)

# STAMPO I RISULTATI
chiq = np.sum( ((y2_masked - linear_2(1.0/x_masked, *popt)) / (dy2_masked))**2 )
ndof = len(1.0/x_masked) - len(popt)
red_chiq = chiq/ndof
perr = np.sqrt(np.diag(pcov))
print('delta_h^2 = %f +- %f' % (popt[0], perr[0]))
print('m= %f +- %f' % (popt[1], perr[1]))
print('red_chiq = %f' % (red_chiq))

# PLOT
pl.figure(2)
xx = np.linspace(np.min(1.0/x_masked), np.max(1.0/x_masked), 500)
y = linear_2(xx, *popt)
pl.errorbar(1.0/x, y2, dy2, ls = '', marker = 'o', fillstyle = 'none')
pl.plot(xx, y)

# L, y, dy = np.genfromtxt('./intersection_N2_K0.4.dat', delimiter = '\t', unpack = True)
#
# x = 1/L**3
#
# popt, pcov = curve_fit(linear_1, x, y, sigma = dy, absolute_sigma = True)
#
# # STAMPO I RISULTATI
# chiq = np.sum( ((y - linear_1(x, *popt)) / (dy))**2 )
# ndof = len(x) - len(popt)
# red_chiq = chiq/ndof
# perr = np.sqrt(np.diag(pcov))
# print('jc = %f +- %f' % (popt[0], perr[0]))
# print('m = %f +- %f' % (popt[1], perr[1]))
# print('red_chiq = %f' % (red_chiq))
#
# # PLOT
# c_style, m_style = set_style()
# pl.figure(2)
# c = next(c_style)
# m = next(m_style)
#
# xx = np.linspace(np.min(x), np.max(x), 500)
# yy = linear_1(xx, *popt)
# pl.errorbar(x, y, dy, ls = '', marker = m, color = c, fillstyle = 'none', label = 'Misure')
# pl.plot(xx, yy, color = c, label = 'Fit')
# pl.xlabel(r'$V^{-1}$')
# pl.ylabel(r'$J_c^{prova}$')
# pl.legend()
# pl.savefig('../grafici/continuo/K_0.4/jc_U_intersect_V.png')

pl.show()
