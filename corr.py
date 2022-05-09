import pylab as pl
import numpy as np
import sys
import os
from tqdm import tqdm
from scipy.optimize import curve_fit

from graph_style import *

def correlazioni(x, cutoff):
    # ARRAY PER I RISULTATI
    kk = np.array([])
    corr = np.array([])
    int_tau = np.array([])

    ave = np.mean(x)

    # LOOP PER LE CORRELAZIONI
    tau = 0
    for k in tqdm(range(0, cutoff)):
        sum2 = 0
        C = 0
        for i in range(0, len(x) - k):
            sum2 = sum2 + (x[i] - ave)*(x[i+k] - ave)
        C = sum2/(len(x) - k)
        tau = tau + C
        kk = np.append(kk, k)
        corr = np.append(corr, C)
        int_tau = np.append(int_tau,tau)
    norm = corr[0]
    corr = corr/norm
    int_tau = int_tau/norm - 1

    return kk, corr, int_tau

def gauss(x, mu, sigma, A):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
    return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

ctrl_obs = True
ctrl_bimodal = False

fname = sys.argv[1]

L, V, N, J, K = np.genfromtxt(fname, dtype = "double", delimiter = "\t", unpack = True, max_rows = 1)
ene_sp, ene_g, ene_dens, susc, G_pm, mu2 = np.genfromtxt(fname, dtype = "double", delimiter = "\t", unpack = True, skip_header = 2)

# SCARTO LA TESTA DEL FILE SCARTO IL PRIMO DECIMO
skip = len(ene_dens)//10
ene_dens = ene_dens[skip:]
ene_sp = ene_sp[skip:]
ene_g = ene_g[skip:]
susc = susc[skip:]
G_pm = G_pm[skip:]
mu2 = mu2[skip:]

cutoff = int(L**2)


k_ene_sp, corr_ene_sp, tau_ene_sp = correlazioni(ene_sp, cutoff)

if (ctrl_obs == True):
    k_ene_g, corr_ene_g, tau_ene_g = correlazioni(ene_g, cutoff)
    k_ene_dens, corr_ene_dens, tau_ene_dens = correlazioni(ene_dens, cutoff)
    k_susc, corr_susc, tau_susc = correlazioni(susc, cutoff)
    k_Gpm, corr_Gpm, tau_Gpm = correlazioni(G_pm, cutoff)
    k_mu2, corr_mu2, tau_mu2 = correlazioni(mu2, cutoff)

pl.figure(1)
pl.suptitle('ene_spin')
pl.subplot(311)
pl.scatter(k_ene_sp, corr_ene_sp)
pl.ylabel('C')
pl.subplot(312)
pl.scatter(k_ene_sp, tau_ene_sp)
pl.ylabel('tau')
pl.subplot(313)
y, x, _ = pl.hist(ene_sp, bins = int(np.sqrt(len(ene_sp))), weights = np.array([1.0/len(ene_sp)]*len(ene_sp)))
pl.ylabel('Freq')

if(ctrl_bimodal == True):
    y, x, _ =pl.hist(ene_g, bins = int(np.sqrt(len(ene_g))), weights = np.array([1.0/len(ene_g)]*len(ene_g)))
    x = (x[1:]+x[:-1])/2
    init = (0.24, 0.04, 0.004, 0.38, 0.04, 0.003)
    popt, pcov = curve_fit(bimodal, x, y, init)
    perr = np.sqrt(np.diag(pcov))
    xx = np.linspace(x.min(), x.max(), 500)
    pl.plot(xx, bimodal(xx, *popt), color = 'black')
    print('mu1 = %f +- %f' % (popt[0], perr[0]))
    print('mu2 = %f +- %f' % (popt[3], perr[3]))
    delta_h = np.abs(popt[0]-popt[3])
    err_delta_h = np.sqrt(perr[0]**2 + perr[3]**2)
    print('delta_h = %f +- %f' % (delta_h, err_delta_h))

    c, m = set_style()
    pl.figure(7)
    pl.xlabel(r'$\epsilon_{g}$')
    pl.ylabel(r'Frequenza')
    pl.plot(xx, bimodal(xx, *popt), color = 'black', label='Fit')
    y, x, _ = pl.hist(ene_g, bins = int(np.sqrt(len(ene_g))), weights = np.array([1.0/len(ene_g)]*len(ene_g)), label = 'L=' + str(int(L)), color = 'gray')
    pl.legend()
    pl.savefig('./grafici/discreto/K_0.275/L_16_ene_g_hist.png')


if (ctrl_obs == True):
    pl.figure(2)
    pl.suptitle('ene_gauge')
    pl.subplot(311)
    pl.scatter(k_ene_g, corr_ene_g)
    pl.ylabel('C')
    pl.subplot(312)
    pl.scatter(k_ene_g, tau_ene_g)
    pl.ylabel('tau')
    pl.xlabel('k')
    pl.subplot(313)
    pl.hist(ene_g, bins = int(np.sqrt(len(ene_g))), weights = np.array([1.0/len(ene_g)]*len(ene_g)))
    pl.ylabel('Freq')

    pl.figure(3)
    pl.suptitle('ene_dens')
    pl.subplot(311)
    pl.scatter(k_ene_dens, corr_ene_dens)
    pl.ylabel('C')
    pl.subplot(312)
    pl.scatter(k_ene_dens, tau_ene_dens)
    pl.ylabel('tau')
    pl.xlabel('k')
    pl.subplot(313)
    pl.hist(ene_dens, bins = int(np.sqrt(len(ene_dens))), weights = np.array([1.0/len(ene_dens)]*len(ene_dens)))
    pl.ylabel('Freq')

    pl.figure(4)
    pl.suptitle('susc')
    pl.subplot(311)
    pl.scatter(k_susc, corr_susc)
    pl.ylabel('C')
    pl.subplot(312)
    pl.scatter(k_susc, tau_susc)
    pl.ylabel('tau')
    pl.xlabel('k')
    pl.subplot(313)
    pl.hist(susc, bins = int(np.sqrt(len(susc))), weights = np.array([1.0/len(susc)]*len(susc)))
    pl.ylabel('Freq')

    pl.figure(5)
    pl.suptitle('G_pm')
    pl.subplot(311)
    pl.scatter(k_Gpm, corr_Gpm)
    pl.ylabel('C')
    pl.subplot(312)
    pl.scatter(k_Gpm, tau_Gpm)
    pl.ylabel('tau')
    pl.xlabel('k')
    pl.subplot(313)
    pl.hist(G_pm, bins = int(np.sqrt(len(G_pm))), weights = np.array([1.0/len(G_pm)]*len(G_pm)))
    pl.ylabel('Freq')

    pl.figure(6)
    pl.suptitle('mu2')
    pl.subplot(311)
    pl.scatter(k_mu2, corr_mu2)
    pl.ylabel('C')
    pl.subplot(312)
    pl.scatter(k_mu2, tau_mu2)
    pl.ylabel('tau')
    pl.xlabel('k')
    pl.subplot(313)
    pl.hist(mu2, bins = int(np.sqrt(len(mu2))), weights = np.array([1.0/len(mu2)]*len(mu2)))
    pl.ylabel('Freq')



pl.show()
