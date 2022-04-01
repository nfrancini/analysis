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

from graph_style import *
from corr_len import *
from susc import *
from FO_intersection import *

# DEFINISCO LE FUNZIONI PER I PLOT
def poly(X, a, b, d, *c):
    k,l = X
    poly = 0
    for i in range(n_term):
        poly = poly + c[i] * ((k-a)*l**(b))**i
    return (poly)*l**d

def poly_analytic_corrections(X, a, b, *c):
    k,l = X
    poly = 0
    poly2 = 0
    for i in range(n_term1):
        poly = poly + c[i+n_term3] * ((k-a)*l**(b))**i
    for i in range(n_term2):
        poly = poly + (c[i+n_term1+n_term3] * ((k-a)*l**(b))**i)*(l**(-omega))
    for i in range(n_term3):
        poly2 = poly2 + c[i]*(k)**i
    return (poly)*l**(2*b) + poly2*l**3

# DEFINISCO I PARAMETRI PER L'ANALISI DATI E I VALORI INIZIALI CHE MI ASPETTO
k_min = 0.05
k_max = 0.09

l_min = 7

n_term_min = 8
n_term_max = 15

kc_init = 0.076
y_init = 1.488

ctrl_data = False
ctrl_K2_g = False
ctrl_K3_g = False
ctrl_FO = False

# CREO FILE UNICO TEMPORANEO PER I DATI
# MEMORIZZO TUTTI I NOMI DEI FILE IN CARTELLA
# SUCCESSIVAMENTE CREO UN FILE TEMPORANEO DOVE METTO INSIEME TUTTI I DATI
basepath = sys.argv[1]
filenames = natsort.natsorted(os.listdir(basepath))
path_filenames=[]
for entry in filenames:
    path_filenames.append(os.path.join(basepath, entry))

with open("./temp.dat", "w") as tempfile:
    for fname in path_filenames:
        with open(fname) as infile:
            for line in infile:
                tempfile.write(line)

# APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
# L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len, K2_g, err_K2_g, K2_sp, err_K2_sp, K3_g, err_K3_g, K3_sp, err_K3_sp = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

# PLOT DI CONTROLLO INIZIALE
if (ctrl_data == True):
    for fname in path_filenames:
        with open(fname) as infile:
            # aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)
            aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len = np.genfromtxt(infile, delimiter ="\t", unpack = True)

            # CUMULANTE K2
            fig_1 = pl.figure(1)
            pl.errorbar(aux_K, aux_K2_g, aux_err_K2_g, ls='', marker = 'o', fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$\kappa$')
            pl.ylabel(r'$K_{2g}$')
            pl.legend()

            # CUMULANTE K2 RISCALATO
            fig_2 = pl.figure(2)
            pl.errorbar((aux_K - kc_init)*(aux_L)**(y_init), aux_K2_g/(aux_L**(2*y_init)), aux_err_K2_g/(aux_L**(2*y_init)) ls='', marker = 'o', fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$(\kappa - \kappa_c)L^{y}$')
            pl.ylabel(r'$K_{2g}/L^{2y}$')
            pl.legend()

            # CUMULANTE K3
            fig_3 = pl.figure(3)
            pl.errorbar(aux_K, aux_K3_g, aux_err_K3_g, ls='', marker = 'o', fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$\kappa$')
            pl.ylabel(r'$K_{3g}$')
            pl.legend()

            # CUMULANTE K3 RISCALATO
            fig_4 = pl.figure(4)
            pl.errorbar((aux_K - kc_init)*(aux_L)**(y_init), aux_K3_g/(aux_L**(3*y_init)), aux_err_K3_g/(aux_L**(3*y_init)) ls='', marker = 'o', fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$(\kappa - \kappa_c)L^{y}$')
            pl.ylabel(r'$K_{3g}/L^{3y}$')
            pl.legend()

# FIT SENZA E CON CORREZIONI PER TRANSIZIONE CONTINUA
# IN QUESTO MODO POSSO CALCOLARE ALCUNI INDICI CRITICI ED IL PUNTO CRITICO
if (ctrl_K2_g == True):
    popt_K2, err_opt_K2, n_term_K2, red_chisq_opt_K2 = K2_analysis(k_min, k_max, l_min, n_term_min, n_term_max, kc_init, y_init)

    omega_min = 0.5
    omega_max = 1.0
    omega_step = 0.05

    eps = 100

    plot_params_K2, n_term1_K2, n_term2_K2, n_term3_K2, omega_plot_K2, mean_Kc_K2, std_Kc_K2, mean_yk_K2, std_yk_K2, mean_chisq_red_K2 = K2_analysis_corrections(k_min, k_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps, kc_init, y_init)

    # STAMPO I RISULTATI
    os.system('cls' if os.name == 'nt' else 'clear')
    print("-----------------------------------------------------")
    print("K2: RISULTATI SENZA CORREZIONE")
    print("-----------------------------------------------------")
    print("Kc = %f +- %f" % (popt_K2[0], err_opt_K2[0]))
    print("y_k = %f +- %f" % (popt_K2[1], err_opt_K2[1]))
    print("theta_2 = %f +- %f" % (popt_K2[2], err_opt_K2[2]))
    print("red_chisq = %f" % (red_chisq_opt_K2))
    print("termini polinomio %d" % (n_term_K2))
    print("-----------------------------------------------------")

    print("-----------------------------------------------------")
    print("K2: RISULTATI CON CORREZIONE")
    print("-----------------------------------------------------")
    print("Kc (medio) = %f +- %f" % (mean_Kc_K2, std_Kc_K2))
    print("y_k (medio) = %f +- %f" % (mean_yk_K2, std_yk_K2))
    # print("theta_2 (medio) = %f +- %f" % (mean_theta2, std_theta2))
    print("red_chisq = %f" % (mean_chisq_red_K2))
    print("omega = %f" % (omega_plot_K2))
    print("termini polinomio senza correzione %d" % (n_term1_K2))
    print("termini polinomio con correzione %d" % (n_term2_K2))
    print("termini polinomio parte di background %d" % (n_term3_K2))
    print("-----------------------------------------------------")

    # PLOT DI VERIFICA
    c_cycle, m_cycle = set_style()

    for fname in path_filenames:
        with open(fname) as infile:
            # aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)
            aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len = np.genfromtxt(infile, delimiter ="\t", unpack = True)

            c = next(c_cycle)
            m = next(m_cycle)

            # CUMULANTE K2
            fig_1 = pl.figure(1)
            pl.errorbar(aux_K, aux_K2_g, aux_err_K2_g, ls='', fillstyle = 'none', color = c, marker = m, label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$\kappa$')
            pl.ylabel(r'$K_{2g}$')
            pl.legend()

            # CUMULANTE K2 RISCALATO
            fig_2 = pl.figure(2)
            pl.errorbar((aux_K - mean_Kc_K2)*(aux_L)**(mean_yk_K2), aux_K2_g/(aux_L**(2*mean_yk_K2)), aux_err_K2_g/(aux_L**(2*mean_yk_K2)) ls='', fillstyle = 'none', color = c, marker = m, label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$(\kappa - \kappa_c)L^{y}$')
            pl.ylabel(r'$K_{2g}/L^{2y}$')
            pl.legend()

if (ctrl_K3_g == True):
    popt_K3, err_opt_K3, n_term_K3, red_chisq_opt_K3 = K3_analysis(k_min, k_max, l_min, n_term_min, n_term_max, kc_init, y_init)

    omega_min = 0.5
    omega_max = 1.0
    omega_step = 0.05

    eps = 100

    plot_params_K3, n_term1_K3, n_term2_K3, n_term3_K3, omega_plot_K3, mean_Kc_K3, std_Kc_K3, mean_yk_K3, std_yk_K3, mean_chisq_red_K3 = K3_analysis_corrections(k_min, k_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps, kc_init, y_init)

    # STAMPO I RISULTATI
    os.system('cls' if os.name == 'nt' else 'clear')
    print("-----------------------------------------------------")
    print("K3: RISULTATI SENZA CORREZIONE")
    print("-----------------------------------------------------")
    print("Kc = %f +- %f" % (popt_K3[0], err_opt_K3[0]))
    print("y_k = %f +- %f" % (popt_K3[1], err_opt_K3[1]))
    print("theta_2 = %f +- %f" % (popt_K3[2], err_opt_K3[2]))
    print("red_chisq = %f" % (red_chisq_opt_K3))
    print("termini polinomio %d" % (n_term_K3))
    print("-----------------------------------------------------")

    print("-----------------------------------------------------")
    print("K3: RISULTATI CON CORREZIONE")
    print("-----------------------------------------------------")
    print("Kc (medio) = %f +- %f" % (mean_Kc_K3, std_Kc_K3))
    print("y_k (medio) = %f +- %f" % (mean_yk_K3, std_yk_K3))
    # print("theta_2 (medio) = %f +- %f" % (mean_theta2, std_theta2))
    print("red_chisq = %f" % (mean_chisq_red_K3))
    print("omega = %f" % (omega_plot_K3))
    print("termini polinomio senza correzione %d" % (n_term1_K3))
    print("termini polinomio con correzione %d" % (n_term2_K3))
    print("termini polinomio parte di background %d" % (n_term3_K3))
    print("-----------------------------------------------------")

    # PLOT DI VERIFICA
    c_cycle, m_cycle = set_style()

    for fname in path_filenames:
        with open(fname) as infile:
            # aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)
            aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len = np.genfromtxt(infile, delimiter ="\t", unpack = True)

            c = next(c_cycle)
            m = next(m_cycle)

            # CUMULANTE K3
            fig_1 = pl.figure(1)
            pl.errorbar(aux_K, aux_K3_g, aux_err_K3_g, ls='', fillstyle = 'none', color = c, marker = m, label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$\kappa$')
            pl.ylabel(r'$K_{3g}$')
            pl.legend()

            # CUMULANTE K3 RISCALATO
            fig_2 = pl.figure(2)
            pl.errorbar((aux_K - mean_Kc_K3)*(aux_L)**(mean_yk_K3), aux_K3_g/(aux_L**(3*mean_yk_K3)), aux_err_K3_g/(aux_L**(3*mean_yk_K3)) ls='', fillstyle = 'none', color = c, marker = m, label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$(\kappa - \kappa_c)L^{y}$')
            pl.ylabel(r'$K_{3g}/L^{3y}$')
            pl.legend()


pl.show()

# ELIMINO IL FILE TEMPORANEO
os.remove("./temp.dat")
