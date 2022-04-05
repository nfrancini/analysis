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
def poly(X, a, b, *c):
    j,l = X
    poly = 0
    for i in range(n_term):
        poly = poly + c[i] * ((j-a)*l**(1/b))**i
    return poly

def poly_corrections(X, a, b, *c):
    j,l = X
    poly = 0
    for i in range(n_term1):
        poly = poly + c[i] * ((j-a)*l**(1/b))**i
    for i in range(n_term2):
        poly = poly + (c[i + n_term1] * ((j-a)*l**(1/b))**i)*(l**(-omega))
    return poly

# DEFINISCO I PARAMETRI PER L'ANALISI DATI E I VALORI INIZIALI CHE MI ASPETTO
j_min = 0.640
j_max = 0.750

l_min = 7
l_max = 35

n_term_min = 10
n_term_max = 15
shift = 0

jc_init = 0.71
nu_init = 0.7117
eta_init = 0.0378

ctrl_data = False
ctrl_continuos = True
ctrl_FO = False
ctrl_savefig = True

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
L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len, K2_g, err_K2_g, K2_sp, err_K2_sp, K3_g, err_K3_g, K3_sp, err_K3_sp = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

# L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

# PLOT DI CONTROLLO INIZIALE
if (ctrl_data == True):
    for fname in path_filenames:
        with open(fname) as infile:
            aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)
            # aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len = np.genfromtxt(infile, delimiter ="\t", unpack = True)

            # LUNGHEZZA DI CORRELAZIONE (PER ADESSO SOLO CON CORREZIONI)
            fig_1 = pl.figure(1)
            pl.errorbar(aux_J, aux_corr_len/aux_L, aux_err_corr_len/aux_L, ls='', marker = 'o', fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$J$')
            pl.ylabel(r'$R_{\xi}$')
            pl.legend()

            # SUSCETTIVITÀ
            fig_2 = pl.figure(2)
            pl.errorbar(aux_corr_len/aux_L, aux_susc*aux_L**(eta_init-2), aux_err_susc*aux_L**(eta_init-2), aux_err_corr_len/aux_L, ls='', marker = 'o', fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$R_{\xi}$')
            pl.ylabel(r'$\chi L^{\eta - 2}$')
            pl.legend()

            # BINDER
            fig_4 = pl.figure(4)
            pl.errorbar(aux_corr_len/aux_L, aux_U, aux_err_U, aux_err_corr_len/aux_L, ls='', marker = 'o', fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$R_{\xi}$')
            pl.ylabel(r'$U$')
            pl.legend()

            fig_5 = pl.figure(5)
            sorted_J = aux_J[aux_J.argsort()]
            sorted_U = aux_U[aux_J.argsort()]
            sorted_err_U = aux_err_U[aux_J.argsort()]
            pl.errorbar(sorted_J, sorted_U, sorted_err_U, ls='-', marker = 'o', fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$J$')
            pl.ylabel(r'$U$')
            pl.legend()

# FIT SENZA E CON CORREZIONI PER TRANSIZIONE CONTINUA
# IN QUESTO MODO POSSO CALCOLARE ALCUNI INDICI CRITICI ED IL PUNTO CRITICO
if (ctrl_continuos == True):
    popt, err_opt, n_term_xi, red_chisq_opt_xi = xi_analysis(j_min, j_max, l_min, n_term_min, n_term_max, jc_init, nu_init)
    popt1, err_opt1, n_term_susc, red_chisq_opt_susc = susc_analysis(j_min, j_max, l_min, n_term_min, n_term_max, eta_init)

    omega_min = 0.9
    omega_max = 1.2
    omega_step = 0.1

    eps = 0.8

    plot_params_xi, n_term1_xi, n_term2_xi, omega_plot_xi, omega_opt_xi, std_omega_xi, mean_Jc, std_Jc, mean_nu, std_nu, mean_chisq_red_xi = xi_analysis_corrections(j_min, j_max, l_min, l_max, n_term_min, n_term_max, shift, omega_min, omega_max, omega_step, eps, popt[0], popt[1])
    plot_params_susc, n_term1_susc, n_term2_susc, omega_plot_susc, omega_opt_susc, std_omega_susc, mean_eta, std_eta, mean_chisq_red_susc = susc_analysis_corrections(j_min, j_max, l_min, n_term_min, n_term_max, shift, omega_min, omega_max, omega_step, eps, popt1[0])

    # STAMPO I RISULTATI
    os.system('cls' if os.name == 'nt' else 'clear')
    print("-----------------------------------------------------")
    print("LUNGHEZZA DI CORRELAZIONE: RISULTATI SENZA CORREZIONE")
    print("-----------------------------------------------------")
    print("Jc = %f +- %f" % (popt[0], err_opt[0]))
    print("nu = %f +- %f" % (popt[1], err_opt[1]))
    print("red_chisq = %f" % (red_chisq_opt_xi))
    print("termini polinomio %d" % (n_term_xi))
    print("-----------------------------------------------------")
    print("-----------------------------------------------------")
    print("LUNGHEZZA DI CORRELAZIONE: RISULTATI CON CORREZIONE")
    print("-----------------------------------------------------")
    print("Jc (medio) = %f +- %f" % (mean_Jc, std_Jc))
    print("nu (medio) = %f +- %f" % (mean_nu, std_nu))
    print("red_chisq = %f" % (mean_chisq_red_xi))
    print("omega_opt = %f +- %f" % (omega_opt_xi, std_omega_xi))
    print("-----------------------------------------------------")
    print("-----------------------------------------------------")
    print("SUSCETTIVITÀ: RISULTATI SENZA CORREZIONE")
    print("-----------------------------------------------------")
    print("eta = %f +- %f" % (popt1[0], err_opt1))
    print("red_chisq = %f" % (red_chisq_opt_susc))
    print("termini polinomio %d" % (n_term_susc))
    print("-----------------------------------------------------")
    print("-----------------------------------------------------")
    print("SUSCETTIVITÀ: RISULTATI CON CORREZIONE")
    print("-----------------------------------------------------")
    print("eta = %f +- %f" % (mean_eta, std_eta))
    print("red_chisq = %f" % (mean_chisq_red_susc))
    print("omega_opt = %f +- %f" % (omega_opt_susc, std_omega_susc))
    # print("termini polinomio senza correzione %d" % (n_term1_susc))
    # print("termini polinomio con correzione %d" % (n_term2_susc))
    print("-----------------------------------------------------")

    # PLOT PER IL CASO CONTINUO
    c_cycle, m_cycle = set_style()

    for fname in path_filenames:
        with open(fname) as infile:
            aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)
            # aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len = np.genfromtxt(infile, delimiter ="\t", unpack = True)

            c = next(c_cycle)
            m = next(m_cycle)

            # LUNGHEZZA DI CORRELAZIONE (PER ADESSO SOLO CON CORREZIONI)
            fig_1 = pl.figure(1)
            n_term1 = n_term1_xi
            n_term2 = n_term2_xi
            omega = omega_plot_xi
            pl.errorbar(aux_J, aux_corr_len/aux_L, aux_err_corr_len/aux_L, ls='', fillstyle = 'none', color = c, marker = m, label = 'L=' + str(int(aux_L[0])))

            # x = np.linspace(j_min, j_max, 500)
            # y = np.array([aux_L[0]]*len(x))
            # pl.plot(x, poly_corrections((x,y), *plot_params_xi), color = c)
            # pl.plot(x, poly((x,y), *popt))

            pl.xlabel(r'$J$')
            pl.ylabel(r'$R_{\xi}$')
            pl.legend()
            if (ctrl_savefig == True):
                pl.savefig('./grafici/continuo/K_0.04/corr_len.png')

            fig_2 = pl.figure(2)
            pl.errorbar((aux_J-mean_Jc)*(aux_L**(1.0/mean_nu)), aux_corr_len/aux_L, aux_err_corr_len/aux_L, ls='', fillstyle = 'none', color = c, marker = m, label = 'L=' + str(int(aux_L[0])))

            pl.xlabel(r'$(J-J_c)L^{1/\nu}$')
            pl.ylabel(r'$R_{\xi}$')
            pl.legend()
            if (ctrl_savefig == True):
                pl.savefig('./grafici/continuo/K_0.04/corr_len_rescaled.png')

            # SUSCETTIVITÀ
            fig_3 = pl.figure(3)
            pl.errorbar(aux_corr_len/aux_L, aux_susc*aux_L**(mean_eta-2), aux_err_susc*aux_L**(mean_eta-2), aux_err_corr_len/aux_L, ls='', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))

            pl.xlabel(r'$R_{\xi}$')
            pl.ylabel(r'$\chi L^{\eta - 2}$')
            pl.legend()
            if (ctrl_savefig == True):
                pl.savefig('./grafici/continuo/K_0.04/susc_rescaled.png')

            # BINDER
            fig_4 = pl.figure(4)
            pl.errorbar(aux_corr_len/aux_L, aux_U, aux_err_U, aux_err_corr_len/aux_L, ls='', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$R_{\xi}$')
            pl.ylabel(r'$U$')
            pl.legend()
            if (ctrl_savefig == True):
                pl.savefig('./grafici/continuo/K_0.04/binder_rescaled.png')

if (ctrl_FO == True):
    inter_min = 0.65
    inter_max = 0.75

    ctrl_plot = True

    mean_intersect, err_intersect, n_term_FO = FO_intersection(path_filenames, j_min, j_max, n_term_min, n_term_max, inter_min, inter_max, ctrl_plot)

    print("-----------------------------------------------------")
    print("INTERSEZIONE DELLE CURVE U(J) NELLA TRANSIZIONE DEL PRIMO ORDINE")
    print("-----------------------------------------------------")
    print("Jc = %f +- %f" % (mean_intersect, err_intersect))
    print("termini polinomio %d" % (n_term_FO))
    print("-----------------------------------------------------")

    # PLOT PER IL LA TRANSIZIONE DEL PRIMO ORDINE
    # IN QUESTO CASO OSSERVO LO SCALING DI U
    # CON LA LUNGHEZZA DI CORRELAZIONE ED
    # IL PLOT DI U RISPETTO A J PER VEDERE IL PUNTO
    # DI INCONTRO
    c_cycle, m_cycle = set_style()

    for fname in path_filenames:
        with open(fname) as infile:
            # aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)
            aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len = np.genfromtxt(infile, delimiter ="\t", unpack = True)

            c = next(c_cycle)
            m = next(m_cycle)

            # BINDER
            fig_1 = pl.figure(1)
            pl.errorbar(aux_corr_len/aux_L, aux_U, aux_err_U, aux_err_corr_len/aux_L, ls='', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$R_{\xi}$')
            pl.ylabel(r'$U$')
            pl.legend()

            fig_2 = pl.figure(2)
            sorted_J = aux_J[aux_J.argsort()]
            sorted_U = aux_U[aux_J.argsort()]
            sorted_err_U = aux_err_U[aux_J.argsort()]
            pl.errorbar(sorted_J, sorted_U, sorted_err_U, ls='-', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$J$')
            pl.ylabel(r'$U$')
            pl.legend()

pl.show()

# ELIMINO IL FILE TEMPORANEO
os.remove("./temp.dat")
