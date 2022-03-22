import numpy as np
import pylab as pl
import sys
import os
from scipy.optimize import curve_fit
from matplotlib.pyplot import cm
from scipy.odr import odrpack
from itertools import cycle
from itertools import compress
import natsort
from tqdm import tqdm

from graph_style import *
from corr_len import *
from susc import *

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

# ANALISI PER LA LUNGHEZZA DI CORRELAZIONE

# PRIMO FIT, SERVE SOPRATTUTTO AD AVERE UNA PRIMA IDEA
# VEDIAMO SOPRATTUTTO QUANTO È GRANDE IL CHI_QUADRO RIDOTTO
# E CAPIAMO SE SERVONO LE CORREZIONI
j_min = 0
j_max = 10

l_min = 7

n_term_min = 11
n_term_max = 20

popt, err_opt, n_term_xi, red_chisq_opt_xi = xi_analysis(j_min, j_max, l_min, n_term_min, n_term_max)

# SECONDO FIT, ENTRANO IN GIOCO LE CORREZIONI
# IL RISULTATO CHE OTTENGO È UN RISULTATO MEDIATO
# SU DIVERSI FIT A CHI2 RIDOTTO ENTRO UN INTERVALLO SPECIFICATO
# L'ERRORE CHE OTTENGO È LA DEVIAZIONE STANDARD DEI PARAMETRI
# E POSSO CONSIDERARLO DI FATTO COME UN ERRORE SISTEMATICO (CREDO)
j_min = 0
j_max = 10

n_term_min = 10
n_term_max = 16

l_min = 7

omega_min = 0.5
omega_max = 1.5
omega_step = 0.1

eps = 0.2

plot_params_xi, n_term1_xi, n_term2_xi, omega_plot_xi, mean_Jc, std_Jc, mean_nu, std_nu, mean_chisq_red_xi = xi_analysis_corrections(j_min, j_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps)

# ANALISI PER LA SUSCETTIVITÀ
# ANCHE IN QUESTO CASO IL PRIMO FIT È
# PER AVERE UN'IDEA DELLA BONTÀ DELLO SCALING E
# DELLA NECESSITÀ DELLE CORREZIONI
j_min = 0
j_max = 10

l_min = 7

n_term_min = 8
n_term_max = 20

popt1, err_opt1, n_term_susc, red_chisq_opt_susc = susc_analysis(j_min, j_max, l_min, n_term_min, n_term_max)

# SECONDO FIT, ENTRANO IN GIOCO LE CORREZIONI
# IL RISULTATO CHE OTTENGO È UN RISULTATO MEDIATO
# SU DIVERSI FIT A CHI2 RIDOTTO ENTRO UN INTERVALLO SPECIFICATO
# L'ERRORE CHE OTTENGO È LA DEVIAZIONE STANDARD DEI PARAMETRI
# E POSSO CONSIDERARLO DI FATTO COME UN ERRORE SISTEMATICO (CREDO)
j_min = 0
j_max = 10

l_min = 7

n_term_min = 11
n_term_max = 20

omega_min = 0.5
omega_max = 1.5
omega_step = 0.1

eps = 0.2

plot_params_susc, n_term1_susc, n_term2_susc, omega_plot_susc, mean_eta, std_eta, mean_chisq_red_susc = susc_analysis_corrections(j_min, j_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps)

# PULISCO IL TERMINALE DALLE BARRE DI CARICAMENTO
os.system('cls' if os.name == 'nt' else 'clear')

# STAMPO I RISULTATI
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
print("-----------------------------------------------------")
# print("omega_plot = %f" %  (omega_plot))
# print("termini polinomio senza correzione %d" % (n_term1_xi))
# print("termini polinomio con correzione %d" % (n_term2_xi))
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
# print("omega = %f" % (omega_opt1))
# print("termini polinomio senza correzione %d" % (n_term1_susc))
# print("termini polinomio con correzione %d" % (n_term2_susc))
print("-----------------------------------------------------")


# PLOT
c_cycle, m_cycle = set_style()

for fname in path_filenames:
    with open(fname) as infile:
        aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2, aux_err_K2, aux_K3, aux_err_K3, aux_K4, aux_err_K4 = np.genfromtxt(infile, delimiter ="\t", unpack = True)

        c = next(c_cycle)
        m = next(m_cycle)
        # CORRELAZIONI (PER ADESSO SOLO CON CORREZIONI)
        fig_1 = pl.figure(1)
        n_term1 = n_term1_xi
        n_term2 = n_term2_xi
        # n_term = n_term_xi
        omega = omega_plot_xi
        pl.errorbar(aux_J, aux_corr_len/aux_L, aux_err_corr_len/aux_L, ls='', fillstyle = 'none', color = c, marker = m, label = 'L=' + str(int(aux_L[0])))

        x = np.linspace(np.min(aux_J), np.max(aux_J), 500)
        y = np.array([aux_L[0]]*len(x))
        pl.plot(x, poly_corrections((x,y), *plot_params_xi), color = c)
        # pl.plot(x, poly((x,y), *popt))

        pl.xlabel(r'$J$')
        pl.ylabel(r'$R_{\xi}$')
        pl.legend()
        pl.savefig('fit_xi_k=0.pdf')

        fig_2 = pl.figure(2)
        pl.errorbar((aux_J-mean_Jc)*(aux_L**(1.0/mean_nu)), aux_corr_len/aux_L, aux_err_corr_len/aux_L, ls='', fillstyle = 'none', color = c, marker = m, label = 'L=' + str(int(aux_L[0])))

        pl.xlabel(r'$(J-J_c)L^{1/\nu}$')
        pl.ylabel(r'$R_{\xi}$')
        pl.legend()
        pl.savefig('scaling_xi_k=0.pdf')

        # SUSCETTIVITÀ
        fig_3 = pl.figure(3)
        pl.errorbar(aux_corr_len/aux_L, aux_susc*aux_L**(mean_eta-2), aux_err_susc*aux_L**(mean_eta-2), aux_err_corr_len/aux_L, ls='', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))

        pl.xlabel(r'$R_{\xi}$')
        pl.ylabel(r'$\chi L^{\eta - 2}$')
        pl.legend()
        pl.savefig('scaling_susc_k=0.pdf')

        # BINDER
        fig_4 = pl.figure(4)
        pl.errorbar(aux_corr_len/aux_L, aux_U, aux_err_U, aux_err_corr_len/aux_L, ls='', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))

        pl.xlabel(r'$R_{\xi}$')
        pl.ylabel(r'$U$')
        pl.legend()
        pl.savefig('scaling_U_k=0.pdf')
pl.show()

# ELIMINO IL FILE TEMPORANEO
os.remove("./temp.dat")
