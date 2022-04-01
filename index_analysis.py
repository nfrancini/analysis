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
from cumulants import *
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

def poly_cumulants(X, a, b, d, *c):
    k,l = X
    poly = 0
    for i in range(n_term):
        poly = poly + c[i] * ((k-a)*l**(b))**i
    return (poly)*l**d

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

# IMPOSTARE QUI QUALE ANALISI SI VUOLE EFFETTUARE
an_xi = False
an_susc = False
an_cumul_K2 = False
an_cumul_K3 = False
an_binder = True
FO_intersect = True

# ANALISI PER LA LUNGHEZZA DI CORRELAZIONE

# PRIMO FIT, SERVE SOPRATTUTTO AD AVERE UNA PRIMA IDEA
# VEDIAMO SOPRATTUTTO QUANTO È GRANDE IL CHI_QUADRO RIDOTTO
# E CAPIAMO SE SERVONO LE CORREZIONI

if(an_xi == True):
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

if(an_xi == True):
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

if(an_susc == True):
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

if(an_susc == True):
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

# ANALISI DEI CUMULANTI
# PER IL MOMENTO SENZA CORREZIONI
# POI INTRODUCO LE CORREZIONI SUCCESSIVAMENTE

if(an_cumul_K2 == True):
    k_min = 0.05
    k_max = 0.10

    l_min = 7

    n_term_min = 8
    n_term_max = 15

    popt_K2, err_opt_K2, n_term_K2, red_chisq_opt_K2 = K2_analysis(k_min, k_max, l_min, n_term_min, n_term_max)

if(an_cumul_K3 == True):
    k_min = 0.05
    k_max = 0.10

    l_min = 7

    n_term_min = 8
    n_term_max = 15

    popt_K3, err_opt_K3, n_term_K3, red_chisq_opt_K3 = K3_analysis(k_min, k_max, l_min, n_term_min, n_term_max)

# ANALISI DEI CUMULANTI
# METTO ANCHE LE CORREZIONI

if(an_cumul_K2 == True):
    k_min = 0.060
    k_max = 0.09

    l_min = 9

    n_term_min = 2
    n_term_max = 8

    omega_min = 0.5
    omega_max = 1.0
    omega_step = 0.05

    eps = 100

    plot_params_K2, n_term1_K2, n_term2_K2, n_term_K3, omega_plot_K2, mean_Kc_K2, std_Kc_K2, mean_yk_K2, std_yk_K2, mean_chisq_red_K2 = K2_analysis_corrections(k_min, k_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps)

if(an_cumul_K3 == True):
    k_min = 0.060
    k_max = 0.09

    l_min = 9

    n_term_min = 2
    n_term_max = 8

    omega_min = 0.5
    omega_max = 1.0
    omega_step = 0.05

    eps = 100

    plot_params_K3, n_term1_K3, n_term2_K3, omega_plot_K3, mean_Kc_K3, std_Kc_K3, mean_yk_K3, std_yk_K3, mean_theta3, std_theta3, mean_chisq_red_K3 = K3_analysis_corrections(k_min, k_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps)


# ANALISI DELLA TRANSIZIONE  DI PRIMO ORDINE
if(FO_intersect ==  True):

    ctrl_plot = False

    j_min = 0.320
    j_max = 0.384

    n_term_min = 2
    n_term_max = 11

    inter_min = 0.34
    inter_max = 0.36

    mean_intersect, err_intersect, n_term_FO = FO_intersection(path_filenames, j_min, j_max, n_term_min, n_term_max, inter_min, inter_max, ctrl_plot)

# PULISCO IL TERMINALE DALLE BARRE DI CARICAMENTO
# os.system('cls' if os.name == 'nt' else 'clear')

# STAMPO I RISULTATI
if (an_xi == True):
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

if(an_susc == True):
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

if(an_cumul_K2 == True):
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
    print("termini polinomio parte di background %d" % (n_term_K2))
    print("-----------------------------------------------------")

if(an_cumul_K3 == True):
    print("-----------------------------------------------------")
    print("K3: RISULTATI SENZA CORREZIONE")
    print("-----------------------------------------------------")
    print("Kc = %f +- %f" % (popt_K3[0], err_opt_K3[0]))
    print("y_k = %f +- %f" % (popt_K3[1], err_opt_K3[1]))
    print("theta_3 = %f +- %f" % (popt_K3[2], err_opt_K3[2]))
    print("red_chisq = %f" % (red_chisq_opt_K3))
    print("termini polinomio %d" % (n_term_K3))
    print("-----------------------------------------------------")

    print("-----------------------------------------------------")
    print("K3: RISULTATI CON CORREZIONE")
    print("-----------------------------------------------------")
    print("Kc (medio) = %f +- %f" % (mean_Kc_K3, std_Kc_K3))
    print("y_k (medio) = %f +- %f" % (mean_yk_K3, std_yk_K3))
    print("theta_3 (medio) = %f +- %f" % (mean_theta3, std_theta3))
    print("red_chisq = %f" % (mean_chisq_red_K3))
    print("termini polinomio senza correzione %d" % (n_term1_K3))
    print("termini polinomio con correzione %d" % (n_term2_K3))
    print("-----------------------------------------------------")

if(FO_intersect == True):
    print("-----------------------------------------------------")
    print("INTERSEZIONE DELLE CURVE U(J) NELLA TRANSIZIONE DEL PRIMO ORDINE")
    print("-----------------------------------------------------")
    print("Jc = %f +- %f" % (mean_intersect, err_intersect))
    print("termini polinomio %d" % (n_term_FO))
    print("-----------------------------------------------------")

# PLOT
c_cycle, m_cycle = set_style()

for fname in path_filenames:
    with open(fname) as infile:
        aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)

        c = next(c_cycle)
        m = next(m_cycle)
        # pl.figure(10)
        # pl.errorbar(aux_K, aux_ene_dens, aux_err_ene_dens, ls='', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
        # pl.xlabel(r'$\kappa$')
        # pl.ylabel(r'$\epsilon$')
        # pl.legend()

        # LUNGHEZZA DI CORRELAZIONE (PER ADESSO SOLO CON CORREZIONI)
        if(an_xi == True):
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
        if(an_susc == True):
            fig_3 = pl.figure(3)
            pl.errorbar(aux_corr_len/aux_L, aux_susc*aux_L**(mean_eta-2), aux_err_susc*aux_L**(mean_eta-2), aux_err_corr_len/aux_L, ls='', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))

            pl.xlabel(r'$R_{\xi}$')
            pl.ylabel(r'$\chi L^{\eta - 2}$')
            pl.legend()
            pl.savefig('scaling_susc_k=0.pdf')

        # CUMULANTI
        # if(an_cumul == True):
        #     fig_8 = pl.figure(8)
        #     pl.errorbar((aux_K-mean_Kc_K2)*(aux_L)**(mean_yk_K2), aux_K2/(aux_L**(2*mean_yk_K2)), aux_err_K2/(aux_L**mean_theta2), ls = '', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
        #     pl.xlabel(r'$(\kappa -\kappa_c)L^y_k$')
        #     pl.ylabel(r'$K_2/L^{\theta_2}$')
        #     pl.legend()
        #
        #     fig_9 = pl.figure(9)
        #     pl.errorbar((aux_K-mean_Kc_K3)*(aux_L)**(mean_yk_K3), aux_K3/(aux_L**mean_theta3), aux_err_K3/(aux_L**mean_theta3), ls = '', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
        #     pl.xlabel(r'$(\kappa - \kappa_c)L^y_k$')
        #     pl.ylabel(r'$K_3/L^{\theta_3}$')
        #     pl.legend()


        # BINDER
        if(an_binder == True):

            fig_4 = pl.figure(4)
            pl.errorbar(aux_corr_len/aux_L, aux_U, aux_err_U, aux_err_corr_len/aux_L, ls='', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
            pl.xlabel(r'$R_{\xi}$')
            pl.ylabel(r'$U$')
            pl.legend()
            # pl.savefig('./grafici/k=0.4/binder_FO_scaling.pdf')

            if(FO_intersect == True):
                fig_7 = pl.figure(7)
                sorted_J = aux_J[aux_J.argsort()]
                sorted_U = aux_U[aux_J.argsort()]
                sorted_err_U = aux_err_U[aux_J.argsort()]
                pl.errorbar(sorted_J, sorted_U, sorted_err_U, ls='-', color = c, marker = m, fillstyle = 'none', label = 'L=' + str(int(aux_L[0])))
                pl.xlabel(r'$J$')
                pl.ylabel(r'$U$')
                pl.legend()
                # pl.savefig('./grafici/k=0.4/binder_intersect.pdf')

pl.show()

# ELIMINO IL FILE TEMPORANEO
os.remove("./temp.dat")
