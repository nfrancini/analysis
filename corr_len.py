# PACCHETTI UTILI
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from tqdm import tqdm
import time

# DEFINISCO LA FUNZIONE DI ANALISI
def xi_analysis(j_min, j_max, l_min, n_term_min, n_term_max):
    # DEFINISCO LE FUNZIONI DI FIT
    def poly(X, a, b, *c):
        j,l = X
        poly = 0
        for i in range(n_term):
            poly = poly + c[i] * ((j-a)*l**(1/b))**i
        return poly

    # APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
    L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

    # TAGLIO LE MISURE ENTRO L'INTERVALLO DESIDERATO
    mask = ((J>=j_min) & (J<=j_max) & (L>=l_min))
    J = J[mask]
    L = L[mask]
    corr_len = corr_len[mask]
    err_corr_len = err_corr_len[mask]

    # CREO ARRAY DI SUPPORTO PER SALVARE I VARI RISULTATI
    Jc = np.array([])
    err_Jc = np.array([])
    nu = np.array([])
    err_nu = np.array([])
    red_chisq = np.array([])
    term = np.array([])
    coeff = []

    c = np.ones(n_term_min)/2.0
    initParam = [0.70, 0.7117, *c]

    # CICLO SUL NUMERO DI TERMINI DEL POLINOMIO
    for n_term in tqdm(range(n_term_min, n_term_max), desc = 'xi_fit_no_corrections'):
        inf_bounds = -np.ones(n_term)
        sup_bounds = np.ones(n_term)
        popt, pcov = curve_fit(poly, (J,L), corr_len/L, p0 = initParam, sigma =err_corr_len/L, absolute_sigma = True, bounds = ((0.650, 0.650, *inf_bounds), (0.750, 0.750, *sup_bounds)))

        initParam = popt
        initParam = np.append(initParam, 0.5)

        chiq = np.sum(((corr_len/L - poly((J,L), *popt)) / (err_corr_len/L))**2)
        ndof = len(corr_len) - len(initParam)
        red_chiq = chiq/ndof
        perr = np.sqrt(np.diag(pcov))
        Jc = np.append(Jc, popt[0])
        err_Jc = np.append(err_Jc, perr[0])
        nu = np.append(nu, popt[1])
        err_nu = np.append(err_nu, perr[1])
        term = np.append(term, n_term)
        red_chisq = np.append(red_chisq, red_chiq)
        coeff.append(popt[2:])

    # PARAMETRI OTTIMALI DA RESTITUIRE
    delta = np.abs(1-red_chisq)
    c = np.array(coeff[np.argmin(delta)])
    popt = [Jc[np.argmin(delta)], nu[np.argmin(delta)], *c]
    err_opt = [err_Jc[np.argmin(delta)], err_nu[np.argmin(delta)]]
    n_term_csi = int(term[np.argmin(delta)])
    red_chisq_opt = red_chisq[np.argmin(delta)]

    return popt, err_opt, n_term_csi, red_chisq_opt

# DEFINISCO LA FUNZIONE DI ANALISI CON CORREZIONI
def xi_analysis_corrections(j_min, j_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps):
    def poly_corrections(X, a, b, *c):
        j,l = X
        poly = 0
        for i in range(n_term1):
            poly = poly + c[i] * ((j-a)*l**(1/b))**i
        for i in range(n_term2):
            poly = poly + (c[i+n_term1] * ((j-a)*l**(1/b))**i)*(l**(-omega))
        return poly

    # APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
    L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

    # TAGLIO LE MISURE ENTRO L'INTERVALLO DESIDERATO
    mask = ((J>=j_min) & (J<=j_max) & (L>=l_min))
    J = J[mask]
    L = L[mask]
    corr_len = corr_len[mask]
    err_corr_len = err_corr_len[mask]

    Jc = np.array([])
    err_Jc = np.array([])
    nu = np.array([])
    err_nu = np.array([])
    red_chisq = np.array([])
    term1 = np.array([])
    term2 = np.array([])
    omg = np.array([])
    coeff = []

    c = np.zeros(n_term_min+n_term_min, dtype = 'float')
    initParam = [0.70, 0.7117, *c]
    for omega in tqdm(np.arange(omega_min, omega_max, omega_step), desc = 'xi_fit_w_corrections'):
        for n_term1 in range(n_term_min, n_term_max):
            for n_term2 in range(n_term_min, n_term_max):
                inf_bounds = -np.ones(n_term1 + n_term2)
                sup_bounds = np.ones(n_term1 + n_term2)
                popt, pcov = curve_fit(poly_corrections, (J,L), corr_len/L, p0 = initParam, sigma = err_corr_len/L, absolute_sigma = True, bounds = ((0.690, 0.700, *inf_bounds), (0.710, 0.720, *sup_bounds)))

                initParam = popt
                initParam = np.append(initParam, 0.0)

                chiq = np.sum(((corr_len/L - poly_corrections((J,L), *popt)) / (err_corr_len/L))**2)
                ndof = len(corr_len) - len(initParam)
                red_chiq = chiq/ndof

                Jc = np.append(Jc, popt[0])
                term1 = np.append(term1, n_term1)
                term2 = np.append(term2, n_term2)
                nu = np.append(nu, popt[1])

                perr = np.sqrt(np.diag(pcov))
                err_Jc = np.append(err_Jc, perr[0])
                err_nu = np.append(err_nu, perr[1])

                red_chisq = np.append(red_chisq, red_chiq)
                omg = np.append(omg, omega)
                coeff.append(popt[2:])

            c = np.zeros(n_term1+n_term_min+1, dtype = 'float')
            initParam = [0.70, 0.7117, *c]

        c = np.zeros(n_term_min+n_term_min, dtype = 'float')
        initParam = [0.70, 0.7117, *c]

    # PARAMETRI OTTIMALI DA RESTITUIRE
    # NEI PASSI SEGUENTI TENGO CONTO
    # DELLE SISTEMATICHE, CIOÈ TENENDO IL CHI^2
    # ENTRO UN CERTO INTERVALLO GUARDO COME VARIANO
    # I PARAMETRI PER DIVERSI GRADI, OMEGA ..
    delta = np.abs(1-red_chisq)
    mask = ((delta<=eps))
    red_chisq_masked = red_chisq[mask]
    Jc_masked = Jc[mask]
    err_Jc_masked = err_Jc[mask]
    nu_masked = nu[mask]
    err_nu_masked = err_nu[mask]

    mean_chisq_red = np.mean(red_chisq_masked)

    mean_Jc = np.mean(Jc_masked)
    std_Jc = (np.max(Jc_masked)-np.min(Jc_masked))/2.0

    mean_nu = np.mean(nu_masked)
    std_nu = (np.max(nu_masked)-np.min(nu_masked))/2.0

    # PARAMETRI PER PLOT, DA NON PRENDERE TROPPO SUL SERIO
    cc = np.array(coeff[np.argmin(delta)])
    plot_params = [Jc[np.argmin(delta)], nu[np.argmin(delta)], *cc]
    n_term1_xi = int(term1[np.argmin(delta)])
    n_term2_xi = int(term2[np.argmin(delta)])
    omega_plot = omg[np.argmin(delta)]

    return plot_params, n_term1_xi, n_term2_xi, omega_plot, mean_Jc, std_Jc, mean_nu, std_nu, mean_chisq_red
