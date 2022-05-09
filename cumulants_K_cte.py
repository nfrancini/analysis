# PACCHETTI UTILI
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from tqdm import tqdm
import time

# DEFINISCO LA FUNZIONE DI ANALISI PER K3
def K3_analysis(j_min, j_max, l_min, n_term_min, n_term_max, jc_init, y_init):
    # DEFINISCO LE FUNZIONI DI FIT
    def poly(X, a, b, d, *c):
        k,l = X
        poly = 0
        for i in range(n_term):
            poly = poly + c[i] * ((k-a)*l**(b))**i
        return (poly)*l**d

    # APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
    L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len, K2_g, err_K2_g, K2_sp, err_K2_sp, K3_g, err_K3_g, K3_sp, err_K3_sp = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

    # TAGLIO LE MISURE ENTRO L'INTERVALLO DESIDERATO
    mask = ((J>=j_min) & (J<=j_max) & (L>=l_min))
    J = J[mask]
    L = L[mask]
    K3_g = K3_g[mask]
    err_K3_g = err_K3_g[mask]

    # RIORDINO IL FILE
    sorted_J = K[J.argsort()]
    sorted_L = L[J.argsort()]
    sorted_K3_g = K3_g[J.argsort()]
    sorted_err_K3_g = err_K3_g[J.argsort()]

    # CREO ARRAY DI SUPPORTO PER SALVARE I VARI RISULTATI
    Jc = np.array([])
    err_Jc = np.array([])
    yk = np.array([])
    err_yk = np.array([])
    theta3 = np.array([])
    err_theta3 = np.array([])
    red_chisq = np.array([])
    term = np.array([])
    coeff = []

    c = np.zeros(n_term_min)
    initParam = [jc_init, y_init, 3*y_init, *c]

    # CICLO SUL NUMERO DI TERMINI DEL POLINOMIO
    for n_term in tqdm(range(n_term_min, n_term_max), desc = 'K3_g_fit_no_corrections'):
        inf_bounds = -np.ones(n_term)*100
        sup_bounds = np.ones(n_term)*100
        popt, pcov = curve_fit(poly, (sorted_J, sorted_L), sorted_K3_g, p0 = initParam, sigma = sorted_err_K3_g, absolute_sigma = True, bounds = ((0.07, 1.3, 3.9, *inf_bounds), (0.08, 1.6, 4.8, *sup_bounds)))

        initParam = popt
        initParam = np.append(initParam, 0.0)

        chiq = np.sum(((sorted_K3_g - poly((sorted_J, sorted_L), *popt)) / (sorted_err_K3_g))**2)
        ndof = len(sorted_K3_g) - len(initParam)
        red_chiq = chiq/ndof
        perr = np.sqrt(np.diag(pcov))
        Jc = np.append(Jc, popt[0])
        err_Jc = np.append(err_Jc, perr[0])
        yk = np.append(yk, popt[1])
        err_yk = np.append(err_yk, perr[1])
        theta3 = np.append(theta3, popt[2])
        err_theta3 = np.append(err_theta3, perr[2])
        term = np.append(term, n_term)
        red_chisq = np.append(red_chisq, red_chiq)
        coeff.append(popt[3:])

    # PARAMETRI OTTIMALI DA RESTITUIRE
    delta = np.abs(1-red_chisq)
    c = np.array(coeff[np.argmin(delta)])
    popt = [Jc[np.argmin(delta)], yk[np.argmin(delta)], theta3[np.argmin(delta)], *c]
    err_opt = [err_Jc[np.argmin(delta)], err_yk[np.argmin(delta)], err_theta3[np.argmin(delta)]]
    n_term_K3_g = int(term[np.argmin(delta)])
    red_chisq_opt = red_chisq[np.argmin(delta)]

    return popt, err_opt, n_term_K3_g, red_chisq_opt

# DEFINISCO LA FUNZIONE DI ANALISI CON CORREZIONI
# PER ORA SOLO CORREZIONI CON SCALING IN -OMEGA
def K3_analysis_corrections(j_min, j_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps, jc_init, y_init, shift):
    def poly_corrections(X, a, b, d, *c):
        k,l = X
        poly = 0
        for i in range(n_term1):
            poly = poly + c[i] * ((k-a)*l**(b))**i
        for i in range(n_term2):
            poly = poly + (c[i+n_term1] * ((k-a)*l**(b))**i)*(l**(-omega))
        return (poly)*l**(d)

    # APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
    L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len, K2_g, err_K2_g, K2_sp, err_K2_sp, K3_g, err_K3_g, K3_sp, err_K3_sp = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

    # TAGLIO LE MISURE ENTRO L'INTERVALLO DESIDERATO
    mask = ((J>=j_min) & (J<=j_max) & (L>=l_min))
    J = J[mask]
    L = L[mask]
    K3_g = K3_g[mask]
    err_K3_g = err_K3_g[mask]

    # RIORDINO IL FILE
    sorted_J = J[J.argsort()]
    sorted_L = L[J.argsort()]
    sorted_K3_g = K3_g[J.argsort()]
    sorted_err_K3_g = err_K3_g[J.argsort()]

    Jc = np.array([])
    err_Jc = np.array([])
    yk = np.array([])
    err_yk = np.array([])
    theta3 = np.array([])
    err_theta3 = np.array([])
    red_chisq = np.array([])
    term1 = np.array([])
    term2 = np.array([])
    omg = np.array([])
    coeff = []

    c = np.zeros(n_term_min+n_term_min-shift, dtype = 'float')
    initParam = [jc_init, y_init, 3*y_init, *c]
    for omega in tqdm(np.arange(omega_min, omega_max, omega_step), desc = 'K3_g_fit_w_corrections'):
        for n_term1 in range(n_term_min, n_term_max):
            for n_term2 in range(n_term_min-shift, n_term_max-shift):
                inf_bounds = -np.ones(n_term1 + n_term2)*100
                sup_bounds = np.ones(n_term1 + n_term2)*100
                popt, pcov = curve_fit(poly_corrections, (sorted_J, sorted_L), sorted_K3_g, p0 = initParam, sigma = sorted_err_K3_g, absolute_sigma = True, bounds = ((0.254, 1.54, 4.6, *inf_bounds), (0.257, 1.63, 5.0, *sup_bounds)))

                initParam = popt
                initParam = np.append(initParam, 0.0)

                chiq = np.sum(((sorted_K3_g - poly_corrections((sorted_J, sorted_L), *popt)) / (sorted_err_K3_g))**2)
                ndof = len(sorted_K3_g) - len(initParam)
                red_chiq = chiq/ndof

                Jc = np.append(Jc, popt[0])
                term1 = np.append(term1, n_term1)
                term2 = np.append(term2, n_term2)
                yk = np.append(yk, popt[1])
                theta3 = np.append(theta3, popt[2])

                perr = np.sqrt(np.diag(pcov))
                err_Jc = np.append(err_Jc, perr[0])
                err_yk = np.append(err_yk, perr[1])
                err_theta3 = np.append(err_theta3, perr[2])

                red_chisq = np.append(red_chisq, red_chiq)
                omg = np.append(omg, omega)
                coeff.append(popt[3:])

            c = np.zeros(n_term1+n_term_min-shift+1, dtype = 'float')
            initParam = [jc_init, y_init, 3*y_init, *c]

        c = np.zeros(n_term_min+n_term_min-shift, dtype = 'float')
        initParam = [jc_init, y_init, 3*y_init, *c]

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
    yk_masked = yk[mask]
    err_yk_masked = err_yk[mask]
    theta3_masked = theta3[mask]
    err_theta3_masked = err_theta3[mask]

    mean_chisq_red = np.mean(red_chisq_masked)

    mean_Jc = np.mean(Jc_masked)
    std_Jc = (np.max(Jc_masked)-np.min(Jc_masked))/2.0

    mean_yk = np.mean(yk_masked)
    std_yk = (np.max(yk_masked)-np.min(yk_masked))/2.0

    mean_theta3 = np.mean(theta3_masked)
    std_theta3 = (np.max(theta3_masked) - np.min(theta3_masked))/2.0

    # PARAMETRI PER PLOT, DA NON PRENDERE TROPPO SUL SERIO
    cc = np.array(coeff[np.argmin(delta)])
    plot_params = [Kc[np.argmin(delta)], yk[np.argmin(delta)], *cc]
    n_term1_K3_g = int(term1[np.argmin(delta)])
    n_term2_K3_g = int(term2[np.argmin(delta)])
    omega_plot = omg[np.argmin(delta)]

    return plot_params, n_term1_K3_g, n_term2_K3_g, omega_plot, mean_Jc, std_Jc, mean_yk, std_yk, mean_theta3, std_theta3, mean_chisq_red
