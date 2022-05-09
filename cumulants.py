# PACCHETTI UTILI
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from tqdm import tqdm
import time

# DEFINISCO LA FUNZIONE DI ANALISI PER K3
def K3_analysis(k_min, k_max, l_min, n_term_min, n_term_max, kc_init, y_init):
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
    mask = ((K>=k_min) & (K<=k_max) & (L>=l_min))
    K = K[mask]
    L = L[mask]
    K3_g = K3_g[mask]
    err_K3_g = err_K3_g[mask]

    # RIORDINO IL FILE
    sorted_K = K[K.argsort()]
    sorted_L = L[K.argsort()]
    sorted_K3_g = K3_g[K.argsort()]
    sorted_err_K3_g = err_K3_g[K.argsort()]

    # CREO ARRAY DI SUPPORTO PER SALVARE I VARI RISULTATI
    Kc = np.array([])
    err_Kc = np.array([])
    yk = np.array([])
    err_yk = np.array([])
    theta3 = np.array([])
    err_theta3 = np.array([])
    red_chisq = np.array([])
    term = np.array([])
    coeff = []

    c = np.zeros(n_term_min)
    initParam = [kc_init, y_init, 3*y_init, *c]

    # CICLO SUL NUMERO DI TERMINI DEL POLINOMIO
    for n_term in tqdm(range(n_term_min, n_term_max), desc = 'K3_g_fit_no_corrections'):
        inf_bounds = -np.ones(n_term)
        sup_bounds = np.ones(n_term)
        popt, pcov = curve_fit(poly, (sorted_K, sorted_L), sorted_K3_g, p0 = initParam, sigma = sorted_err_K3_g, absolute_sigma = True, bounds = ((0.295, 1.45, 3.9, *inf_bounds), (0.305, 1.55, 4.8, *sup_bounds)))

        initParam = popt
        initParam = np.append(initParam, 0.0)

        chiq = np.sum(((sorted_K3_g - poly((sorted_K, sorted_L), *popt)) / (sorted_err_K3_g))**2)
        ndof = len(sorted_K3_g) - len(initParam)
        red_chiq = chiq/ndof
        perr = np.sqrt(np.diag(pcov))
        Kc = np.append(Kc, popt[0])
        err_Kc = np.append(err_Kc, perr[0])
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
    popt = [Kc[np.argmin(delta)], yk[np.argmin(delta)], theta3[np.argmin(delta)], *c]
    err_opt = [err_Kc[np.argmin(delta)], err_yk[np.argmin(delta)], err_theta3[np.argmin(delta)]]
    n_term_K3_g = int(term[np.argmin(delta)])
    red_chisq_opt = red_chisq[np.argmin(delta)]

    return popt, err_opt, n_term_K3_g, red_chisq_opt

# DEFINISCO LA FUNZIONE DI ANALISI CON CORREZIONI
# PER ORA SOLO CORREZIONI CON SCALING IN -OMEGA
def K3_analysis_corrections(k_min, k_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps, kc_init, y_init, shift):
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
    mask = ((K>=k_min) & (K<=k_max) & (L>=l_min))
    K = K[mask]
    L = L[mask]
    K3_g = K3_g[mask]
    err_K3_g = err_K3_g[mask]

    # RIORDINO IL FILE
    sorted_K = K[K.argsort()]
    sorted_L = L[K.argsort()]
    sorted_K3_g = K3_g[K.argsort()]
    sorted_err_K3_g = err_K3_g[K.argsort()]

    Kc = np.array([])
    err_Kc = np.array([])
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
    initParam = [kc_init, y_init, 3*y_init, *c]
    for omega in tqdm(np.arange(omega_min, omega_max, omega_step), desc = 'K3_g_fit_w_corrections'):
        for n_term1 in range(n_term_min, n_term_max):
            for n_term2 in range(n_term_min-shift, n_term_max-shift):
                print(omega, n_term1, n_term2)
                inf_bounds = -np.ones(n_term1 + n_term2)*100
                sup_bounds = np.ones(n_term1 + n_term2)*100
                popt, pcov = curve_fit(poly_corrections, (sorted_K, sorted_L), sorted_K3_g, p0 = initParam, sigma = sorted_err_K3_g, absolute_sigma = True, bounds = ((0.295, 1.45, 4.0, *inf_bounds), (0.304, 1.52, 4.65, *sup_bounds)))

                initParam = popt
                initParam = np.append(initParam, 0.0)

                chiq = np.sum(((sorted_K3_g - poly_corrections((sorted_K, sorted_L), *popt)) / (sorted_err_K3_g))**2)
                ndof = len(sorted_K3_g) - len(initParam)
                red_chiq = chiq/ndof

                Kc = np.append(Kc, popt[0])
                term1 = np.append(term1, n_term1)
                term2 = np.append(term2, n_term2)
                yk = np.append(yk, popt[1])
                theta3 = np.append(theta3, popt[2])

                perr = np.sqrt(np.diag(pcov))
                err_Kc = np.append(err_Kc, perr[0])
                err_yk = np.append(err_yk, perr[1])
                err_theta3 = np.append(err_theta3, perr[2])

                red_chisq = np.append(red_chisq, red_chiq)
                omg = np.append(omg, omega)
                coeff.append(popt[3:])

            c = np.zeros(n_term1+n_term_min-shift+1, dtype = 'float')
            initParam = [kc_init, y_init, 3*y_init, *c]

        c = np.zeros(n_term_min+n_term_min-shift, dtype = 'float')
        initParam = [kc_init, y_init, 3*y_init, *c]

    # PARAMETRI OTTIMALI DA RESTITUIRE
    # NEI PASSI SEGUENTI TENGO CONTO
    # DELLE SISTEMATICHE, CIOÈ TENENDO IL CHI^2
    # ENTRO UN CERTO INTERVALLO GUARDO COME VARIANO
    # I PARAMETRI PER DIVERSI GRADI, OMEGA ..
    delta = np.abs(1-red_chisq)
    mask = ((delta<=eps))
    red_chisq_masked = red_chisq[mask]
    Kc_masked = Kc[mask]
    err_Kc_masked = err_Kc[mask]
    yk_masked = yk[mask]
    err_yk_masked = err_yk[mask]
    theta3_masked = theta3[mask]
    err_theta3_masked = err_theta3[mask]

    mean_chisq_red = np.mean(red_chisq_masked)

    mean_Kc = np.mean(Kc_masked)
    std_Kc = (np.max(Kc_masked)-np.min(Kc_masked))/2.0

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

    return plot_params, n_term1_K3_g, n_term2_K3_g, omega_plot, mean_Kc, std_Kc, mean_yk, std_yk, mean_theta3, std_theta3, mean_chisq_red

# DEFINISCO LA FUNZIONE DI ANALISI PER K4
def K2_analysis(k_min, k_max, l_min, n_term_min, n_term_max, kc_init, y_init):
    # DEFINISCO LE FUNZIONI DI FIT
    def poly(X, a, b, d, *c):
        k,l = X
        poly = 0
        for i in range(n_term):
            poly = poly + c[i] * ((k-a)*l**(b))**i
        return (poly)*l**d

    # def poly(X, a, b, *c):
    #     k,l = X
    #     poly = 0
    #     for i in range(n_term):
    #         poly = poly + c[i] * ((k-a)*l**(b))**i
    #     return (poly)*l**(2*b)

    # APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
    L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len, K2_g, err_K2_g, K2_sp, err_K2_sp, K3_g, err_K3_g, K3_sp, err_K3_sp = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

    # TAGLIO LE MISURE ENTRO L'INTERVALLO DESIDERATO
    mask = ((K>=k_min) & (K<=k_max) & (L>=l_min))
    K = K[mask]
    L = L[mask]
    K2_g = K2_g[mask]
    err_K2_g = err_K2_g[mask]

    # RIORDINO IL FILE
    sorted_K = K[K.argsort()]
    sorted_L = L[K.argsort()]
    sorted_K2_g = K2_g[K.argsort()]
    sorted_err_K2_g = err_K2_g[K.argsort()]

    # CREO ARRAY DI SUPPORTO PER SALVARE I VARI RISULTATI
    Kc = np.array([])
    err_Kc = np.array([])
    yk = np.array([])
    err_yk = np.array([])
    theta2 = np.array([])
    err_theta2 = np.array([])
    red_chisq = np.array([])
    term = np.array([])
    coeff = []

    c = np.zeros(n_term_min)/2.0
    initParam = [kc_init, y_init, y_init*2, *c]
    # initParam = [0.07, 1.48, *c]

    # CICLO SUL NUMERO DI TERMINI DEL POLINOMIO
    for n_term in tqdm(range(n_term_min, n_term_max), desc = 'K2_g_fit_no_corrections'):
        inf_bounds = -np.ones(n_term)*100
        sup_bounds = np.ones(n_term)*100
        popt, pcov = curve_fit(poly, (sorted_K,sorted_L), sorted_K2_g, p0 = initParam, sigma = sorted_err_K2_g, absolute_sigma = True, bounds = ((0.08, 1.40, 2.8, *inf_bounds), (0.1, 1.60, 3.2, *sup_bounds)))
        # popt, pcov = curve_fit(poly, (sorted_K,sorted_L), sorted_K2_g, p0 = initParam, sigma = sorted_err_K2_g, absolute_sigma = True, bounds = ((0.05, 1.40, *inf_bounds), (0.1, 1.80, *sup_bounds)))

        initParam = popt
        initParam = np.append(initParam, 0.0)

        chiq = np.sum(((sorted_K2_g - poly((sorted_K,sorted_L), *popt)) / (sorted_err_K2_g))**2)
        ndof = len(sorted_K2_g) - len(initParam)
        red_chiq = chiq/ndof
        perr = np.sqrt(np.diag(pcov))
        Kc = np.append(Kc, popt[0])
        err_Kc = np.append(err_Kc, perr[0])
        yk = np.append(yk, popt[1])
        err_yk = np.append(err_yk, perr[1])
        theta2 = np.append(theta2, popt[2])
        err_theta2 = np.append(err_theta2, perr[2])
        term = np.append(term, n_term)
        red_chisq = np.append(red_chisq, red_chiq)
        coeff.append(popt[3:])

    # PARAMETRI OTTIMALI DA RESTITUIRE
    delta = np.abs(1-red_chisq)
    c = np.array(coeff[np.argmin(delta)])
    popt = [Kc[np.argmin(delta)], yk[np.argmin(delta)], theta2[np.argmin(delta)], *c]
    err_opt = [err_Kc[np.argmin(delta)], err_yk[np.argmin(delta)], err_theta2[np.argmin(delta)]]
    n_term_K2_g = int(term[np.argmin(delta)])
    red_chisq_opt = red_chisq[np.argmin(delta)]

    return popt, err_opt, n_term_K2_g, red_chisq_opt

# DEFINISCO LA FUNZIONE DI ANALISI CON CORREZIONI
def K2_analysis_corrections(k_min, k_max, l_min, n_term_min, n_term_max, omega_min, omega_max, omega_step, eps, kc_init, y_init, shift):
    def poly_corrections(X, a, b, d, *c):
        k,l = X
        poly = 0
        for i in range(n_term1):
            poly = poly + c[i] * ((k-a)*l**(b))**i
        for i in range(n_term2):
            poly = poly + (c[i+n_term1] * ((k-a)*l**(b))**i)*(l**(-omega))
        return (poly)*l**d

    # def poly_corrections(X, a, b, *c):
    #     k,l = X
    #     poly = 0
    #     poly2 = 0
    #     for i in range(n_term1):
    #         poly = poly + c[i+n_term3] * ((k-a)*l**(b))**i
    #     for i in range(n_term2):
    #         poly = poly + (c[i+n_term1+n_term3] * ((k-a)*l**(b))**i)*(l**(-omega))
    #     for i in range(n_term3):
    #         poly2 = poly2 + c[i]*(k)**i
    #     return (poly)*l**(2*b) + poly2*l**3

    # APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
    L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len, K2_g, err_K2_g, K2_sp, err_K2_sp, K3_g, err_K3_g, K3_sp, err_K3_sp = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

    # TAGLIO LE MISURE ENTRO L'INTERVALLO DESIDERATO
    mask = ((K>=k_min) & (K<=k_max) & (L>=l_min))
    K = K[mask]
    L = L[mask]
    K2_g = K2_g[mask]
    err_K2_g = err_K2_g[mask]

    # RIORDINO IL FILE
    sorted_K = K[K.argsort()]
    sorted_L = L[K.argsort()]
    sorted_K2_g = K2_g[K.argsort()]
    sorted_err_K2_g = err_K2_g[K.argsort()]

    Kc = np.array([])
    err_Kc = np.array([])
    yk = np.array([])
    err_yk = np.array([])
    theta2 = np.array([])
    err_theta2 = np.array([])
    red_chisq = np.array([])
    term1 = np.array([])
    term2 = np.array([])
    term3 = np.array([])
    omg = np.array([])
    coeff = []

    c = np.zeros(3*n_term_min - 2*shift, dtype = 'float')
    initParam = [kc_init, y_init, y_init*2, *c]
    # initParam = [kc_init, y_init, *c]
    for omega in tqdm(np.arange(omega_min, omega_max, omega_step), desc = 'K2_g_fit_w_corrections'):
        for n_term1 in range(n_term_min, n_term_max):
            for n_term2 in range(n_term_min-shift, n_term_max-shift):
                for n_term3 in range(n_term_min-shift, n_term_max-shift):
                    print(omega, n_term1, n_term2, n_term3)

                    inf_bounds = -np.ones(n_term3 + n_term1 + n_term2)*100
                    sup_bounds = np.ones(n_term3 + n_term1 + n_term2)*100

                    popt, pcov = curve_fit(poly_corrections, (sorted_K, sorted_L), sorted_K2_g, p0 = initParam, sigma = sorted_err_K2_g, absolute_sigma = True, bounds = ((kc_init-0.005, 1.4, 2.8, *inf_bounds), (kc_init+0.005, 1.6, 4.0, *sup_bounds)))

                    initParam = popt
                    initParam = np.append(initParam, 0.0)

                    chiq = np.sum(((sorted_K2_g - poly_corrections((sorted_K,sorted_L), *popt)) / (sorted_err_K2_g))**2)
                    ndof = len(sorted_K2_g) - len(initParam)
                    red_chiq = chiq/ndof

                    Kc = np.append(Kc, popt[0])
                    term1 = np.append(term1, n_term1)
                    term2 = np.append(term2, n_term2)
                    term3 = np.append(term3, n_term3)
                    yk = np.append(yk, popt[1])
                    theta2 = np.append(theta2, popt[2])

                    perr = np.sqrt(np.diag(pcov))
                    err_Kc = np.append(err_Kc, perr[0])
                    err_yk = np.append(err_yk, perr[1])
                    err_theta2 = np.append(err_theta2, perr[2])

                    red_chisq = np.append(red_chisq, red_chiq)
                    omg = np.append(omg, omega)
                    coeff.append(popt[3:])
                    # coeff.append(popt[2:])

                c = np.zeros(n_term2 + n_term1 + n_term_min - shift +1, dtype = 'float')
                initParam = [kc_init, y_init, y_init*2, *c]
                # initParam = [kc_init, y_init, *c]

            c = np.zeros(n_term_min + n_term_min + n_term1 - 2*shift + 1, dtype = 'float')
            initParam = [kc_init, y_init, y_init*2, *c]
            # initParam = [kc_init, y_init, *c]

        c = np.zeros(3*n_term_min - 2*shift, dtype = 'float')
        initParam = [kc_init, y_init, y_init*2, *c]
        # initParam = [kc_init, y_init, *c]

    # PARAMETRI OTTIMALI DA RESTITUIRE
    # NEI PASSI SEGUENTI TENGO CONTO
    # DELLE SISTEMATICHE, CIOÈ TENENDO IL CHI^2
    # ENTRO UN CERTO INTERVALLO GUARDO COME VARIANO
    # I PARAMETRI PER DIVERSI GRADI, OMEGA ..
    delta = np.abs(1-red_chisq)
    mask = ((delta<=eps))
    red_chisq_masked = red_chisq[mask]
    Kc_masked = Kc[mask]
    err_Kc_masked = err_Kc[mask]
    yk_masked = yk[mask]
    err_yk_masked = err_yk[mask]
    theta2_masked = theta2[mask]
    err_theta2_masked = err_theta2[mask]

    mean_chisq_red = np.mean(red_chisq_masked)

    mean_Kc = np.mean(Kc_masked)
    std_Kc = (np.max(Kc_masked)-np.min(Kc_masked))/2.0

    mean_yk = np.mean(yk_masked)
    std_yk = (np.max(yk_masked)-np.min(yk_masked))/2.0

    mean_theta2 = np.mean(theta2_masked)
    std_theta2 = (np.max(theta2_masked) - np.min(theta2_masked))/2.0

    # PARAMETRI PER PLOT, DA NON PRENDERE TROPPO SUL SERIO
    cc = np.array(coeff[np.argmin(delta)])
    plot_params = [Kc[np.argmin(delta)], yk[np.argmin(delta)], *cc]
    n_term1_K2_g = int(term1[np.argmin(delta)])
    n_term2_K2_g = int(term2[np.argmin(delta)])
    n_term3_K2_g = int(term3[np.argmin(delta)])
    omega_plot = omg[np.argmin(delta)]

    pl.errorbar(sorted_K, sorted_K2_g, sorted_err_K2_g, ls = '', marker = 'o', fillstyle='none')
    n_term1 = n_term1_K2_g
    n_term2 = n_term2_K2_g
    n_term3 = n_term3_K2_g
    omega = omega_plot
    xx = np.linspace(min(sorted_K), max(sorted_K),500)
    yy = np.array([sorted_L[0]]*len(xx))
    pl.plot(xx, poly_corrections((xx, yy), *plot_params))

    return plot_params, n_term1_K2_g, n_term2_K2_g, n_term3_K2_g, omega_plot, mean_Kc, std_Kc, mean_yk, std_yk, mean_chisq_red, mean_theta2, std_theta2
