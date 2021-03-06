import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from scipy.odr import odrpack
from tqdm import tqdm


# DEFINISCO LA FUNZIONE DI ANALISI
def susc_analysis(j_min, j_max, l_min, n_term_min, n_term_max, eta_init):
    # DEFINISCO LA FUNZIONE DI FIT
    def poly_odr(pars, x):
        poly = 0
        for i in range(n_term):
            poly = poly + pars[1+i] * (x[0])**i
        return (x[1]**(2-pars[0]))*poly

    # APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
    L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len, K2_g, err_K2_g, K2_sp, err_K2_sp, K3_g, err_K3_g, K3_sp, err_K3_sp = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

    # PROVO IL FIT PER LA SUSCETTIVITÃ€, USO ODR PER TENERE IN CONTO
    # DEGLI ERRORI SULLA LUNGHEZZA DI CORRELAZIONE CHE FA DA VARIABILE INDIPENDENTE
    mask = ((J>=j_min) & (J<=j_max) & (L>=l_min))
    corr_len = corr_len[mask]
    L = L[mask]
    susc = susc[mask]
    err_susc = err_susc[mask]
    err_corr_len = err_corr_len[mask]
    J = J[mask]

    sorted_L = L[J.argsort()]
    sorted_corr_len = corr_len[J.argsort()]
    sorted_err_corr_len = err_corr_len[J.argsort()]
    sorted_susc = susc[J.argsort()]
    sorted_err_susc = err_susc[J.argsort()]

    # DEFINISCO ARRAY DI SUPPORTO PER I RISULTATI E PER
    # FISSARE I VARI ERRORI
    machineEps = np.finfo('float').eps
    realEps = np.sqrt(machineEps)

    eta = np.array([])
    err_eta = np.array([])
    red_chisq = np.array([])
    term = np.array([])
    coeff = []

    a = np.full_like(sorted_corr_len, 1, dtype = 'int')
    b = np.full_like(sorted_L, 0, dtype = 'int')

    # CICLO SUL NUMERO DI TERMINI DEL POLINIMIO
    for n_term in tqdm(range(n_term_min, n_term_max), desc = 'susc_fit_no_corrections'):
        c = np.ones(n_term)
        init_odr = [eta_init, *c]
        x = np.row_stack( (sorted_corr_len/sorted_L, sorted_L) )
        sigma_x = np.row_stack( (sorted_err_corr_len/sorted_L, sorted_L*realEps)) # ERRORI FITTIZI SULLA L

        fixed = np.row_stack((a , b) )          # USO QUESTO PER FISSARE GLI ERRORI DI L A ZERO

        model = odrpack.Model(poly_odr)
        data = odrpack.RealData(x, sorted_susc, sx = sigma_x, sy = sorted_err_susc, fix = fixed)
        odr = odrpack.ODR(data, model, beta0 = init_odr, ifixx = fixed)
        # odr.set_iprint(init = 2, final = 2)
        out = odr.run()

        ndof = len(sorted_susc) - len(init_odr)
        red_chisq = np.append(red_chisq, (out.sum_square)/ndof)
        eta = np.append(eta, out.beta[0])
        err_eta = np.append(err_eta, np.sqrt(out.cov_beta.diagonal())[0])
        term = np.append(term, n_term)
        coeff.append(out.beta[1:])

    # GRANDEZZE DA RESTITUIRE
    delta = np.abs(1-red_chisq)
    n_term_susc = int(term[np.argmin(delta)])
    c = np.array(coeff[np.argmin(delta)])
    popt = [eta[np.argmin(delta)], *c]
    err_opt = err_eta[np.argmin(delta)]
    red_chisq_opt = red_chisq[np.argmin(delta)]

    return popt, err_opt, n_term_susc, red_chisq_opt

def susc_analysis_corrections(j_min, j_max, l_min, n_term_min, n_term_max, shift, omega_min, omega_max, omega_step, eps, eta_init):
    # DEFINISCO LA FUNZIONE DI FIT
    def poly_odr_corrections(pars, x):
        poly = 0
        for i in range(n_term1):
            poly = poly + pars[1+i] * (x[0])**i
        for i in range(n_term2):
            poly = poly + (x[1]**(-omega)) * pars[1+i+n_term1] * (x[0])**i
        return (x[1]**(2-pars[0]))*poly

    # APRO I VARI FILE DA ANALIZZARE E CREO UN UNICO BLOCCO
    L, J, K, ene_sp, err_ene_sp, ene_g, err_ene_g, ene_dens, err_ene_dens, susc, err_susc, G_pm, err_G_pm, C, err_C, U, err_U, corr_len, err_corr_len, K2_g, err_K2_g, K2_sp, err_K2_sp, K3_g, err_K3_g, K3_sp, err_K3_sp = np.genfromtxt("./temp.dat", delimiter ="\t", unpack = True)

    # PROVO IL FIT PER LA SUSCETTIVITÃ€, USO ODR PER TENERE IN CONTO
    # DEGLI ERRORI SULLA LUNGHEZZA DI CORRELAZIONE CHE FA DA VARIABILE INDIPENDENTE
    mask = ((J>=j_min) & (J<=j_max) & (L>=l_min))
    corr_len = corr_len[mask]
    L = L[mask]
    susc = susc[mask]
    err_susc = err_susc[mask]
    err_corr_len = err_corr_len[mask]
    J = J[mask]

    sorted_L = L[J.argsort()]
    sorted_corr_len = corr_len[J.argsort()]
    sorted_err_corr_len = err_corr_len[J.argsort()]
    sorted_susc = susc[J.argsort()]
    sorted_err_susc = err_susc[J.argsort()]

    machineEps = np.finfo('float').eps
    realEps = np.sqrt(machineEps)

    eta = np.array([])
    err_eta = np.array([])
    red_chisq = np.array([])
    term1 = np.array([])
    term2 = np.array([])
    omg = np.array([])
    coeff = []

    a = np.full_like(sorted_corr_len, 1, dtype = 'int')
    b = np.full_like(sorted_L, 0, dtype = 'int')

    for omega in tqdm(np.arange(omega_min, omega_max, omega_step), desc = 'susc_fit_w_corrections'):
        for n_term1 in range(n_term_min, n_term_max):
            for n_term2 in range(n_term_min -shift, n_term_max-shift):
                c = np.zeros(n_term1 + n_term2)
                init_odr = [eta_init, *c]
                fixed = np.row_stack((a , b) )          # USO QUESTO PER FISSARE GLI ERRORI DI L A ZERO
                model = odrpack.Model(poly_odr_corrections)
                x = np.row_stack((sorted_corr_len/sorted_L, sorted_L))
                sigma_x = np.row_stack((sorted_err_corr_len/sorted_L, sorted_L*realEps))
                data = odrpack.RealData(x, sorted_susc, sx = sigma_x, sy = sorted_err_susc, fix = fixed)
                odr = odrpack.ODR(data, model, beta0 = init_odr, ifixx = fixed)
                out = odr.run()
                ndof = len(sorted_susc) - len(init_odr)
                red_chisq = np.append(red_chisq, (out.sum_square)/ndof)
                eta = np.append(eta, out.beta[0])
                err_eta = np.append(err_eta, np.sqrt(out.cov_beta.diagonal())[0])
                term1 = np.append(term1, n_term1)
                term2 = np.append(term2, n_term2)
                coeff.append(out.beta[1:])
                omg = np.append(omg, omega)

    # GRANDEZZE DA RESTITUIRE
    delta = np.abs(1 - red_chisq)
    mask = ((delta<=eps))
    red_chisq_masked = red_chisq[mask]
    eta_masked = eta[mask]
    err_eta_masked = err_eta[mask]
    omega_masked = omg[mask]

    mean_chisq_red = np.mean(red_chisq_masked)

    mean_eta = np.mean(eta_masked)
    std_eta = (np.max(eta_masked) - np.min(eta_masked))/2.0

    omega_opt = np.mean(omega_masked)
    std_omega = (np.max(omega_masked)-np.min(omega_masked))/2.0

    # PARAMETRI PER PLOT, DA NON PRENDERE TROPPO SUL SERIO
    n_term1_susc = term1[np.argmin(delta)]
    n_term2_susc = term2[np.argmin(delta)]
    cc = np.array(coeff[np.argmin(delta)])
    plot_params = [eta[np.argmin(delta)], *cc]
    omega_plot = omg[np.argmin(delta)]

    return plot_params, n_term1_susc, n_term2_susc, omega_plot, omega_opt, std_omega, mean_eta, std_eta, mean_chisq_red
