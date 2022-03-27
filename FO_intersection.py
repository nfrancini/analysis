import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from tqdm import tqdm
import time

def FO_fit(x, y, dy, j_min, j_max, n_term1, n_term2):
    def poly(z, *c):
        poly = 0
        for i in range(n_term):
            poly = poly + c[i] * z**i
        return poly

    # TAGLIO LE MISURE ENTRO L'INTERVALLO DESIDERATO
    mask = ((x>=j_min) & (x<=j_max))
    x = x[mask]
    y = y[mask]
    dy = dy[mask]
    cc = []
    red_chisq = np.array([])
    term = np.array([])

    initParam = np.ones(n_term1)
    for n_term in tqdm(range(n_term1, n_term2)):
        popt, pcov = curve_fit(poly, x, y, p0=initParam, sigma = dy, absolute_sigma=True)

        initParam = popt
        initParam = np.append(initParam, 1.0)

        chiq = np.sum(( (y - poly(x, *popt)) / dy )**2)
        ndof = len(y) - len(initParam)
        red_chiq = chiq/ndof
        red_chisq = np.append(red_chisq, red_chiq)
        term = np.append(term, n_term)
        cc.append(popt[0:])

    delta = np.abs(1-red_chisq)
    p = cc[np.argmin(delta)]
    red_chisq_opt = red_chisq[np.argmin(delta)]
    term_opt = term[np.argmin(delta)]

    return p, red_chisq_opt, term_opt

def FO_intersection(path_filenames, j_min, j_max, n_term_min, n_term_max, inter_min, inter_max, ctrl_plot):
    coeffs = []

    # TROVO I VARI COEFFICIENTI PER OGNI TAGLIA
    # POI CERCO LA RADICE A DUE A DUE
    n_term1 = n_term_min
    n_term2 = n_term_max
    for fname in path_filenames:
        with open(fname) as infile:
            aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2, aux_err_K2, aux_K3, aux_err_K3, aux_K4, aux_err_K4 = np.genfromtxt(infile, delimiter ="\t", unpack = True)

            sorted_J = aux_J[aux_J.argsort()]
            sorted_U = aux_U[aux_J.argsort()]
            sorted_err_U = aux_err_U[aux_J.argsort()]

            p, red_chisq_opt_FO, n_term_opt = FO_fit(sorted_J, sorted_U, sorted_err_U, j_min, j_max, n_term1, n_term2)
            coeffs.append(p)

            n_term1 = int(n_term_opt)
            n_term2 = int(n_term_opt+1)

            # FARE PLOT DI CONTROLLO
            if(ctrl_plot == True):
                curve = np.polynomial.Polynomial(p)
                x = np.linspace(min(sorted_J), max(sorted_J), 500)
                pl.errorbar(sorted_J, sorted_U, sorted_err_U, ls = '', fillstyle = 'none', label = 'L='+str(int(aux_L[0])))
                pl.plot(x, curve(x))
                pl.legend()
                pl.show

    intersection = np.array([])

    for i in range(len(coeffs)-1):
        print(i)
        roots = np.polynomial.polynomial.Polynomial.roots(np.polynomial.Polynomial(coeffs[i+1]) - np.polynomial.Polynomial(coeffs[i]))
        roots = np.real(roots[np.isreal(roots)])
        mask = ((roots>=inter_min) & (roots<=inter_max))
        masked_roots = roots[mask]
        intersection = np.append(intersection, masked_roots)

    mean_intersect = np.mean(intersection)
    err_intersect = (max(intersection) - min(intersection)) / 2.0

    return mean_intersect, err_intersect, n_term_opt
