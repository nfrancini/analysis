import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from tqdm import tqdm

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

    initParam = np.zeros(n_term1)
    for n_term in tqdm(range(n_term1, n_term2)):
        popt, pcov = curve_fit(poly, x, y, p0=initParam, sigma = dy, absolute_sigma=True)
        # popt, stats = np.polynomial.polynomial.polyfit(x, y, n_term-1, full = True, w = 1.0/(dy*100))
        # initParam = popt//
        initParam =  np.zeros(n_term+1)

        chiq = np.sum(( (y - poly(x, *popt)) / dy )**2)
        ndof = len(y) - len(initParam)
        red_chiq = chiq/ndof
        # red_chiq = stats[0]/ndof
        red_chisq = np.append(red_chisq, red_chiq)
        term = np.append(term, n_term)
        cc.append(popt[0:])

    delta = np.abs(1-red_chisq)
    p = cc[np.argmin(delta)]
    red_chisq_opt = red_chisq[np.argmin(delta)]
    term_opt = term[np.argmin(delta)]

    return p, red_chisq_opt, term_opt

def FO_intersection(path_filenames, j_min, j_max, l_min, l_max, n_term_min, n_term_max, inter_min, inter_max, ctrl_plot):
    coeffs = []

    # TROVO I VARI COEFFICIENTI PER OGNI TAGLIA
    # POI CERCO LA RADICE A DUE A DUE
    n_term1 = n_term_min
    n_term2 = n_term_max
    for fname in path_filenames:
        with open(fname) as infile:
            aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)

            if((aux_L[0]>=l_min) & (aux_L[0]<=l_max)):
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
                    x = np.linspace(j_min, j_max, 500)
                    pl.errorbar(sorted_J, sorted_U, sorted_err_U, ls = '', fillstyle = 'none', marker = 'o', label = 'L='+str(int(aux_L[0])))
                    pl.plot(x, curve(x))
                    pl.xlim(j_min, j_max)
                    pl.legend()

    if(ctrl_plot == True):
        pl.show()

    intersection = np.array([])

    for i in range(len(coeffs)-1):
        roots = np.polynomial.polynomial.Polynomial.roots(np.polynomial.Polynomial(coeffs[i+1]) - np.polynomial.Polynomial(coeffs[i]))
        roots = np.real(roots[np.isreal(roots)])
        mask = ((roots>=inter_min) & (roots<=inter_max))
        masked_roots = roots[mask]
        intersection = np.append(intersection, masked_roots)

    mean_intersect = np.mean(intersection)
    err_intersect = (max(intersection) - min(intersection)) / 2.0

    return mean_intersect, err_intersect, n_term_opt

# def FO_quadratic(path_filenames, j_min, j_max, l_min, l_max, ctrl_plot):
#     def quadratic(x, a, b, c):
#         return a*(x-b)**2 + c
#
#     for fname in path_filenames:
#         with open(fname) as infile:
#             aux_L, aux_J, aux_K, aux_ene_sp, aux_err_ene_sp, aux_ene_g, aux_err_ene_g, aux_ene_dens, aux_err_ene_dens, aux_susc, aux_err_susc, aux_G_pm, aux_err_G_pm, aux_C, aux_err_C, aux_U, aux_err_U, aux_corr_len, aux_err_corr_len, aux_K2_g, aux_err_K2_g, aux_K2_sp, aux_err_K2_sp, aux_K3_g, aux_err_K3_g, aux_K3_sp, aux_err_K3_sp = np.genfromtxt(infile, delimiter ="\t", unpack = True)
#             # print(aux_C, len(aux_C))
#             # print(aux_err_C, len(aux_err_C))
#             if((aux_L[0]>=l_min) & (aux_L[0]<=l_max)):
#                 sorted_J = aux_J[aux_J.argsort()]
#                 sorted_U = aux_U[aux_J.argsort()]
#                 sorted_err_U = aux_err_U[aux_J.argsort()]
#                 sorted_C = aux_C[aux_J.argsort()]
#                 sorted_err_C = aux_err_C[aux_J.argsort()]
#
#                 mask = ((sorted_J >= j_min) & (sorted_J <= j_max))
#                 sorted_J = sorted_J[mask]
#                 sorted_U = sorted_U[mask]
#                 sorted_C = sorted_C[mask]
#                 sorted_err_U = sorted_err_U[mask]
#                 sorted_err_C = sorted_err_C[mask]
#
#                 # TROVO IL MASSIMO E LO METTO NEI PARAMETRI INIZIALI
#                 C_crit = np.max(sorted_C)
#                 J_crit = sorted_J[np.argmax(sorted_C)]
#                 a_init = (sorted_C[0] - C_crit) / (sorted_J[0] - J_crit)**2
#
#                 initParam = [a_init, J_crit, C_crit]
#
#                 # PRENDO SOLO I SEI PUNTI PRIMA E DOPO IL PICCO
#                 cut_J = sorted_J[np.argmax(sorted_C)-1 : np.argmax(sorted_C) + 2]
#                 cut_C = sorted_C[np.argmax(sorted_C)-1 : np.argmax(sorted_C) + 2]
#                 cut_err_C = sorted_err_C[np.argmax(sorted_C)-1 : np.argmax(sorted_C) + 2]
#
#                 # PROVO IL FIT QUADRATICO
#                 popt, pcov = curve_fit(quadratic, cut_J, cut_C, p0 = initParam, sigma = cut_err_C, absolute_sigma = True)
#                 chiq = np.sum(( (cut_C - quadratic(cut_J, *popt)) / cut_err_C )**2)
#                 ndof = len(cut_C) - len(initParam)
#                 perr = np.sqrt(np.diag(pcov))
#                 red_chiq = chiq/ndof
#                 print(red_chiq)
#                 print("a = %f +- %f" % (popt[0], perr[0]))
#                 print("b = %f +- %f" % (popt[1], perr[1]))
#                 print("c = %f +- %f" % (popt[2], perr[2]))
#
#                 # FARE PLOT DI CONTROLLO
#                 if(ctrl_plot == True):
#                     x = np.linspace(np.min(cut_J), np.max(cut_J), 500)
#                     y = quadratic(x, *popt)
#                     pl.errorbar(cut_J, cut_C, cut_err_C, ls = '', fillstyle = 'none', marker = 'o', label = 'L='+str(int(aux_L[0])))
#                     pl.plot(x,y)
#                     pl.legend()
#     return popt
