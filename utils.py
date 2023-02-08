import numpy as np
from scipy import constants as const
from scipy import special as sp
import math
from numba import njit


def h_conj(state):
    return np.conj(state.T)


def g(l, r):
    g = 1 - l / r
    # if g == 0:
    #     return 1e-15
    # elif g != 0:
    return g


def N(a, l, lmbda):
    return a**2 / (l * lmbda)


def A(l, r_2):
    return 2 - 4 * l / r_2


def tau_c(l, gamma):
    c = const.c
    return 2 * l / (c * gamma)


def linewidth(l, gamma):
    tau = tau_c(l, gamma)
    return 1 / (2*np.pi*tau)


def eigfreq(g_1, g_2, l, n, m=0):
    arccos = np.arccos(np.sqrt(g_1 * g_2))
    nu = const.c / (2 * l) * (n + (1 + m) / np.pi * arccos)
    return nu


def mode_spacing(g_1, g_2, l):
    nu_0 = eigfreq(g_1, g_2, l, n=0, m=0)
    nu_1 = eigfreq(g_1, g_2, l, n=1, m=0)
    nu_0m = eigfreq(g_1, g_2, l, n=0, m=0)
    nu_1m = eigfreq(g_1, g_2, l, n=0, m=1)
    return nu_1 - nu_0, nu_1m - nu_0m

def misaligment(g_1, g_2, l, lmbda):
    term_1 = (1 + g_1 * g_2) / np.sqrt((1 - g_1 * g_2)**3)
    term_2 = np.abs(g_1 + g_2) / np.sqrt(g_1 * g_2)
    return np.sqrt(np.pi * l / lmbda * term_1 * term_2)


def mirror_1(g_1, g_2, l, lmbda):
    term_1 = np.sqrt(np.sqrt(g_2 / (g_1 * (1 - g_1*g_2))))
    return np.sqrt(l * lmbda / np.pi) * term_1


def mirror_2(g_1, g_2, l, lmbda):
    term_1 = np.sqrt(np.sqrt(g_1 / (g_2 * (1 - g_1*g_2))))
    return np.sqrt(l * lmbda / np.pi) * term_1


def beam_waist(g_1, g_2, l, lmbda):
    term_3 = (g_1*g_2*(1-g_1*g_2)/(g_1+g_2-2*g_1*g_2)**2)**(1/4)
    return np.sqrt(l * lmbda / np.pi) * term_3


def lorenz_PDF(freq, c_freq, gamma):
    spectrum = np.zeros((len(freq), len(c_freq)))
    sp = np.zeros(len(freq))
    for i in range(len(c_freq)):
        spec = gamma**2 / ((freq-c_freq[i])**2 + gamma**2)
        spectrum[:, i] = spec
        sp += spec
    return spectrum, sp


def LG_beam_far(rho, phi, w, l, p):

    C_pl = np.sqrt((2 * math.factorial(p)) /
                   (np.pi * math.factorial(p + abs(l))))
    L_pl = sp.eval_genlaguerre(p, abs(l), 2 * rho**2 / w**2)

    U = C_pl / w * (np.sqrt(2) * rho / w)**(abs(l)) * np.exp(- rho**2 / w**2) * \
        L_pl * np.exp(1j * l * phi)

    return U / np.max(U)


def func_1(rho_2, rho_1, k, L, g_2, R_n, l):
    dI = R_n * sp.jv(l, (k*rho_1*rho_2) / L) * \
        np.exp(-1j*k/(2*L)*g_2*rho_2**2) * rho_2
    return dI


def func_2(rho_2, rho_1, k, L, g_1, R_n, l, tilt_delta, a):
    d_r = L * tilt_delta / a
    dI = sp.jv(l, (k*rho_1*rho_2) / L) * \
        np.exp(-1j*k/(2*L)*g_1*(rho_1+d_r)**2) * R_n * (rho_1)
    return dI


def beam_rayleigh(w0, lmbda):
    z_r = np.pi * w0**2 / lmbda
    ans = 'Beam Rayleigh length is: ' + str(np.round(z_r*1e2, 2)) + ' cm'
    return ans


def beam_waist_at_z(z, w0, lmbda):
    z_r = np.pi * w0**2 / lmbda
    w = w0 * np.sqrt(1 + (z / z_r)**2)
    ans = 'Beam waist at z = ' + str(z*1e2) + ' cm is: ' + str(np.round(w*1e6, 2)) + ' um'
    return ans


def beam_waist_after_lens(f, w0, lmbda):
    f = f*1e-3
    w = lmbda * f / np.pi / w0
    ans = 'Beam waist after lens is: ' + str(np.round(w*1e6,2)) + ' um' 
    return ans


def lens_focal_dist(w, w0, lmbda):
    f = w * w0 * np.pi / lmbda
    ans = 'Lens focal dist is: ' + str(np.round(f*1e3,2)) + ' mm' 
    return ans

### OLD FUNCTION FOR LOSSES CALCULATIONS ###

# def calc_old():
#     for j in range(iters):
#         for ss in range(len(rho_2)):
#             I_r[ss] = simpson(func_2(rho_2[ss], rho_1, k, l,
#                               g_1, R_i_left[:, j], n, tilt_delta, a), rho_1)

#         R_i_right[:, j] = 1j**(n+1) * k / l * \
#             np.exp(-1j*k / (2*l)*g_2*rho_2**2) * I_r
#         gamma_i_right[j] = simpson(np.abs(R_i_right[:, j])**2, rho_2)
#         R_i_right[:, j] /= np.max(np.abs(R_i_right[:, j]))

#         for ii in range(len(rho_1)):
#             I_l[ii] = simpson(func_1(rho_2, rho_1[ii], k, l,
#                               g_2, R_i_right[:, j], n), rho_2)

#         R_i_left[:, j+1] = 1j**(n+1) * k / l * \
#             np.exp(-1j*k/(2*l)*g_1*(rho_1+d_r)**2) * I_l
#         gamma_i_left[j] = simpson(np.abs(R_i_left[:, j+1])**2, rho_1)
#         R_i_left[:, j+1] /= np.max(np.abs(R_i_left[:, j+1]))
#     return R_i_right, R_i_left, gamma_i_right, gamma_i_left
