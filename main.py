<<<<<<< HEAD
import numpy as np
import scipy.constants as const
from utils import *
import math


### GLOBAL VARIABLES ###
GAMMA = 2 * np.pi * 6e6
GAMMA_g = GAMMA / 3
GAMMA_g1 = 2 * GAMMA / 3

DELTA = 2 * np.pi * 740e6
delta = 0.1e6

OMEGA_b = 2 * np.pi * 30e6
OMEGA_r = 2 * np.pi * 100e6

### BASIS STATES ###
g = np.zeros((4, 1), dtype=complex)
g[0, 0] = 1

r = np.zeros((4, 1), dtype=complex)
r[1, 0] = 1

p = np.zeros((4, 1), dtype=complex)
p[2, 0] = 1

g1 = np.zeros((4, 1), dtype=complex)
g1[3, 0] = 1

### HAMILTONIAN ###
H = OMEGA_r / 2 * (g @ h_conj(p) + p @ h_conj(g)) + \
    OMEGA_b / 2 * (p @ h_conj(r) + r @ h_conj(p)) - \
        DELTA * (p * h_conj(p)) - delta * (r @ h_conj(r))


def ham():
    return H


def func(t, y):
    rho = np.reshape(y, (4,4))
    L = GAMMA_g / 2 * (2 * g @ h_conj(p) @ rho @ p @ h_conj(g) - \
                p @ h_conj(p) @ rho - rho @ p @ h_conj(p)) + \
        GAMMA_g1 / 2 * (2 * g1 @ h_conj(p) @ rho @ p @ h_conj(g1) -
                p @ h_conj(p) @ rho - rho @ p @ h_conj(p))

    dfdt =  (H @ rho - rho @ H) / (1j) + L
    dfdt_fl = np.matrix.flatten(dfdt)
    # print(dfdt)
    return dfdt_fl
=======
import numpy as np
import scipy.constants as const
from utils import *
import math


### GLOBAL VARIABLES ###
GAMMA = 2 * np.pi * 6e6
GAMMA_g = GAMMA / 3
GAMMA_g1 = 2 * GAMMA / 3

DELTA = 2 * np.pi * 740e6
delta = 0.1e6

OMEGA_b = 2 * np.pi * 100e6
OMEGA_r = 2 * np.pi * 100e6

### BASIS STATES ###
g = np.zeros((4, 1), dtype=complex)
g[0, 0] = 1

r = np.zeros((4, 1), dtype=complex)
r[1, 0] = 1

p = np.zeros((4, 1), dtype=complex)
p[2, 0] = 1

g1 = np.zeros((4, 1), dtype=complex)
g1[3, 0] = 1

### HAMILTONIAN ###
H = OMEGA_r / 2 * (g @ h_conj(p) + p @ h_conj(g)) + \
    OMEGA_b / 2 * (p @ h_conj(r) + r @ h_conj(p)) - \
        DELTA * (p * h_conj(p)) - delta * (r @ h_conj(r))


def ham():
    return H


def func(t, y):
    rho = np.reshape(y, (4,4))
    L = GAMMA_g / 2 * (2 * g @ h_conj(p) @ rho @ p @ h_conj(g) - \
                p @ h_conj(p) @ rho - rho @ p @ h_conj(p)) + \
        GAMMA_g1 / 2 * (2 * g1 @ h_conj(p) @ rho @ p @ h_conj(g1) -
                p @ h_conj(p) @ rho - rho @ p @ h_conj(p))

    dfdt =  (H @ rho - rho @ H) / (1j) + L
    dfdt_fl = np.matrix.flatten(dfdt)
    # print(dfdt)
    return dfdt_fl
>>>>>>> 6b24c187b39760fbbf057dde81250bdfb2f0bb43
