

def calc_old():
    for j in range(iters):
        for ss in range(len(rho_2)):
            I_r[ss] = simpson(func_2(rho_2[ss], rho_1, k, l,
                              g_1, R_i_left[:, j], n, tilt_delta, a), rho_1)

        R_i_right[:, j] = 1j**(n+1) * k / l * \
            np.exp(-1j*k / (2*l)*g_2*rho_2**2) * I_r
        gamma_i_right[j] = simpson(np.abs(R_i_right[:, j])**2, rho_2)
        R_i_right[:, j] /= np.max(np.abs(R_i_right[:, j]))

        for ii in range(len(rho_1)):
            I_l[ii] = simpson(func_1(rho_2, rho_1[ii], k, l,
                              g_2, R_i_right[:, j], n), rho_2)

        R_i_left[:, j+1] = 1j**(n+1) * k / l * \
            np.exp(-1j*k/(2*l)*g_1*(rho_1+d_r)**2) * I_l
        gamma_i_left[j] = simpson(np.abs(R_i_left[:, j+1])**2, rho_1)
        R_i_left[:, j+1] /= np.max(np.abs(R_i_left[:, j+1]))
    return R_i_right, R_i_left, gamma_i_right, gamma_i_left
