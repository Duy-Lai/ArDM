# quenching functions, written by Alex Stauffer
# adapted by Duy Lai, 09.2022

import numpy as np

# Lindhard nuclear quenching
def f_Lindhard_Nuclear_Quenching(E_NR):
    # Parameters
    Ar_A = 40  # Atomic weight A, unit-less
    Ar_Z = 18  # Atomic number Z, unit-less
    Xe_A = 131  # Atomic weight A, unit-less
    Xe_Z = 54  # Atomic number Z, unit-less
    # unit-less efficiency of transferring the kinetic energy of NRs to electric energy -> S1L
    k_Ar = 0.133 * Ar_Z ** (2 / 3) * Ar_A ** (-1 / 2)  # for Ar specifically; unit-less
    k_Xe = 0.133 * Xe_Z ** (2 / 3) * Xe_A ** (-1 / 2)  # for Xe specifically; unit-less
    epsilon_Ar = 11.5 * E_NR * Ar_Z ** (-7 / 3)  # E_NR in keV -> Epsilon unit-less
    epsilon_Xe = 11.5 * E_NR * Xe_Z ** (-7 / 3)  # E_NR in keV -> Epsilon unit-less
    del_epsilon_del_E_nr_Ar = 11.5 * Ar_Z ** (-7 / 3)  # [1/keV]
    del_epsilon_del_E_nr_Xe = 11.5 * Xe_Z ** (-7 / 3)  # [1/keV]
    g_of_epsilon_Ar = 3 * epsilon_Ar ** 0.15 + 0.7 * epsilon_Ar ** 0.6 + epsilon_Ar  # unit-less
    g_of_epsilon_Xe = 3 * epsilon_Xe ** 0.15 + 0.7 * epsilon_Xe ** 0.6 + epsilon_Xe  # unit-less
    del_g_del_epsilon_Ar = 3 * 0.15 * epsilon_Ar ** (-0.85) + 0.7 * 0.6 * epsilon_Ar ** (-0.4) + 1  # unit-less
    del_g_del_epsilon_Xe = 3 * 0.15 * epsilon_Xe ** (-0.85) + 0.7 * 0.6 * epsilon_Xe ** (-0.4) + 1  # unit-less
    q_nc_Ar = k_Ar * g_of_epsilon_Ar / (1 + k_Ar * g_of_epsilon_Ar)  # unit-less
    q_nc_Xe = k_Xe * g_of_epsilon_Xe / (1 + k_Xe * g_of_epsilon_Xe)  # unit-less
    del_q_nc_del_g_Ar = k_Ar / ((1 + k_Ar * g_of_epsilon_Ar) ** 2)  # unit-less
    del_q_nc_del_g_Xe = k_Xe / ((1 + k_Xe * g_of_epsilon_Xe) ** 2)  # unit-less
    del_q_nc_del_E_nr_Ar = del_q_nc_del_g_Ar * del_g_del_epsilon_Ar * del_epsilon_del_E_nr_Ar  # [1/keV]
    del_q_nc_del_E_nr_Xe = del_q_nc_del_g_Xe * del_g_del_epsilon_Xe * del_epsilon_del_E_nr_Xe  # [1/keV]
    # del_E_ee_del_E_nr_Ar = q_nc_Ar * (1 + E_NR / q_nc_Ar * del_q_nc_del_E_nr_Ar)  # unit-less
    # del_E_ee_del_E_nr_Xe = q_nc_Xe * (1 + E_NR / q_nc_Xe * del_q_nc_del_E_nr_Xe)  # unit-less
    return [q_nc_Ar, q_nc_Xe, del_q_nc_del_E_nr_Ar, del_q_nc_del_E_nr_Xe]

q_bi_excitonic_Ar = 0.6  # unit-less, Hitachi2019

# Drift field recombination quenching
def f_Drift_Field_Recombination_Quenching(E_NR, F_drift, L_eff, del_L_eff_del_E_nr=np.empty(0)):
    # Kimura2019 & Washimi2018 (for S1L)
    # F_drift [V/cm]
    # E_NR [keV]
    # L_eff = Lindhard * Bi-excitonic [unit-less]
    alpha_0 = 1.0  # unit-less
    D_alpha = 8.9e-4  # [cm/V]
    alpha = alpha_0 * np.exp(-D_alpha * F_drift)  # unit-less
    gamma = 1.15  # [(V/cm) ** delta]
    delta = 0.576  # unit-less
    eta = gamma * (F_drift ** (-delta))  # unit-less
    W_value_Ar = 19.5  # [eV]
    N_i = (E_NR * 1000) / W_value_Ar * 1 / (alpha + 1) * L_eff  # unit-less
    Rstar = 1 - np.log(1 + N_i * eta) / (N_i * eta)  # unit-less
    q_rec_Ar = (alpha + Rstar) / (alpha + 1)  # unit-less
    if del_L_eff_del_E_nr.any():
        del_N_i_del_E_nr = L_eff / (W_value_Ar / 1000.0) / (alpha + 1) * (1 + E_NR / L_eff * del_L_eff_del_E_nr)  # [1/keV]
        del_Rstar_del_N_i = np.log(1 + N_i * eta) / N_i ** 2 / eta - 1 / (N_i + N_i ** 2 * eta)  # unit-less
        del_Rstar_del_E_nr = del_Rstar_del_N_i * del_N_i_del_E_nr  # [1/keV]
        del_q_rec_del_E_nr_Ar = 1 / (alpha + 1) * del_Rstar_del_E_nr  # [1/kev]
        return q_rec_Ar, del_q_rec_del_E_nr_Ar  # unit-less, [1/kev]
    else:
        return q_rec_Ar  # unit-less

# converting keV to photoelectrons
def conversion_keV_pe(E_nr):
    light_yield = 1.33
    F_drift = 230.0
    L_eff = f_Lindhard_Nuclear_Quenching(E_nr)[0] * q_bi_excitonic_Ar
    return light_yield * L_eff * E_nr * f_Drift_Field_Recombination_Quenching(E_nr, F_drift, L_eff)