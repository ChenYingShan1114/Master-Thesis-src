import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import gridspec
from scipy.optimize import root_scalar, curve_fit

um = 1e-6
nm = 1e-9
fs = 1e-15
q_e = 1.602e-19
m_e = 9.109e-31
epsilon_0 = 8.854e-12
mu_0 = 4 * np.pi * 1e-7
hbar = 6.626e-34 / (2 * np.pi)
c = 299792458
eV = 1.602e-19

q_H_in = 95 #169
PATH = './' +str(q_H_in) + '_png'
w_0 = 21.22*np.sqrt(2) * um #24.75
z_f = 0 * um

lambda_0 = 405 * nm
omega_d = 2 * np.pi * c / lambda_0
k_0 = 2 * np.pi / lambda_0
I_p2 = 54.41776 * eV

import matplotlib as mpl
def format(fig):
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['xtick.labelsize'] = 19
    plt.rcParams['ytick.labelsize'] = 19
    plt.rcParams['font.size'] = 19
    plt.rcParams['figure.figsize'] = [5.6*6, 4*3]
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 6
    plt.rcParams['legend.fontsize'] = 15
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.linewidth'] = 1.5

def ax_format(ax, xmaj, xmin, ymaj, ymin):
    ax.xaxis.set_tick_params(which='major', size=5, width=1, direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=3, width=1, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=5, width=1, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=3, width=1, direction='in', right='on')
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(xmaj))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xmin))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ymaj))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ymin))

# analysis
#-------------------------------------------------------------------------
# HHG_DipolePhase:
#    Calculate the intrinsic dipole phase of the HHG.
#
#    HHG_DipolePhase(q_H,E_d,I_p,omega_d);
#
#       q_H: harmonic order
#       E_d: driving electric field amplitude (V/m)
#       I_p: material ionization potential (J)
#       omega_d: laser angular frequency (rad/sec)
#
#    output: [Phi_dipole_l  Phi_dipole_s]
#       Phi_dipole_l: long-trajectory dipole phase (rad)
#       Phi_dipole_s: short-trajectory dipole phase (rad)
#
# Author: Hsu-hsin Chu (2021/7/11)
# MATLAB to python Adaptor: Ying-Shan Chen
#-------------------------------------------------------------------------
def HHG_DipolePhase(q_H,E_d,I_p,omega_d):
    # Find the ionization time t_0_q for a given harmonic order q_H and laser
    # amplitude E_d: t_0_q(q_H,E_d)
    #    long trajectory:  t_0 = 0 ~ 0.05T
    #    short trajectory: t_0 = 0.05T ~ 0.25T
    t_0_q_l = IonizationTime_long(q_H,E_d,I_p,omega_d)   # unit: sec
    t_0_q_s = IonizationTime_short(q_H,E_d,I_p,omega_d)  # unit: sec

    # Find the recombination time t_r_q for a given harmonic order q_H and laser
    # amplitude E_d: t_r_q(q_H,E_d)
    t_r_q_l = RecombinationTime(omega_d,t_0_q_l)  # unit: sec
    t_r_q_s = RecombinationTime(omega_d,t_0_q_s)  # unit: sec

    # Find the dipole phase for a given harmonic order q_H as a function of laser
    # amplitude E_d: Phi_dipole(E_d)
    Phi_dipole_l = q_H*omega_d*t_r_q_l - (I_p/hbar)*(t_r_q_l-t_0_q_l) - (q_e**2*E_d**2)/(2*hbar*m_e*omega_d**2) * ((-1/2/omega_d)*np.sin(omega_d*t_r_q_l)*np.cos(omega_d*t_r_q_l) + (2/omega_d)*np.cos(omega_d*t_r_q_l)*np.sin(omega_d*t_0_q_l) - (3/2/omega_d)*np.sin(omega_d*t_0_q_l)*np.cos(omega_d*t_0_q_l) + (1/2+np.sin(omega_d*t_0_q_l)**2)*(t_r_q_l-t_0_q_l))
    Phi_dipole_s = q_H*omega_d*t_r_q_s - (I_p/hbar)*(t_r_q_s-t_0_q_s) - (q_e**2*E_d**2)/(2*hbar*m_e*omega_d**2) * ((-1/2/omega_d)*np.sin(omega_d*t_r_q_s)*np.cos(omega_d*t_r_q_s) + (2/omega_d)*np.cos(omega_d*t_r_q_s)*np.sin(omega_d*t_0_q_s) - (3/2/omega_d)*np.sin(omega_d*t_0_q_s)*np.cos(omega_d*t_0_q_s) + (1/2+np.sin(omega_d*t_0_q_s)**2)*(t_r_q_s-t_0_q_s))
    return Phi_dipole_l, Phi_dipole_s

#-----------------------------------------------------------------------
#Find the recombination time t_r(omega_d,t_0) (unit: sec)
#-----------------------------------------------------------------------
# electron position x(t_0,t)
def x_fun(t_fs,t_0):
    t_0_fs     = t_0/fs;          # unit: fs
    omega_d_fs = omega_d*fs;      # unit: red/fs
    return np.cos(omega_d_fs*t_fs) - np.cos(omega_d_fs*t_0_fs) + omega_d_fs * np.sin(omega_d_fs*t_0_fs)*(t_fs - t_0_fs)

def RecombinationTime(omega_d,t_0):
    # unit conversion
    t_0_fs     = t_0/fs;          # unit: fs
    omega_d_fs = omega_d*fs;      # unit: red/fs
    T_fs = 2*np.pi/omega_d_fs;       # period, unit: fs

    # recombination time, determined from x(t_0,t_r) = 0.
    sol = root_scalar(x_fun, args = (t_0), bracket = [t_0_fs+0.001,T_fs])
    t_r_fs = sol.root
    return t_r_fs * fs;                       # unit: sec

#-----------------------------------------------------------------------
# Find the ionization time of the "LONG" trajectory as a function of the
# harmonic order q_H and the laser amplitude E_d:
#    t_0_q_long(q_H,E_d)  unit: sec
# long trajectory:  t_0 = 0 ~ 0.05T
# -----------------------------------------------------------------------
# photon energy U(t_0,I) (unit: J) - Error function (unit: eV)
def U_fun(t_0_fs,E_d,I_p):
    omega_d_fs = omega_d*fs
    return (I_p + q_e**2*E_d**2/2/m_e/omega_d**2 * (np.sin(omega_d*RecombinationTime(omega_d,t_0_fs*fs)) - np.sin(omega_d_fs*t_0_fs))**2 - q_H*hbar*omega_d)/eV

def IonizationTime_long(q_H,E_d,I_p,omega_d):

    # angular frequency and period
    omega_d_fs = omega_d*fs      # unit: rad/fs
    T_fs = 2*np.pi/omega_d_fs       # period, unit: fs

    # Find the ionization time t_0_q_fs for q_H-th harmonic
    # long trajectory:  t_0 = 0 ~ 0.05T
    sol = root_scalar(U_fun, args = (E_d,I_p), bracket = [0, 0.05*T_fs])
    t_0_q_fs = sol.root
    return t_0_q_fs * fs

# -----------------------------------------------------------------------
# Find the ionization time of the "SHORT" trajectory as a function of the
# harmonic order q_H and the laser amplitude E_d:
#    t_0_q_short(q_H,E_d)  unit: sec
# short trajectory:  t_0 = 0.05T ~ 0.25T
# -----------------------------------------------------------------------
def IonizationTime_short(q_H,E_d,I_p,omega_d):

    # angular frequency and period
    omega_d_fs = omega_d*fs      # unit: rad/fs
    T_fs = 2*np.pi/omega_d_fs       # period, unit: fs

    # Find the ionization time t_0_q_fs for q_H-th harmonic
    # short trajectory:  t_0 = 0.05T ~ 0.25T
    sol = root_scalar(U_fun, args = (E_d,I_p), bracket = [0.05*T_fs, 0.249*T_fs])
    t_0_q_fs = sol.root
    return t_0_q_fs * fs;                          # unit: sec

# function
def Plasma_frequency(plasma_density):
    return np.sqrt(plasma_density * q_e * q_e / m_e / epsilon_0)
def Refractive_index_of_plasma(plasma_density, frequency):
    return np.sqrt(1 - (Plasma_frequency(plasma_density) / frequency)**2)
def Rayleigh_length(wavelength):
    return np.pi * w_0 * w_0 / wavelength
def Gouy_phase_shift(wavelength, position):
    return -np.arctan((position - z_f) / Rayleigh_length(wavelength))
def Curvature_Radius(wavelength, position):
    return ((position - z_f)**2 + Rayleigh_length(wavelength)**2) / (position - z_f)

# Load data and convert it to SI units
ave_I_inter_all = np.load(PATH + '/npy/i0_2d.npy') 
ne_all = np.load(PATH + '/npy/ne_2d.npy')
ng_all = np.load(PATH + '/npy/n0_2d.npy')
ni_1_all = np.load(PATH + '/npy/n1_2d.npy')
ni_2_all = np.load(PATH + '/npy/n2_2d.npy')
z_all = np.load(PATH + '/npy/z_1d.npy') * 1e-6
r_all = np.load(PATH + '/npy/r_1d.npy') * 1e-6

NZ_total = z_all.size
NR = r_all.size
DZ = (z_all[1] - z_all[0])
DR = (r_all[1] - r_all[0])

num = 31   # 'num' should be an odd number. This parameter determines the range of orders that the code evaluates in the spectrum. 
bandwidth = np.zeros(num)
q_H_list = np.zeros(num)
for ii in range(num):
    q_H = q_H_in + 2 * ii - (num - 1)
    omega_q = q_H * omega_d
    lambda_q = lambda_0 / q_H
    q_H_list[ii] = q_H

    front = 60
    back = 243

    # Calculate plasma dispersion.
    phi_plasma_all = q_H * omega_d / c * np.cumsum((Refractive_index_of_plasma(ne_all, omega_d) - Refractive_index_of_plasma(ne_all, omega_q)) * DZ, axis = 0)
    phi_plasma_all = phi_plasma_all[:, :] - phi_plasma_all[front, :]

    # Calculate Gouy phase.
    phi_gouy_all = np.repeat((q_H * Gouy_phase_shift(lambda_0, z_all) - Gouy_phase_shift(lambda_q, z_all)), NR).reshape(NZ_total, NR)
    phi_gouy_all = phi_gouy_all[:, :] - phi_gouy_all[front, :]

    # Calculate dipole phase.
    # Calculate the laser intensity threshold of the target order HHG generated by He^{1+}.
    lambda_cutoff = lambda_0 * 1e6 / q_H
    E_cutoff = 1.2398 / lambda_cutoff
    U_p2 = (E_cutoff - I_p2 / eV) / 3.17
    I_threshold_2 = U_p2 / 9.33 / lambda_0 / 1e6 / lambda_0 / 1e6 * 1e14

    # Calculate the sampling dipole phase.
    E_field_all = np.sqrt(2 * mu_0 * c * ave_I_inter_all)
    E_threshold_2 = np.sqrt(2 * mu_0 * c * I_threshold_2 * 1e4)   
    N_dipole = 2000
    E_dipole = np.linspace(0.9 * np.min(E_field_all), 1.1 * np.max(E_field_all), N_dipole)
    dE = E_dipole[1] - E_dipole[0]
    dipole_s2 = np.zeros(N_dipole)
    dipole_l2 = np.zeros(N_dipole)
    for i in range(N_dipole):
        if E_dipole[i] > E_threshold_2:
            dipole_l2[i], dipole_s2[i] = HHG_DipolePhase(q_H, E_dipole[i], I_p2, omega_d)
        elif E_dipole[i] <= E_threshold_2:
            dipole_l2[i], dipole_s2[i] = HHG_DipolePhase(q_H, E_threshold_2, I_p2, omega_d)

    # Use interpolation to find the dipole phase which is related to laser intensity.
    phi_dipole_l2_all = np.zeros(ave_I_inter_all.shape)
    for j in range(NR):
        for k in range(NZ_total):
            x_grid = math.floor((E_field_all[k, j] - E_dipole[0]) / dE)
            x_rate = (E_field_all[k, j] - E_dipole[x_grid]) / dE
            phi_dipole_l2_all[k, j] = x_rate * (dipole_l2[x_grid+1] - dipole_l2[x_grid]) + dipole_l2[x_grid]
            if ave_I_inter_all[k, j] < I_threshold_2 * 1e4:
                ni_2_all[k, j] = 0
    phi_dipole_l2_all = phi_dipole_l2_all[:, :] - phi_dipole_l2_all[front, :]
    phi_total_2 = phi_gouy_all + phi_plasma_all + phi_dipole_l2_all

    # Reconstruct driving laser.
    r_grid, z_grid = np.meshgrid(r_all, z_all)
    phi_0 = k_0 * r_grid * r_grid / 2 / Curvature_Radius(lambda_0, z_grid) + np.repeat(Gouy_phase_shift(lambda_0, z_all), NR).reshape(NZ_total, NR) + 2 * np.pi / lambda_0 * np.cumsum(Refractive_index_of_plasma(ne_all, omega_d) * DZ, axis = 0)
    Am = np.sqrt(2 * mu_0 * c * ave_I_inter_all)
    A_d = Am * np.exp(phi_0 * 1j)
    
    # Calculate local harmonic field.
    phi_LH_2 = q_H * phi_0 + phi_dipole_l2_all
    source2 = ni_2_all
    E_LH_2 = source2 * Am **5 * np.exp(phi_LH_2 * 1j)
    E_LH_p = source2 * Am **5

    # Calculate propagation phase.
    refractive_index_of_plasma_omega_q = Refractive_index_of_plasma(ne_all, omega_q)
    phi_prop = np.zeros(ne_all.shape)
    for i in range(NZ_total- 1):
        phi_prop[i, :] = (k_0 * r_all * r_all / 2 / Curvature_Radius(lambda_q, z_all[i+1]) - k_0 * r_all * r_all / 2 / Curvature_Radius(lambda_q, z_all[i])) + Gouy_phase_shift(lambda_q, z_all[i+1]) - Gouy_phase_shift(lambda_q, z_all[i]) + 2 * np.pi * q_H / lambda_0 * 0.5 * (refractive_index_of_plasma_omega_q[i, :] + refractive_index_of_plasma_omega_q[i+1, :]) * DZ
    phi_prop[-1, :] = phi_prop[-2, :]

    # Calculate propagation HHG field.
    E_H_2 = np.zeros(ne_all.shape, dtype = complex)
    E_H_p = np.zeros(ne_all.shape, dtype = complex)
    E_H_2[0, :] = E_LH_2[0, :]
    E_H_p[0, :] = E_LH_p[0, :]
    for i in range(1, NZ_total):
        E_H_2[i, :] = E_H_2[i-1, :] * np.exp(phi_prop[i-1, :] * 1j) + E_LH_2[i, :]
        E_H_p[i, :] = E_H_p[i-1, :] + E_LH_p[i, :]
        
    # Calculate HHG yield.
    I_H_2 = 0.5 * epsilon_0 * c * E_H_2 * np.conjugate(E_H_2)
    I_H_p = 0.5 * epsilon_0 * c * E_H_p * np.conjugate(E_H_p)

    # Compare the dipole phase of short trajectory and long trajectory.
    fig1 = plt.figure(figsize=(4, 4), dpi = 500)
    format(fig1)

    ax1 = fig1.add_subplot(1, 1, 1)
    ax_format(ax1, 0.5e16, 0.1e16, 100, 20)
    ax1.plot(E_dipole*E_dipole/2/mu_0/c/1e4, dipole_s2, color = 'royalblue', label = 'short trajectory')
    ax1.plot(E_dipole*E_dipole/2/mu_0/c/1e4, dipole_l2, color = 'seagreen', label = 'long trajectory')
    ax1.set_xlabel(r'$I$ $\mathrm{\left( \frac{W}{cm^2}\right)}$')
    ax1.set_ylabel(r'$\Phi_{\mathrm{dipole}}$ $\mathrm{\left( \frac{rad}{\pi}\right)}$')
    plt.xlim(1.1e16, 2.5e16)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PATH + '/png/dipole.png')
    plt.close()
    
    # Plot the phase-matching conditions.
    fig1 = plt.figure(figsize=(7, 11), dpi = 500)
    format(fig1)

    ax2 = fig1.add_subplot(5, 1, 1)
    ax_format(ax2, 500, 100, 20, 4)
    cm = plt.cm.get_cmap('viridis')
    a = ax2.imshow(phi_gouy_all[front:back, :].transpose() / np.pi, extent = [z_all[front] / um, z_all[back] / um, r_all[0] / um, r_all[r_all.size-1] / um], cmap = cm)
    ax2.axis('tight')
    cbar = plt.colorbar(a)
    cbar.set_label(r'$\Delta \Phi_\mathrm{Gouy}$ $\mathrm{\left(\frac{rad}{\pi}\right)}$')
    ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    ax2.set_ylabel('$r$ $\mathrm{(\mu m)}$')
    if q_H_in == 95:
        ax2.set_xlim(7000, 10000)
    elif q_H_in == 169:
        ax2.set_xlim(8000, 9500)
    ax2.set_ylim(-25, 25)

    ax2 = fig1.add_subplot(5, 1, 2)
    ax_format(ax2, 500, 100, 20, 4)
    cm = plt.cm.get_cmap('viridis')
    a = ax2.imshow(phi_plasma_all[front:back, :].transpose() / np.pi, extent = [z_all[front] / um, z_all[back] / um, r_all[0] / um, r_all[r_all.size-1] / um], cmap = cm)
    ax2.axis('tight')
    cbar = plt.colorbar(a)
    cbar.set_label(r'$\Delta \Phi_\mathrm{plasma}$ $\mathrm{\left(\frac{rad}{\pi}\right)}$')
    ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    ax2.set_ylabel('$r$ $\mathrm{(\mu m)}$')
    if q_H_in == 95:
        ax2.set_xlim(7000, 10000)
    elif q_H_in == 169:
        ax2.set_xlim(8000, 9500)
    ax2.set_ylim(-25, 25)

    ax2 = fig1.add_subplot(5, 1, 3)
    ax_format(ax2, 500, 100, 20, 4)
    cm = plt.cm.get_cmap('viridis')
    a = ax2.imshow(phi_dipole_l2_all[front:back, :].transpose() / np.pi, extent = [z_all[front] / um, z_all[back] / um, r_all[0] / um, r_all[r_all.size-1] / um], cmap = cm)
    ax2.axis('tight')
    cbar = plt.colorbar(a)
    cbar.set_label(r'$\Delta \Phi_\mathrm{dipole}$ $\mathrm{\left(\frac{rad}{\pi}\right)}$')
    ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    ax2.set_ylabel('$r$ $\mathrm{(\mu m)}$')
    if q_H_in == 95:
        ax2.set_xlim(7000, 10000)
    elif q_H_in == 169:
        ax2.set_xlim(8000, 9500)
    ax2.set_ylim(-25, 25)

    ax2 = fig1.add_subplot(5, 1, 4)
    ax_format(ax2, 500, 100, 20, 4)
    cm = plt.cm.get_cmap('viridis')
    a = ax2.imshow(phi_total_2[front:back, :].transpose() / np.pi, extent = [z_all[front] / um, z_all[back] / um, r_all[0] / um, r_all[r_all.size-1] / um], cmap = cm)
    ax2.axis('tight')
    cbar = plt.colorbar(a)
    cbar.set_label(r'$\Delta \Phi_\mathrm{total}$ $\mathrm{\left(\frac{rad}{\pi}\right)}$')
    ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    ax2.set_ylabel('$r$ $\mathrm{(\mu m)}$')
    if q_H_in == 95:
        ax2.set_xlim(7000, 10000)
    elif q_H_in == 169:
        ax2.set_xlim(8000, 9500)
    ax2.set_ylim(-25, 25)

    plt.tight_layout()
    plt.savefig(PATH + '/png/pm2D_'+str(q_H)+'.png')
    plt.close()
    
    # Plot optimized HHG yield and perfect phase-matching condition.
    fig1 = plt.figure(figsize=(12, 8), dpi = 500)
    format(fig1)
    spec = gridspec.GridSpec(ncols=3, nrows=3,width_ratios=[7, 0.5, 3])

    ax2 = fig1.add_subplot(spec[0])
    ax_format(ax2, 500, 100, 20, 4)
    cm = plt.cm.get_cmap('viridis')
    a = ax2.imshow(np.real(I_H_2[front:back, :]).transpose() / np.max(np.real(I_H_p[front:back, int(0.5 * NR)])), extent = [z_all[front] / um, z_all[back] / um, r_all[0] / um, r_all[r_all.size-1] / um], cmap = cm, vmin=0, vmax=1)
    ax2.axis('tight')
    ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    ax2.set_ylabel('$r$ $\mathrm{(\mu m)}$')
    if q_H_in == 95:
        ax2.set_xlim(7000, 10000)
    elif q_H_in == 169:
        ax2.set_xlim(8000, 9500)
    ax2.set_ylim(-25, 25)

    colorAx = fig1.add_subplot(spec[1])
    cbar = plt.colorbar(a, cax=colorAx)
    cbar.set_label(r'$I_\mathrm{H}$ (a.u.)')
   
    y_value = np.real(I_H_2[back, :]) / np.max(np.real(I_H_p[front:back, int(0.5 * NR)]))
    y_peak_half = 0.5 * np.max(y_value)
    for i in range(0, int(0.5 * y_value.size)):
        if y_value[i] < y_peak_half:
            p1 = i + 1
    for i in range(int(0.5 * y_value.size), y_value.size):
        if y_value[i] > y_peak_half:
            p2 = i
    
    ax2 = fig1.add_subplot(spec[2])
    ax_format(ax2, 0.5, 0.1, 20, 4)
    ax2.plot(np.real(I_H_2[back, :]) / np.max(np.real(I_H_p[front:back, int(0.5 * NR)])), r_all / um, color = 'goldenrod')
    ax2.set_xlabel(r'$I_\mathrm{H}$ (a.u.)')
    ax2.set_ylabel('$r$ $\mathrm{(\mu m)}$')
    ax2.set_xlim(-0.05, 1.05)
    ax2.set_ylim(-25, 25)

    ax2 = fig1.add_subplot(spec[3])
    ax_format(ax2, 500, 100, 20, 4)
    cm = plt.cm.get_cmap('viridis')
    a = ax2.imshow(np.real(I_H_p[front:back, :]).transpose() / np.max(np.real(I_H_p[front:back, int(0.5 * NR)])), extent = [z_all[front] / um, z_all[back] / um, r_all[0] / um, r_all[r_all.size-1] / um], cmap = cm, vmin=0, vmax=1)
    ax2.axis('tight')
    ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    ax2.set_ylabel('$r$ $\mathrm{(\mu m)}$')
    if q_H_in == 95:
        ax2.set_xlim(7000, 10000)
    elif q_H_in == 169:
        ax2.set_xlim(8000, 9500)
    ax2.set_ylim(-25, 25)

    colorAx = fig1.add_subplot(spec[4])
    cbar = plt.colorbar(a, cax=colorAx)
    cbar.set_label(r'$I_\mathrm{H, perfect}$ (a.u.)')
   
    ax2 = fig1.add_subplot(spec[5])
    ax_format(ax2, 0.5, 0.1, 20, 4)
    ax2.plot(np.real(I_H_p[back, :]) / np.max(np.real(I_H_p[front:back, int(0.5 * NR)])), r_all / um, color = 'olivedrab')
    ax2.set_xlabel(r'$I_\mathrm{H, perfect}$ (a.u.)')
    ax2.set_ylabel('$r$ $\mathrm{(\mu m)}$')
    ax2.set_xlim(-0.05, 1.05)
    ax2.set_ylim(-25, 25)

    ax2 = fig1.add_subplot(spec[6])
    ax_format(ax2, 500, 100, 0.5, 0.1)
    ax2.plot(z_all[front:back] / um, np.real(I_H_2[front:back, int(0.5 * NR)]) / np.max(np.real(I_H_p[front:back, int(0.5 * NR)])), color = 'goldenrod', label = r'$I_\mathrm{H}$')
    ax2.plot(z_all[front:back] / um, np.real(I_H_p[front:back, int(0.5 * NR)]) / np.max(np.real(I_H_p[front:back, int(0.5 * NR)])), color = 'olivedrab', label = r'$I_\mathrm{H, perfect}$')
    ax2.legend()
    ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    if q_H_in == 95:
        ax2.set_xlim(7000, 10000)
    elif q_H_in == 169:
        ax2.set_xlim(8000, 9500)

    plt.tight_layout()
    plt.savefig(PATH + '/png/result_'+str(q_H)+'.png')
    plt.close()

    # Plot the phase-matching condition on-axis.
    fig1 = plt.figure(figsize=(8, 6), dpi = 500)
    format(fig1)

    ax1 = fig1.add_subplot(2, 1, 1)
    ax_format(ax1, 500, 100, 20, 4)
    ax1.plot(z_all[front:back] / um, phi_plasma_all[front:back, int(0.5 * NR)] / np.pi, label = 'plasma dispersion', color = 'limegreen')
    ax1.plot(z_all[front:back] / um, phi_gouy_all[front:back, int(0.5 * NR)] / np.pi, label = 'Gouy phase', color = 'dodgerblue')
    ax1.plot(z_all[front:back] / um, phi_dipole_l2_all[front:back, int(0.5 * NR)]/np.pi, label='dipole phase',  color = 'orangered')
    ax1.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    ax1.set_ylabel(r'$\Delta \Phi$ $\mathrm{\left(\frac{rad}{\pi}\right)}$')
    ax1.grid()
    ax1.legend()

    ax2 = fig1.add_subplot(2, 1, 2)
    ax_format(ax2, 500, 100, 0.3, 0.1)
    ax2.plot(z_all[front:back] / um, phi_total_2[front:back, int(0.5 * NR)]/np.pi, label='total phase', color = 'goldenrod')
    ax2.set_ylabel(r'$\Delta \Phi_\mathrm{total}$ $\mathrm{\left(\frac{rad}{\pi}\right)}$')
    ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    ax2.grid()
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(PATH + '/png/pm1D_'+ str(q_H) +'.png')
    plt.close()

    bandwidth[ii] = np.real(I_H_2[-1, int(0.5 * NR)])

# Plot frequency spectrum.
fig1 = plt.figure(figsize=(5 ,5), dpi = 500)
format(fig1)

ax1 = fig1.add_subplot(1, 1, 1)
ax_format(ax1, 10, 2, 0.5, 0.1)
ax1.plot(q_H_list, bandwidth / np.max(bandwidth), '--o', label = 'He 2+', color = 'goldenrod')
ax1.set_xlabel('harmonic order')
ax1.set_ylabel(r'$I_\mathrm{H}$ (a.u.)')

plt.tight_layout()
plt.savefig(PATH + '/png/bandwidth_'+str(q_H_in)+'.png')
plt.close()

np.save(PATH + '/npy/q_H_list_'+str(q_H_in), q_H_list)
np.save(PATH + '/npy/bandwidth_'+str(q_H_in), bandwidth)

