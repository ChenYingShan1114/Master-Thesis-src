import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
from mpl_toolkits import mplot3d

t0 = 50
w0 = 21.22 * np.sqrt(2) * 1e-6 # 24.75
tau = 30 * 1e-15
lambda0 = 405e-9
k = 2 * np.pi / lambda0
w_FWHM = 1.18 * w0
# xi0 = p.E * 1e-3 / 1.06 / 1.13 / w_FWHM / w_FWHM / tau_FWHM
xi0 = 2.2e+16 #4.66e+16
snap_num = 80 + 1
case_number = 95 #169
path_read = './'+str(case_number)+'_data/'
filename = 'i0.2d'
PATH = './'+str(case_number)+'_png/'

# The proporty of figure (By morris Huang)
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

def ax_format_3D(ax, xmaj, xmin, ymaj, ymin, zmaj, zmin):
    ax.xaxis.set_tick_params(which='major', size=5, width=1, direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=3, width=1, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=5, width=1, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=3, width=1, direction='in', right='on')
    ax.zaxis.set_tick_params(which='major', size=5, width=1, direction='in', right='on')
    ax.zaxis.set_tick_params(which='minor', size=3, width=1, direction='in', right='on')
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(xmaj))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xmin))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ymaj))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ymin))
    ax.zaxis.set_major_locator(mpl.ticker.MultipleLocator(zmaj))
    ax.zaxis.set_minor_locator(mpl.ticker.MultipleLocator(zmin))
    
# Read peak intenisy 'i0', plasma 'ne', gas 'n0', and ions 'n1', 'n2' in the whole simulation region.
# Those data is output from Paul Gibbon's code directly.
# After this part, we get two 1d arrays about z-axis and r-axis coordinate in SI unit,
#                         five 2d arrays about intensity and gases distribution in SI unit,
#                         and the grid number and grid size (NZ, NR, dz, dr).
r_list = []
z_list = []
i0_list = []
ne_list = []
n0_list = []
n1_list = []
n2_list = []
with open(path_read+filename) as f:
    for line in f.readlines():
        s = line.split('  ')
        if (len(s) == 8 and s[0] == ''):
            z_list.append(float(s[1]))
            r_list.append(float(s[2]))
            i0_list.append(float(s[3]))
            ne_list.append(float(s[4]))
            n0_list.append(float(s[5]))
            n1_list.append(float(s[6]))
            n2_list.append(float(s[7]))
        elif (len(s) == 7 and s[0] != ''):
            z_list.append(float(s[0]))
            r_list.append(float(s[1]))
            i0_list.append(float(s[2]))
            ne_list.append(float(s[3]))
            n0_list.append(float(s[4]))
            n1_list.append(float(s[5]))
            n2_list.append(float(s[6]))        
        elif (len(s) == 7 and s[0] == ''):
            s1 = s[1].split(' ')
            z_list.append(float(s1[0]))
            r_list.append(float(s1[1]))
            i0_list.append(float(s[2]))
            ne_list.append(float(s[3]))
            n0_list.append(float(s[4]))
            n1_list.append(float(s[5]))
            n2_list.append(float(s[6]))   
        elif (len(s) == 6):
            s1 = s[0].split(' ')
            z_list.append(float(s1[1]))
            r_list.append(float(s1[2]))
            i0_list.append(float(s[1]))
            ne_list.append(float(s[2]))
            n0_list.append(float(s[3]))
            n1_list.append(float(s[4]))
            n2_list.append(float(s[5]))   
        else:
            print('沒考慮到的情況')
            print(s)
            print(len(s))

z = np.array(z_list)
r = np.array(r_list)
i0 = np.array(i0_list)
ne  = np.array(ne_list)
n0 = np.array(n0_list)
n1 = np.array(n1_list)
n2 = np.array(n2_list)

for i in range(z.size):
    if (z[i] != z[0]):
        nr = i
        break        
nz = int(z.size/nr)

r_1d = np.zeros(nr)
z_1d = np.zeros(nz)
i0_2d = np.zeros([nz, nr])
ne_2d = np.zeros([nz, nr])
n0_2d = np.zeros([nz, nr])
n1_2d = np.zeros([nz, nr])
n2_2d = np.zeros([nz, nr])
for i in range(nz):
    z_1d[i] = z[nr*i] * 1000    # unit: um
for i in range(nr):
    r_1d[i] = r[i]              # unit: um
print('simulation region')
print('nz = ', nz)
print('nr = ', nr)
print('dz = ', z_1d[1]-z_1d[0], ' um')
print('dr = ', r_1d[1]-r_1d[0], ' um')

a = 0
for i in range(nz):
    for j in range(nr):
        i0_2d[i][j] = i0[a] * 1e4   # unit: I/m^2
        ne_2d[i][j] = ne[a] * 1e6   # unit: m^-3
        n0_2d[i][j] = n0[a] * 1e6
        n1_2d[i][j] = n1[a] * 1e6
        n2_2d[i][j] = n2[a] * 1e6
        a = a + 1

# Read beam size w(z) in the whole simulation region.
# After this part, we get two 1d arrays about z and w(z).
filename = 'rmsb0.xy'
x_list = []
y_list = []
with open(path_read+filename) as f:
    for line in f.readlines():
        s = line.split('  ')
        s1 = s[0].split(' ')
        if len(s1) == 2:
            x_list.append(float(s1[0]))
            y_list.append(float(s1[1]))
        if len(s1) == 3:
            x_list.append(float(s1[1]))
            y_list.append(float(s1[2]))
rmsbz = np.array(x_list)
rmsb = np.array(y_list)

# Read and plot the laser intensity and electron density data in the moving window.
energy = np.zeros(snap_num)
ave_z_I_all = np.zeros(snap_num)
ave_I_all = np.zeros([snap_num, nr])
ave_I_all_center = np.zeros([snap_num, nr])
for num in range(snap_num):
    # Read laser intensity
    if num < 10:
        filename = 'rhoe0'+str(num)+'.snap'
    elif num >=10 and num < 100:
        filename = 'rhoe'+str(num)+'.snap'
    rhoe_2d_list = []
    data_list = []
    z_list = []
    r_list = []
    switch = 0
    with open(path_read+filename) as f:
        for line in f.readlines():
            if line=='\n':
                rhoe_2d_list.append(data_list)
                data_list = []
                switch = 1
            else:
                s = line.split('  ')
                data_list.append(float(s[-1]))
                if switch == 0:
                    r_list.append(float(s[-2]))
                if float(s[-2]) == 0:
                    z_list.append(float(s[-6]))
    z_window = np.array(z_list) - t0
    r_window = np.array(r_list)
    rhoe = np.array(rhoe_2d_list)

    # Read electron density
    if num < 10:
        filename = 'pulse0'+str(num)+'.snap'
    elif num >=10 and num < 100:
        filename = 'pulse'+str(num)+'.snap'
    pulse_2d_list = []
    data_list = []
    with open(path_read+filename) as f:
        for line in f.readlines():
            if line=='\n':
                pulse_2d_list.append(data_list)
                data_list = []
            else:
                data_list.append(float(line))
    pulse = np.array(pulse_2d_list)
    
    # Calculate average intensity between FWHM
    pulse_central = pulse[:, int(0.5*(r_window.size))]
    pulse_central_max = np.max(pulse_central)
    nom_pulse_central = pulse_central / pulse_central_max
    FWHM_nom_pulse_central = abs(nom_pulse_central - 0.5)
    sort_FWHM_nom_pulse_central = sorted(FWHM_nom_pulse_central)
    p1 = np.min(np.where(FWHM_nom_pulse_central == sort_FWHM_nom_pulse_central[0]))
    p2 = np.max(np.where(FWHM_nom_pulse_central == sort_FWHM_nom_pulse_central[1]))
    FWHM_nom_pulse_central_si_1 = int(min(np.array([p1, p2])))
    FWHM_nom_pulse_central_si_2 = int(max(np.array([p1, p2])))
    ave_I = np.mean(pulse[FWHM_nom_pulse_central_si_1:FWHM_nom_pulse_central_si_2, :], 0)
    ave_I_all[num, :] = ave_I
    z_I_max = np.argmax(pulse[:, int(0.5*(r_window.size))])
    ave_z_I_all[num] = z_window[z_I_max]
    ave_I_all_center[num] = np.max(pulse[:, int(0.5*(r_window.size))])

    # Integral the laser intensity to get laser energy.
    tmp = 0
    for i in range(z_window.size-1):
        dz = (z_window[i+1] - z_window[i]) * 1e-6
        dt = dz / 299792458
        for j in range(r_window.size-1):
            dr = (r_window[j+1] - r_window[j]) * 1e-6
            tmp = tmp + (0.25 * (pulse[i, j] + pulse[i+1, j] + pulse[i, j+1]+ pulse[i+1, j+1])) * 1e4 * dt * dr * abs(0.5 * 1e-6 * (r_window[j+1] + r_window[j])) * np.pi
    energy[num] = tmp

    # Plot the laser intensity and electron density in the moving window.
    fig1 = plt.figure(figsize=(4, 6), dpi = 500)
    format(fig1)

    pulse_2d = fig1.add_subplot(2, 1, 1)
    ax_format(pulse_2d, 50, 10, 50, 10)
    cm = plt.cm.get_cmap('viridis')
    a = pulse_2d.imshow(pulse.transpose(), extent =[z_window[0], z_window[z_window.size-1], r_window[0], r_window[r_window.size-1]], cmap = cm)
    cbar = plt.colorbar(a)
    cbar.set_label(r'$I_\mathrm{L}$ $\mathrm{\left({\dfrac{W}{cm^2}}\right)}$')
    pulse_2d.axis('tight')
    pulse_2d.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    pulse_2d.set_ylabel('$r$ $\mathrm{(\mu m)}$')

    rhoe_2d = fig1.add_subplot(2, 1, 2)
    ax_format(rhoe_2d, 50, 10, 50, 10)
    cm = plt.cm.get_cmap('plasma')
    a = rhoe_2d.imshow(rhoe.transpose(), extent =[z_window[0], z_window[z_window.size-1], r_window[0], r_window[r_window.size-1]], cmap = cm)
    cbar = plt.colorbar(a)
    cbar.set_label(r'$n_\mathrm{e}$ $\mathrm{\left({cm^{-3}}\right)}$')
    rhoe_2d.axis('tight')
    rhoe_2d.set_xlabel('$z$ $\mathrm{(\mu m)}$')
    rhoe_2d.set_ylabel('$r$ $\mathrm{(\mu m)}$')

    plt.tight_layout()
    if num < 10:
        plt.savefig(PATH + '/png_sim/0'+str(num)+'_2d.png')
    elif num >=10 and num < 100:
        plt.savefig(PATH + '/png_sim/'+str(num)+'_2d.png')
    plt.close()

# z_1d       : the z-axis in the whole simulation region
# ave_z_I_all: the z-axis of laser average intensity
# Align the grid between laser and gas, and perform cubic interpolation of laser intensity.
for i in range(z_1d.size):
    if z_1d[i] < ave_z_I_all[-1]: bd_index = i
z_1d = z_1d[0:bd_index]
i0_2d = i0_2d[0:bd_index, :]
n0_2d = n0_2d[0:bd_index, :]
ne_2d = ne_2d[0:bd_index, :]
n1_2d = n1_2d[0:bd_index, :]
n2_2d = n2_2d[0:bd_index, :]
nz = bd_index

ave_I_inter_all = np.zeros(i0_2d.shape)
for i in range(nr):
    inter = interp1d(ave_z_I_all, ave_I_all[:, i], kind = 'cubic')
    ave_I_inter_all[:, i] = inter(z_1d)
    for k in range(nz):
        if ave_I_inter_all[k, i] < 0:
            ave_I_inter_all[k, i] = 0

# The theoretical solution of 3D Gaussian laser.
rho = r_window[int(0.5*(r_window.size))] * 1e-6
z0 = 0e-6
s = np.pi * w0 * w0 / lambda0
w = w0 * np.sqrt(1 + ((ave_z_I_all * 1e-6 - z0) / s)**2)
R = ((ave_z_I_all * 1e-6 - z0)**2 + s * s) / (ave_z_I_all * 1e-6 - z0)
phi = -np.arctan((ave_z_I_all * 1e-6 - z0) / s)
a0 = np.sqrt( 2 * xi0 / 299792458 / 8.854187817e-12)
psi = a0 * w0 / w * np.exp(-rho**2 / w**2) * np.exp(1j * (k * (ave_z_I_all * 1e-6 + rho**2 / 2 / R) + phi))
intensity = 0.5 * 299792458 * 8.854187817e-12 * psi * np.conjugate(psi)

# Plot the laser intensity, beam size and energy. Compare the result between numerical simulation and theory.
fig1 = plt.figure(figsize=(6, 8), dpi = 500)
format(fig1)

ax1 = fig1.add_subplot(3, 1, 1)
ax_format(ax1, 1000, 200, 0.1e16, 0.02e16)
#ax_format(ax1, 1000, 200, 0.25 * (np.max(np.real(intensity)) - np.min(np.real(intensity))), 0.05 * (np.max(np.real(intensity)) - np.min(np.real(intensity))))
ax1.plot(ave_z_I_all, np.real(intensity), linewidth=3, color='darkslategray', label = 'Theory')
ax1.plot(ave_z_I_all, ave_I_all_center[:, int(0.5*(r_window.size))], 'o', markerfacecolor='none', markersize=5, color='yellowgreen', label = r'$I_\mathrm{peak}$')
#ax1.plot(ave_z_I_all, ave_I_all[:, int(0.5*(r_window.size))], 'o', markersize=1, color='mediumblue', label = r'$\bar{I}_\mathrm{L}$')
ax1.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax1.set_ylabel(r'$I$ $\mathrm{\left(\frac{W}{cm^2}\right)}$')
ax1.grid()
ax1.legend()

ax2 = fig1.add_subplot(3, 1, 2)
ax_format(ax2, 1000, 200, 5, 1)
ax2.plot(ave_z_I_all, w / 1e-6, linewidth=3, color='darkslategray', label = 'Theory')
ax2.plot(rmsbz[::2], rmsb[::2] * np.sqrt(2), 'o', markerfacecolor='none', markersize=5, color='yellowgreen', label = 'Beam size')
ax2.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax2.set_ylabel('$w$ $\mathrm{(\mu m)}$')
ax2.grid()
ax2.legend()

ax3 = fig1.add_subplot(3, 1, 3)
ax_format(ax3, 1000, 200, 0.5 * (np.max(energy) - np.min(energy)) / 1e-3, 0.1 * (np.max(energy) - np.min(energy)) / 1e-3)
ax3.hlines(1e4 * xi0 * w0 * w0 * tau * np.pi**(3/2) / 2 / 1e-3, ave_z_I_all[0], ave_z_I_all[-1], linewidth=3, color='darkslategray', label = 'Theory')
ax3.plot(ave_z_I_all, energy / 1e-3, 'o', markerfacecolor='none', markersize=5, color='yellowgreen', label = 'Energy')
ax3.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax3.set_ylabel('$E$ $\mathrm{(mJ)}$')
ax3.legend()

plt.tight_layout()
plt.savefig(PATH + '/png_sim/check.png')
plt.close()

# Plot the laser average intensity and gases density in the whole simulation region.
fig1 = plt.figure(figsize=(6, 10.5), dpi = 500)
format(fig1)

ax1 = fig1.add_subplot(5,1,1)
ax_format(ax1, 1000, 200, 50, 10)
cm = plt.cm.get_cmap('pink')
a = ax1.imshow(ave_I_inter_all.transpose() * 1e4 / 1e4, extent = [z_1d[0], z_1d[z_1d.size-1], r_1d[0], r_1d[r_1d.size-1]], cmap = cm)
ax1.axis('tight')
cbar = plt.colorbar(a)
cbar.set_label(r'$\bar{I}_\mathrm{L}$ $\mathrm{\left({\dfrac{W}{cm^2}}\right)}$')
ax1.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax1.set_ylabel('$r$ $\mathrm{(\mu m)}$')

ax3 = fig1.add_subplot(5,1,2)
ax_format(ax3, 1000, 200, 50, 10)
cm = plt.cm.get_cmap('viridis')
a = ax3.imshow(n0_2d.transpose() / 1e6, extent = [z_1d[0], z_1d[z_1d.size-1], r_1d[0], r_1d[r_1d.size-1]], cmap = cm)
ax3.axis('tight')
cbar = plt.colorbar(a)
cbar.set_label(r'$n_\mathrm{0}$ $\mathrm{\left({cm^{-3}}\right)}$')
ax3.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax3.set_ylabel('$r$ $\mathrm{(\mu m)}$')

ax7 = fig1.add_subplot(5,1,3)
ax_format(ax7, 1000, 200, 50, 10)
cm = plt.cm.get_cmap('viridis')
a = ax7.imshow(n1_2d.transpose() / 1e6, extent = [z_1d[0], z_1d[z_1d.size-1], r_1d[0], r_1d[r_1d.size-1]], cmap = cm)
ax7.axis('tight')
cbar = plt.colorbar(a)
cbar.set_label(r'$n_\mathrm{1}$ $\mathrm{\left({cm^{-3}}\right)}$')
ax7.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax7.set_ylabel('$r$ $\mathrm{(\mu m)}$')

ax9 = fig1.add_subplot(5,1,4)
ax_format(ax9, 1000, 200, 50, 10)
cm = plt.cm.get_cmap('viridis')
a = ax9.imshow(n2_2d.transpose() / 1e6, extent = [z_1d[0], z_1d[z_1d.size-1], r_1d[0], r_1d[r_1d.size-1]], cmap = cm)
ax9.axis('tight')
cbar = plt.colorbar(a)
cbar.set_label(r'$n_\mathrm{2}$ $\mathrm{\left({cm^{-3}}\right)}$')
ax9.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax9.set_ylabel('$r$ $\mathrm{(\mu m)}$')

ax5 = fig1.add_subplot(5,1,5)
ax_format(ax5, 1000, 200, 50, 10)
cm = plt.cm.get_cmap('viridis')
a = ax5.imshow(ne_2d.transpose() / 1e6, extent = [z_1d[0], z_1d[z_1d.size-1], r_1d[0], r_1d[r_1d.size-1]], cmap = cm)
ax5.axis('tight')
cbar = plt.colorbar(a)
cbar.set_label(r'$n_\mathrm{e}$ $\mathrm{\left({cm^{-3}}\right)}$')
ax5.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax5.set_ylabel('$r$ $\mathrm{(\mu m)}$')

plt.tight_layout()
plt.savefig(PATH + '/png/intensity.png')
plt.close()

# Plot the laser average intensity of each moving window.
fig = plt.figure(figsize = [10, 10], dpi = 500)
format(fig)

ax = plt.axes(projection ='3d')
ax_format_3D(ax, 1000, 500, 100, 20, 0.5e16, 0.1e16)

for i in range(ave_z_I_all.size):
    z_same = np.zeros(r_1d.size)
    z_same[:] = ave_z_I_all[i]
    ax.plot(z_same, r_1d, ave_I_all[i,:], linewidth = 1, color = 'brown')

ax.set_xlabel('$z$ $\mathrm{(\mu m)}$')
ax.set_ylabel('$r$ $\mathrm{(\mu m)}$')
ax.set_zlabel(r'$\bar{I}_\mathrm{L}$ $\mathrm{\left({\dfrac{W}{cm^2}}\right)}$')

plt.savefig(PATH + '/png/ave_intensity.png')
plt.close()

# Save data
np.save(PATH + '/npy/z_1d', z_1d)
np.save(PATH + '/npy/r_1d', r_1d)
np.save(PATH + '/npy/i0_2d', ave_I_inter_all * 1e4)
np.save(PATH + '/npy/ne_2d', ne_2d) # unit: m
np.save(PATH + '/npy/n0_2d', n0_2d) # unit: m
np.save(PATH + '/npy/n1_2d', n1_2d) # unit: m
np.save(PATH + '/npy/n2_2d', n2_2d) # unit: m
