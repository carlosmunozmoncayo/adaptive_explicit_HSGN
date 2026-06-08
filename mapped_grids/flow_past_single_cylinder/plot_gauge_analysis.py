import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#Load data skipping forst 8 lines
ref_data_list = []
for gaugeno in range(1, 9):
    filename = f'experimental_data/wg{gaugeno}.txt'
    data = np.loadtxt(filename)
    ref_data_list.append(data)

fig, axs = plt.subplots(
    2, 3, figsize=(18, 8),
    sharex='col',                              # each column shares x (cols have different xlims)
    sharey='row',                              # each row shares y
    gridspec_kw={'wspace': 0, 'hspace': 0},    # zero gap
)

rRlist_top = [4.5, 2.61, 1.66, 1.03]
rRlist_bot = [2.92, 2.29, 1.35, 1.03]
lw = 2
coeff = np.sqrt(9.81)  # To match time

for i in range(3):
    gauge_test = np.loadtxt(f'_output/gauge0000{i+1}.txt', skiprows=2)
    gtime = gauge_test[:, 1]
    geta = gauge_test[:, 2] - 1
    gauge_SWE = np.loadtxt(f'_output_SWE/gauge0000{i+1}.txt', skiprows=2)
    gtime_SWE = gauge_SWE[:, 1]
    geta_SWE = gauge_SWE[:, 2] - 1

    axs[0, i].plot(gtime_SWE*coeff - 7.3, geta_SWE, label='SWE', color='blue', linestyle='dashed', linewidth=lw)
    axs[0, i].plot(gtime*coeff - 7.3, geta, label='HSGN', color='green', linewidth=lw)
    axs[0, i].scatter(ref_data_list[i][:, 0], ref_data_list[i][:, 1],
                      label='Experimental data',# color='red', marker='x')
                      color='k', marker='o', facecolors='none', edgecolors='k')
    axs[0, i].set_xlim(0, 20)

    title_top = f'WG{i+1} ' + r'($\theta=0^\circ$, $r/R=$' + f'{rRlist_top[i]})'
    axs[0, i].text(0.575, 0.1, title_top,
                transform=axs[0, i].transAxes,
                ha='right', va='top', fontsize=15,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.85))
    axs[0,i].set_xlim(0, 20)
    axs[0,i].xaxis.set_major_locator(MaxNLocator(prune='both'))  # avoid x-labels colliding at column border
    axs[0,i].yaxis.set_major_locator(MaxNLocator(prune='both'))  # avoid y-labels colliding at row border

for i in range(3):
    idx = 5 + i
    gauge_test = np.loadtxt(f'_output/gauge0000{idx}.txt', skiprows=2)
    gtime = gauge_test[:, 1]
    geta = gauge_test[:, 2] - 1
    gauge_SWE = np.loadtxt(f'_output_SWE/gauge0000{idx}.txt', skiprows=2)
    gtime_SWE = gauge_SWE[:, 1]
    geta_SWE = gauge_SWE[:, 2] - 1

    axs[1, i].plot(gtime_SWE*coeff - 7.3, geta_SWE, label='SWE', color='blue', linestyle='dashed', linewidth=lw)
    axs[1, i].plot(gtime*coeff - 7.3, geta, label='HSGN', color='green', linewidth=lw)
    axs[1, i].scatter(ref_data_list[idx-1][:, 0], ref_data_list[idx-1][:, 1],
                      label='Experimental data',# color='red', marker='x')
                      color='k', marker='o', facecolors='none', edgecolors='k')
    axs[1, i].set_xlim(0, 20)
    axs[1, i].set_xlabel(r'Time ($t\sqrt{g/h_0}$)', fontsize=15)

    title_bot = f'WG{idx} ' + r'($\theta=150^\circ$, $r/R=$' + f'{rRlist_bot[i]})'
    axs[1, i].text(0.65, 0.1, title_bot,
                        transform=axs[1, i].transAxes,
                        ha='right', va='top', fontsize=15,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.85))
    axs[1,i].set_xlim(0, 20)
    axs[1,i].xaxis.set_major_locator(MaxNLocator(prune='both'))  # avoid x-labels colliding at column border
    axs[1,i].yaxis.set_major_locator(MaxNLocator(prune='both'))  # avoid y-labels colliding at row border

for i in range(3):
    axs[0, i].set_ylim(-0.15, 0.5)
    axs[1, i].set_ylim(-0.06, 0.37)

axs[0, 0].set_ylabel(r'Surface Elevation ($\eta/h_0$)',fontsize=15)
axs[1, 0].set_ylabel(r'Surface Elevation ($\eta/h_0$)',fontsize=15)

# Single shared legend on top
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(labels),
           bbox_to_anchor=(0.5, 1.0), fontsize=15, frameon=True)

fig.subplots_adjust(wspace=0, hspace=0, top=0.94)  # keep gaps at 0, leave room for legend

plt.savefig('Figures/soliton_cylinder_gauge_comparison.pdf', bbox_inches='tight')
plt.show()