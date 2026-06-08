import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from clawpack import pyclaw
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
import clawpack.pyclaw.gauges as gauges
from matplotlib.ticker import MaxNLocator


get_n_frames = lambda path: len([name for name in os.listdir(path) if (os.path.isfile(os.path.join(path, name))
                                                                       and name.startswith('fort.q'))])

color_HSGN = 'green'
color_mbHSGN = "#FF6600"
color_SGN = 'magenta'
color_SWE = 'blue'
color_experiment = 'black'

pathSGN = '_output_SGN'
pathSWE = '_output_SWE'
pathHSGN = '_output_HSGN'

n_frames = get_n_frames(pathSGN)

xp_gauges = [-150, -80., -60, -50, -40, -32, -31.2, -30.8, -30,
             -28, -20, -10, -5]

gauge_info = [(-80,   1,(10,60)),
              (-60,   3,(20,70)),
              (-50,   5,(30,80)),
              (-40,   7,(40,90)),
              (-32,   9,(40,90)),
              (-31.2,11,(40,90)),
              (-30.8, 2,(40,90)),
              (-30,   4,(40,90)),
              (-28,   6,(40,90)),
              (-20,   8,(50,100)),
              (-10,  10,(60,110)),
              (-5,   12,(60,110))]

fig, ax = plt.subplots(
    2, 2, figsize=(16, 7),
    sharex='col',                              # each column shares x (cols have different xlims)
    sharey='row',                              # each row shares y
    gridspec_kw={'wspace': 0, 'hspace': 0},    # zero gap
)
axflat = ax.flatten()

gauge_info_subset = [gauge_info[2], gauge_info[3], gauge_info[6], gauge_info[-3]]
add_data = True
experimental_data = np.loadtxt('data_wavegauge.csv', skiprows=2, delimiter=',')

for i, (xp_g, _, (xlim1, xlim2)) in enumerate(gauge_info_subset):
    gaugeno = int(-xp_g*10)
    print(f"Gauge {gaugeno} at x = {xp_g} m")
    gaugeSGN = gauges.GaugeSolution(gaugeno, pathSGN)
    gaugembSGN = gauges.GaugeSolution(gaugeno, pathHSGN)
    gaugeSWE = gauges.GaugeSolution(gaugeno, pathSWE)
    tSGN = gaugeSGN.t
    etaSGN = gaugeSGN.q[-1, :]
    tHSGN = gaugembSGN.t
    etaHSGN = gaugembSGN.q[-1, :]
    tSWE = gaugeSWE.t
    etaSWE = gaugeSWE.q[-1, :]

    axflat[i].plot(tSGN, 100*etaSGN, color=color_SGN, label="SGN", lw=2)
    axflat[i].plot(tHSGN, 100*etaHSGN, color=color_HSGN, label="HSGN", lw=2)
    axflat[i].plot(tSWE, 100*etaSWE, color=color_SWE, linestyle='--', label="SWE", lw=2)
    if add_data:
        axflat[i].plot(experimental_data[:, 0], experimental_data[:, i+2],
                       color=color_experiment, linestyle=':', label="Experimental data", lw=2)

    if i==0 or i==2:
        xlim1, xlim2 = 30, 90
    else:
        xlim1, xlim2 = 40, 100
    if i<2:
        ylim1, ylim2 = -6, 11
    else:
        ylim1, ylim2 = -5, 14
        
    axflat[i].set_xlim(xlim1, xlim2)
    axflat[i].set_ylim(ylim1, ylim2)
    axflat[i].xaxis.set_major_locator(MaxNLocator(prune='both'))  # avoid x-labels colliding at column border
    axflat[i].yaxis.set_major_locator(MaxNLocator(prune='both'))  # avoid y-labels colliding at row border

    # Title inside the axes (upper-right corner)
    axflat[i].text(0.97, 0.95, f"Gauge at x = {xp_g} m",
                   transform=axflat[i].transAxes,
                   ha='right', va='top', fontsize=15,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.85))

axflat[0].set_ylabel('Surface elevation (cm)', fontsize=15)
axflat[2].set_ylabel('Surface elevation (cm)', fontsize=15)
axflat[2].set_xlabel('Time (s)', fontsize=15)
axflat[3].set_xlabel('Time (s)', fontsize=15)

# Single shared legend on top
handles, labels = axflat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(labels),
           bbox_to_anchor=(0.5, 1.0), fontsize=15, frameon=True)

fig.subplots_adjust(wspace=0, hspace=0, top=0.92)  # keep gaps at 0, leave room for legend

fname = 'Matsuyama_gauge_comparison.pdf'
plt.savefig(fname, bbox_inches='tight')
print('Created %s' % fname)