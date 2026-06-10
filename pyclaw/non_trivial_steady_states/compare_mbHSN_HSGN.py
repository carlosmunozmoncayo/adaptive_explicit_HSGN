import numpy as np
from clawpack.pyclaw import Solution
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

frame = 0
exact_sol = Solution()
exact_sol.read(frame=0,path='_output_HSGN')

x = exact_sol.grid.x.centers
exact_h = exact_sol.q[0,:]
bath = exact_sol.aux[0,:]
exact_eta = exact_h + bath

# Read and compare at 3 different times and plot in subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
fig.subplots_adjust(wspace=0, top=0.88)  # remove gap between axes, leave room for legend

frames = [50, 100, 200]
times = [5, 10, 20]
fs1 =15
fs2 =13
plt.rcParams.update({'font.size': fs1})
for i, frame in enumerate(frames):
    solmbHSGN = Solution()
    solHSGN = Solution()
    solmbHSGN.read(frame=frame, path='_output_mbHSGN')
    solHSGN.read(frame=frame, path='_output_HSGN')

    hmbHSGN = solmbHSGN.q[0, :]
    hHSGN = solHSGN.q[0, :]

    etambHSGN = hmbHSGN + bath
    etaHSGN = hHSGN + bath

    axs[i].plot(x, etambHSGN, label='mbHSGN', color="#FF6600", lw=3)
    axs[i].plot(x, etaHSGN, label='HSGN', color='green', lw=3)
    axs[i].plot(x, exact_eta, label='Exact SGN', color='black', linestyle=':', lw=3)
    axs[i].set_xlabel(r'$x \, \left(m\right)$', fontsize=fs2)
    axs[i].set_xlim(-7.5, 7.5)

    # textbox with the time instead of a title
    axs[i].text(0.05, 0.95, f't = {times[i]} s', transform=axs[i].transAxes,
                va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8))
    axs[i].xaxis.set_major_locator(MaxNLocator(prune='both'))

axs[0].set_ylabel(r'$\eta\, \left(m\right)$', fontsize=fs2)


handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.0))
plt.savefig('mbHSGN_vs_HSGN_surface_elevation.pdf')
plt.show()