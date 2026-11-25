import numpy as np
import matplotlib.pyplot as plt

proc_number = [1,2,4,8,16,32,64,80,96]
proc_number = np.array(proc_number)


mbHSGN = np.array([
    5.419,
    5.219,
    5.889,
    5.680,
    5.641,
    6.590,
    7.522,
    9.033,
    10.674
])

swe = np.array([
    1.608,
    1.422,
    1.610,
    1.610,
    1.618,
    1.995,
    2.363,
    2.871,
    3.259
])


HSGN = np.array([
    7.027,
    7.039,
    7.619,
    7.458,
    7.671,
    8.633,
    9.752,
    11.575,
    14.177,
  ])

SGN= np.array([3.279,
3.164,
3.612,
4.893,
6.841,
11.444,
28.701,
70.283,
95.891,
])

fs = 15
plt.rcParams.update({'font.size': fs})

colors = {
    'HSGN': 'green',
    'mbHSGN': 'tab:orange',
    'SWE': 'blue',
    'SGN': 'magenta'
}
linestyles = {
    'HSGN': '-',
    'mbHSGN': '-',
    'SWE': '-',
    'SGN': '-'
}
markers = {
    'HSGN': 'o',
    'mbHSGN': 's',
    'SWE': 'D',
    'SGN': 'P'
}

parallel_efficiency_SGN    = 100 * SGN[0]   / SGN
parallel_efficiency_HSGN   = 100 * HSGN[0]  / HSGN
parallel_efficiency_mbHSGN = 100 * mbHSGN[0] / mbHSGN
parallel_efficiency_SWE    = 100 * swe[0]   / swe

fig, ax = plt.subplots(figsize=(8, 5))

# ── left axis: parallel efficiency ──────────────────────────────────────────
ax.plot(proc_number[:-3], parallel_efficiency_HSGN[:-3],   label='HSGN',   color=colors['HSGN'],   linestyle=linestyles['HSGN'],   marker=markers['HSGN'])
ax.plot(proc_number[:-3], parallel_efficiency_mbHSGN[:-3], label='mbHSGN', color=colors['mbHSGN'], linestyle=linestyles['mbHSGN'], marker=markers['mbHSGN'])
ax.plot(proc_number[:-3], parallel_efficiency_SWE[:-3],    label='SWE',    color=colors['SWE'],    linestyle=linestyles['SWE'],    marker=markers['SWE'])
ax.plot(proc_number[:-3], parallel_efficiency_SGN[:-3],    label='SGN',    color=colors['SGN'],    linestyle=linestyles['SGN'],    marker=markers['SGN'])

ax.set_xlabel('Number of processes')
ax.set_ylabel('Parallel efficiency (%)')
ax.set_xscale('log', base=2)
ax.axvspan(64, 96, alpha=0.2, color='gray')
handles1, labels1 = ax.get_legend_handles_labels()
ax.legend(handles1 , labels1)
ax.set_ylim(0, 110)
ax.grid(True, axis='y', linestyle='--', alpha=0.7, zorder=-10)
plt.savefig('weak_scaling_fixed_grid_comparison.pdf')
plt.show()