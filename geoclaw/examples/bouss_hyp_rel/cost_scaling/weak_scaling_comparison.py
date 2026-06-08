import numpy as np
import matplotlib.pyplot as plt

proc_number = [1,2,4,8,16,32,64,80,96]
proc_number = np.array(proc_number)


mbHSGN = np.array([
    6.691,
    6.674,
    7.025,
    6.794,
    7.115,
    7.330,
    9.670,
    18.313,
    17.499
])

swe = np.array([
    1.521,
    1.413,
    1.445,
    1.501,
    1.571,
    1.727,
    2.252,
    3.650,
    4.025
])


HSGN = np.array([
    9.762,
10.401,
10.508,
10.357,
 10.757,
 11.239,
 14.031,
 25.091,
 27.662,
  ])

SGN= np.array([3.129,
3.193,
3.431,
4.464,
6.822,
10.959,
19.677,
59.926,
95.329,
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
# plt.savefig('weak_scaling_fixed_grid_comparison.pdf')
plt.show()