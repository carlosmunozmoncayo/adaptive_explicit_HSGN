import numpy as np
import matplotlib.pyplot as plt

proc_number = [1,2,4,8,16,32,64,80,96]
proc_number = np.array(proc_number)


HSGN = np.array([
647.056,
327.662,
162.885,
81.367,
 42.725,
 22.426,
 12.834,
 13.060,
 12.384])

mbHSGN = np.array([
  437.545,
  219.638,
  109.067,
  56.814,
  28.917,
  15.249,
  9.354,
  8.924,
  8.866
])

swe = np.array([
  92.626,
  45.495,
  23.351,
  12.019,
  6.483,
  3.714,
  2.376,
  2.403,
  2.425
])

SGN = np.array([
528.048,
282.934,
169.412,
140.273,
 67.137,
 40.577,
 31.259,
 62.744,
 88.235,
])  

fs = 15
plt.rcParams.update({'font.size': fs})

colors ={
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

parallel_efficiency_HSGN = 100*proc_number[0]/(HSGN/HSGN[0])/proc_number
parallel_efficiency_mbHSGN = 100*proc_number[0]/(mbHSGN/mbHSGN[0])/proc_number
parallel_efficiency_SWE = 100*proc_number[0]/(swe/swe[0])/proc_number
parallel_efficiency_SGN = 100*proc_number[0]/(SGN/SGN[0])/proc_number


fig, ax = plt.subplots(figsize=(8,5))
ax.plot(proc_number[:-3], parallel_efficiency_HSGN[:-3], label='HSGN', color=colors['HSGN'], linestyle=linestyles['HSGN'], marker=markers['HSGN'])
ax.plot(proc_number[:-3], parallel_efficiency_mbHSGN[:-3], label='mbHSGN', color=colors['mbHSGN'], linestyle=linestyles['mbHSGN'], marker=markers['mbHSGN'])
ax.plot(proc_number[:-3], parallel_efficiency_SWE[:-3], label='SWE', color=colors['SWE'], linestyle=linestyles['SWE'], marker=markers['SWE'])
ax.plot(proc_number[:-3], parallel_efficiency_SGN[:-3], label='SGN', marker=markers['SGN'], color=colors['SGN'], linestyle=linestyles['SGN'])

ax.set_xlabel('Number of processes')
ax.set_ylabel('Parallel efficiency(%)')
#G background to the right of 64
ax.axvspan(64, 96, alpha=0.2, color='gray')
ax.set_xscale('log', base=2)
ax.grid(True, axis='y', linestyle='--', alpha=0.7, zorder=-10)
ax.set_ylim(0, 110)
ax.legend()
plt.savefig('strong_scaling_fixed_grid_comparison.pdf')