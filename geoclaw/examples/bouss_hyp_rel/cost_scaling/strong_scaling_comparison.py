import numpy as np
import matplotlib.pyplot as plt

proc_number = [1,2,4,8,16,32,64,80,96]
proc_number = np.array(proc_number)


HSGN = np.array([
449.768,
234.634,
121.439,
60.111,
 31.092,
 17.594,
 9.522,
 9.430,
 9.432])

mbHSGN = np.array([
  340.748,
  173.303,
  87.857,
  44.727,
  23.218,
  12.354,
  7.302,
  7.372,
  7.486
])

swe = np.array([
  92.763,
  48.475,
  24.081,
  12.218,
  6.661,
  3.626,
  2.230,
  2.289,
  2.345
])

SGN = np.array([
508.408,
262.973,
159.813,
131.322,
 63.898,
 38.987,
 34.393,
 63.968,
 84.130,
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