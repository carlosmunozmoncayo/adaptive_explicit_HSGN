import numpy as np
import matplotlib.pyplot as plt

alphas = [3.0, 4.242640687119286, 6.000000000000002, 8.485281374238571]
mxs = [25, 50, 100, 200, 400]
mxs = np.array(mxs)

times_SWE =  np.array([ 0.448 ,  0.868 ,  2.128  , 9.102 ,  52.446 ])
#Row 1: mx=25, Row 2: mx=50, Row 3: mx=100, Row 4: mx=200, Row 5: mx=400
times_mbHSGN = np.array([ 
  [ 1.458  , 1.971      ,2.730  , 3.755 ],
  [ 4.631  , 6.221      ,9.146   ,12.177 ],
  [ 19.080 ,  25.167    ,34.863  , 48.770 ],
  [ 112.730,   154.711  , 214.596,   308.301 ],
  [ 809.529,   1115.018 ,  1566.265,   2197.213 ]
])  

times_SGN = np.array([1.248,   5.148,   27.000,   361.834,  5727.394])

color_SGN = 'magenta'
marker_SGN = 'P'
#Slowly changing colors for mbHSGN lines
colors_mbHSGN = plt.cm.turbo(np.linspace(0.2, 0.9, len(alphas)))
colors_mbHSGN = ["#8c4a00", "#cc7000", "#ff9500", "#ffb84d"][::-1]
markers_mbHSGN = ['o', 's', 'D', '^']
labels_power = [r'$\alpha=3$', r'$\alpha=3\cdot 2^{0.5}$', r'$\alpha=3\cdot 2^{1}$', r'$\alpha=3\cdot 2^{1.5}$']

fs = 14
plt.rcParams.update({'font.size': fs})

fig, ax = plt.subplots(figsize=(8, 5))
for i, alpha in enumerate(alphas):
    ax.plot(mxs, times_mbHSGN[:,i], marker=markers_mbHSGN[i], label=f'mbHSGN ('+r'$\alpha=$'+f'{alpha:.2f})', color=colors_mbHSGN[i])

ax.plot(mxs, times_SGN, marker=marker_SGN,
         label='SGN', color=color_SGN, linestyle='-.')

#Plot also SWE
ax.plot(mxs, np.array(times_SWE), marker='x', label='SWE', color='blue', linestyle='--')

#Plot reference lines for power-law scalings
quadratic = 1.e-5*mxs**2.5
ax.plot(mxs[2:], quadratic[2:], label=r'$\mathcal{O}(N^{2.5})$', linestyle='--', color='k')

cubic = 1.e-5*mxs**3
ax.plot(mxs[2:], cubic[2:], label=r'$\mathcal{O}(N^{3})$', linestyle='-', color='k')

quartic = 4.e-7*mxs**4
ax.plot(mxs[2:], quartic[2:], label=r'$\mathcal{O}(N^{4})$', linestyle=':', color='k')


ax.set_xlabel(r'$N \text{ (number of cells per dimension)}$')
ax.set_ylabel(r'Wall-clock time')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xticks(mxs)
ax.set_xticklabels(mxs)
ax.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(ncols=2)
plt.savefig('mbHSGN_SWE_ratio_plot.pdf', bbox_inches='tight')
plt.show()