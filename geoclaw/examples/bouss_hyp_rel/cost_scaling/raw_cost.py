import matplotlib.pyplot as plt
import numpy as np

labels = ['Radial flat (s)', 'Tohoku tsunami (min)']
solvers = ['SGN', 'HSGN', 'mbHSGN', 'SWEs']

fs = 14
plt.rcParams.update({'font.size': fs})

# SGN is split into two components (half each for this example)
sgn_comp1 = [42.1-34.214, 159.2-123.5]
sgn_comp2 = [34.214, 123.5]

data = {
    'HSGN':   [45.6, 118.3],
    'mbHSGN': [30.8,  73.3],
    'SWEs':   [ 2.8,  19.5],
}

colors = {
    'SGN': 'tab:magenta',
    'HSGN': 'green',
    'mbHSGN': 'tab:orange',
    'SWEs': 'blue'
}

x = np.arange(len(labels))
n = len(solvers)  # still 4 bars per group
width = 0.2

fig, ax = plt.subplots(figsize=(6, 4))

# Plot SGN as a stacked bar (position i=0)
ax.grid(True, axis='y', linestyle='--', alpha=0.7, zorder=-10)
offset = (0 - n / 2 + 0.5) * width
ax.bar(x + offset, sgn_comp1, width, label='SGN (explicit part)', color='magenta')
ax.bar(x + offset, sgn_comp2, width, bottom=sgn_comp1, label='SGN (linear solve)', color='purple')

# Plot the rest
for i, solver in enumerate(list(data.keys()), start=1):
    offset = (i - n / 2 + 0.5) * width
    ax.bar(x + offset, data[solver], width, label=solver, color=colors[solver])

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('Wall clock time', fontsize=fs)

ax.legend()
plt.tight_layout()
plt.savefig('timing_comparison.pdf')
plt.show()