import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.ticker import MaxNLocator

def get_experiment_data(filename):
    data =[]
    with open(filename, 'r') as f:
        for line in f:
            line = line.split(',')[0:2]
            # line = ','.join(line)
            # print(line)
            #Instead of printing, save the line in a list
            array = np.array(line, dtype=float)
            data.append(array)
    data = np.array(data)
    return data

colors = {'SWE':"blue",'SGN':'magenta','HSGN':"#01700B"}
linestyles = {'SWE':'--','SGN':'-.','HSGN':'-'}
lw = 3

fig, axs = plt.subplots(
    2, 3, figsize=(20, 8),
    sharex=True,
    sharey='row',                              # each row shares y → interior y-labels auto-dropped
    gridspec_kw={'wspace': 0, 'hspace': 0},    # zero gap between panels
)
Froude_numbers = ["1.104d0", "1.1702d0"]
alphas = [3.0]
gauge_nos = [5, 4, 3]
offsets = [0, 0, 0]

for i, Froude in enumerate(Froude_numbers):
    for j, gauge_no in enumerate(gauge_nos):
        gauge_data_SGN = np.loadtxt('_output_automated/Froude_{}/SGN/gauge{:05d}.txt'.format(Froude, gauge_no), skiprows=2)
        filename = f'experimental_data/Gauge_{gauge_no}_Fr{Froude[:-2]}.csv'
        exp_data = get_experiment_data(filename)
        ordered_data = exp_data[np.argsort(exp_data[:,0])]
        filtered_data = savgol_filter(ordered_data[:,1], 10, 3)
        ordered_data[:,1] = filtered_data

        

        gauge_data_SWE = np.loadtxt('_output_automated/Froude_{}/SWE/gauge{:05d}.txt'.format(Froude, gauge_no), skiprows=2)
        axs[i, j].plot(gauge_data_SWE[:,1]+offsets[j], gauge_data_SWE[:,2],
                       label='SWE', color=colors['SWE'], linestyle=linestyles['SWE'], lw=lw)

        
        for alpha in alphas:
            gauge_data_EDC = np.loadtxt('_output_automated/Froude_{}/EDC_alpha_{}/gauge{:05d}.txt'.format(Froude, alpha, gauge_no), skiprows=2)
            axs[i, j].plot(gauge_data_EDC[:,1]+offsets[j], gauge_data_EDC[:,2],
                           label='HSGN '+ r"$\alpha=$"+str(int(alpha)),
                           color=colors['HSGN'], linestyle=linestyles['HSGN'], lw=lw)
        axs[i, j].plot(gauge_data_SGN[:,1]+offsets[j], gauge_data_SGN[:,2],
                       label='SGN', color=colors['SGN'], linestyle=linestyles['SGN'], lw=lw)

        # Title inside the axes (upper-right corner)
        if j<2:  
            axs[i, j].text(0.97, 0.95, f'Gauge {gauge_no}, Fr = {Froude[:-2]}',
                        transform=axs[i, j].transAxes,
                        ha='right', va='top', fontsize=15,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.85))
        else:
            axs[i, j].text(0.5, 0.95, f'Gauge {gauge_no}, Fr = {Froude[:-2]}',
                                    transform=axs[i, j].transAxes,
                                    ha='right', va='top', fontsize=15,
                                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.85))

        axs[1, j].set_xlabel('Time (s)', fontsize=15)
        axs[i, j].set_xlim(0, 17.5)
        axs[i, j].xaxis.set_major_locator(MaxNLocator(prune='both'))  # avoid x-labels colliding at column borders
        axs[i, j].yaxis.set_major_locator(MaxNLocator(prune='both'))  # avoid y-labels colliding at row border

        axs[i, j].scatter(exp_data[:,0], exp_data[:,1], color='k', label='Experimental data', s=15,marker='o', facecolors='none', edgecolors='k')

axs[0, 0].set_ylabel('Water height (m)', fontsize=15)
axs[1, 0].set_ylabel('Water height (m)', fontsize=15)

for j in range(3):
    axs[0, j].set_ylim(0.247, 0.315)
    axs[1, j].set_ylim(0.245, 0.368)

# Single shared legend on top
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=len(labels),
           bbox_to_anchor=(0.5, 1.0), fontsize=15, frameon=True)

# plt.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0, top=0.93)  # keep gaps at 0, leave room for legend
plt.savefig('Figures/comparison_gauges_experiment.pdf', bbox_inches='tight')