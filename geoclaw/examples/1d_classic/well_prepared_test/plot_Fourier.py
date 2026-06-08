import numpy as np
import matplotlib.pyplot as plt
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
from clawpack.pyclaw import Solution
#Import gaussian filter for smoothing
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

mapc2p = make_mapc2p('celledges.data')


frame = 25

sol_HSGN = Solution()
sol_HSGN.read(frame, path='_output_mbHSGN_alpha3_WP', file_format='ascii', file_prefix='fort')
eta_HSGN = sol_HSGN.q[-1,:]
mapc2p = make_mapc2p('celledges.data')
xc_HSGN = mapc2p[0](sol_HSGN.grid.x.centers)

sol_SGN = Solution()
sol_SGN.read(frame, path='_output_SGN', file_format='ascii', file_prefix='fort')
eta_SGN = sol_SGN.q[-1,:]
xc_SGN = mapc2p[0](sol_SGN.grid.x.centers)


#Fourier analysis of the solution at t=25s
fft_HSGN = np.fft.fft(eta_HSGN)
fft_SGN = np.fft.fft(eta_SGN)
freq = np.fft.fftfreq(len(xc_HSGN), d=xc_HSGN[1] - xc_HSGN[0])
sorted_indices = np.argsort(freq)
scaled_freq = freq*100


#Median filter for trendline
from scipy.signal import medfilt

colors = {'mbHSGN_alpha3_WP': 'green', 'mbHSGN_alpha3_noWP': 'black', 'SGN': 'magenta'}


errors = np.abs(fft_HSGN - fft_SGN)/np.abs(fft_SGN)
median_filtered_errors = medfilt(errors, kernel_size=89)

interp_data = (median_filtered_errors[sorted_indices])[np.where((scaled_freq[sorted_indices] >= 0) & (scaled_freq[sorted_indices] < 1005))]
interp_freq = (scaled_freq[sorted_indices])[np.where((scaled_freq[sorted_indices] >= 0) & (scaled_freq[sorted_indices] < 1005))]

def stretched_exp(x, A, k, beta):
    return A * (1 - np.exp(-(k * x)**beta))

popt, _ = curve_fit(stretched_exp, interp_freq, interp_data, p0=[np.max(interp_data),1,1], maxfev=100000)
fitted_curve = stretched_exp(interp_freq, *popt)

fig, ax = plt.subplots(1,2, figsize=(15,3))

ax[0].semilogy(scaled_freq[sorted_indices], np.abs(fft_HSGN)[sorted_indices], label='HSGN WP', color=colors['mbHSGN_alpha3_WP'])
ax[0].semilogy(scaled_freq[sorted_indices], np.abs(fft_SGN)[sorted_indices], label='SGN', color=colors['SGN'], linestyle=':')
ax[0].set_xlim(0, 1.5)
ax[0].set_ylim(1.e-5,1.e4)
ax[0].set_xlabel(r'$kH$')
ax[0].set_title('Fourier spectrum of surface elevation '+r'$|\hat{\eta}(kH)|$')
ax[0].axvspan(1, 2, color='lightgray', alpha=0.5)
ax[0].text(0.35, 0.1, 'Shallow water regime', transform=ax[0].transAxes, fontsize=12, ha='center', va='center', 
           bbox=dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=1))
ax[0].legend()


#Light grey background at the right of x=1
ax[1].axvspan(1, 2, color='lightgray', alpha=0.5)
ax[1].text(0.35, 0.1, 'Shallow water regime', transform=ax[1].transAxes, fontsize=12, ha='center', va='center', 
           bbox=dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=1))


ax[1].semilogy(scaled_freq[sorted_indices], errors[sorted_indices],label=r'$\frac{|\hat{\eta}_{HSGN}(kH) - \hat{\eta}_{SGN}(kH)|}{|\hat{\eta}_{SGN}(kH)|}$')
ax[1].semilogy(interp_freq, fitted_curve, color='red', label=r'$A\left[1-\exp(-\gamma\, kH)^{\beta}\right]$', linestyle='-.')
ax[1].set_xlabel(r'$kH$')
ax[1].set_title('Relative Fourier error')
ax[1].set_xlim(0, 1.5)
ax[1].set_ylim(1e-4, 1e3)
ax[1].legend(loc='lower right')
plt.tight_layout()
plt.show()
plt.savefig('Figures/Fourier_analysis.pdf')