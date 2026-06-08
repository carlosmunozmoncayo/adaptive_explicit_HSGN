import numpy as np
import matplotlib.pyplot as plt
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
from clawpack.pyclaw import Solution
#Import gaussian filter for smoothing
from scipy.ndimage import gaussian_filter1d

mapc2p = make_mapc2p('celledges.data')

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

colors = {'mbHSGN_alpha3_WP': 'green', 'mbHSGN_alpha3_noWP': 'black', 'SGN': 'magenta'}


fig, ax = plt.subplots(1,2, figsize=(15,4), gridspec_kw={'width_ratios': [1, 3]})
#The inital condition
frame = 0
sol_WP = Solution()
sol_WP.read(frame, path='_output_mbHSGN_alpha3_WP', file_format='ascii', file_prefix='fort')
xc_WP = mapc2p[0](sol_WP.grid.x.centers)
sol_noWP = Solution()
sol_noWP.read(frame, path='_output_mbHSGN_alpha3_noWP', file_format='ascii', file_prefix='fort')
xc_noWP = mapc2p[0](sol_noWP.grid.x.centers)
sol_SGN = Solution()
sol_SGN.read(frame, path='_output_SGN', file_format='ascii', file_prefix='fort')
xc_SGN = mapc2p[0](sol_SGN.grid.x.centers)

eta_WP = sol_WP.q[-1,:]
ax[0].plot(xc_WP, eta_WP, label='HSGN WP', color=colors['mbHSGN_alpha3_WP'])
eta_noWP = sol_noWP.q[-1,:]
ax[0].plot(xc_noWP, eta_noWP, linestyle='--', label='HSGN IP', color=colors['mbHSGN_alpha3_noWP'])
eta_SGN = sol_SGN.q[-1,:]
ax[0].plot(xc_SGN, eta_SGN, label='SGN', color=colors['SGN'], linestyle=':')
ax[0].set_xlim(-600, 600)
ax[0].set_ylabel(r'$\eta\, (m)$')
ax[0].set_xlabel(r'$x\, (m)$')
ax[0].set_title(f"t={sol_WP.t:.0f} s")
ax[0].legend()
# ax[0].legend(loc='lower right')
#The effect
frame = 25
sol_WP = Solution()
sol_WP.read(frame, path='_output_mbHSGN_alpha3_WP', file_format='ascii', file_prefix='fort')
xc_WP = mapc2p[0](sol_WP.grid.x.centers)
sol_noWP = Solution()
sol_noWP.read(frame, path='_output_mbHSGN_alpha3_noWP', file_format='ascii', file_prefix='fort')
xc_noWP = mapc2p[0](sol_noWP.grid.x.centers)
sol_SGN = Solution()
sol_SGN.read(frame, path='_output_SGN', file_format='ascii', file_prefix='fort')
xc_SGN = mapc2p[0](sol_SGN.grid.x.centers)

eta_WP = sol_WP.q[-1,:]
ax[1].plot(xc_WP, eta_WP, label='HSGN WP', color=colors['mbHSGN_alpha3_WP'])
eta_noWP = sol_noWP.q[-1,:]
ax[1].plot(xc_noWP, eta_noWP, linestyle='--', label='HSGN IP ', color=colors['mbHSGN_alpha3_noWP'])
eta_SGN = sol_SGN.q[-1,:]
ax[1].plot(xc_SGN, eta_SGN, label='SGN', color=colors['SGN'], linestyle=':')
ax[1].set_xlim(0, 7500)
ax[1].legend(loc='lower right')
ax[1].set_title(f"t={sol_WP.t:.0f} s")
# ax[1].set_ylabel(r'$\eta\, (m)$')
ax[1].set_xlabel(r'$x\, (m)$')

# --- Inset axis ---
x1, x2 = 3000, 7500
y1, y2 = -0.015, 0.015

# axins = inset_axes(ax, width="50%", height="40%", loc='upper right')
# Create a new axes positioned outside the main plot (in figure fraction coords)
# [left, bottom, width, height] in figure fractions
axins = fig.add_axes([0.62, 0.55, 0.325, 0.325]) 

axins.plot(xc_WP, eta_WP, color=colors['mbHSGN_alpha3_WP'])
axins.plot(xc_noWP, eta_noWP, linestyle='--', color=colors['mbHSGN_alpha3_noWP'])
axins.plot(xc_SGN, eta_SGN, color=colors['SGN'], linestyle=':')

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.tick_params(labelsize=7)

# Draw box on main axis and connecting lines to inset
mark_inset(ax[1], axins, loc1=3, loc2=4, fc="none", ec="black", lw=0.8, linestyle='-')
plt.tight_layout()
plt.show()
# plt.savefig('Figures/WP_surface_elevation_1D.pdf')

xlims = [(-600,600),(0,7500),(0,7500)]
fig, ax = plt.subplots(1,3, figsize=(15,3))


for i, frame in enumerate([0, 10, 25]):
    sol_WP = Solution()
    sol_WP.read(frame, path='_output_mbHSGN_alpha3_WP', file_format='ascii', file_prefix='fort')
    xc_WP = mapc2p[0](sol_WP.grid.x.centers)
    sol_noWP = Solution()
    sol_noWP.read(frame, path='_output_mbHSGN_alpha3_noWP', file_format='ascii', file_prefix='fort')
    xc_noWP = mapc2p[0](sol_noWP.grid.x.centers)
    sol_SGN = Solution()
    sol_SGN.read(frame, path='_output_SGN', file_format='ascii', file_prefix='fort')
    xc_SGN = mapc2p[0](sol_SGN.grid.x.centers)

    #The initial data
    eta_WP = sol_WP.q[-1,:] #This is the water depth for all
    h_WP = sol_WP.q[0,:] #This is the water depth for all --- IGNORE ---
    #Only p changes between the well-prepared and not well-prepared (zero here) cases, so we can only take the well-prepared case for the initial data
    p_WP = sol_WP.q[3,:]/h_WP

    h_noWP = sol_noWP.q[0,:]
    p_noWP = sol_noWP.q[3,:]/h_noWP
    eta_noWP = sol_noWP.q[-1,:]

    # p_SGN = get_p_SGN('_output_SGN', frame)

    # ax[i].plot(xc_SGN, p_SGN, label=r'$p\, (m^2/s^2)$ SGN', color = 'magenta', linestyle=':')
    ax[i].plot(xc_WP, p_WP, label=r'HSGN WP', color = colors['mbHSGN_alpha3_WP'], linestyle='-')
    ax[i].plot(xc_noWP, p_noWP, linestyle='--', label=r'HSGN IP', color = 'black')
    ax[i].legend()
    ax[i].set_xlim(xlims[i])
    if i>0:
        ax[i].set_ylim(-5, 5)
    ax[i].set_title(f"t={sol_WP.t:.0f} s")
    
    ax[i].set_xlabel(r'$x\, (m)$')
ax[0].set_ylabel(r'$p\, (m^2/s^2)$')
plt.tight_layout()
plt.show()
# plt.savefig('Figures/WP_nonhydrostatic_pressure_1D.pdf')