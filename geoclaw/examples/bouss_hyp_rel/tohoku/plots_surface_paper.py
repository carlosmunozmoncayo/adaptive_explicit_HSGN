import numpy as np
import matplotlib.pyplot as plt
from clawpack.pyclaw import Solution
from clawpack.visclaw import gridtools 
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
from matplotlib.ticker import FuncFormatter   # put up with your other imports


plt.rcParams["patch.antialiased"] = False


# -------------------------------------------------
# Colormap + normalization
# -------------------------------------------------
linthresh = 0.2  # 2 cm sensitivity 

norm = mcolors.SymLogNorm(
    linthresh=linthresh,
    vmin=-2,
    vmax=2
)


topo_levels = [-np.inf, 50, 300, np.inf]
topo_colors = ["#D2C3A3", "#9DAA7A", "#7A8F5F"]
# topo_colors = ["red","red","red"]

path_mbHSGN = '_output_UCSB_mbHSGN'
path_HSGN = '_output_UCSB_HSGN'
path_SWE = '_output_UCSB_SWE'
path_SGN = '_output_UCSB_SGN'

def read_all_patches(frame_nos, remove_overlaps=False,
                     format='ascii'):
    #Will return a list of lists of patches, where each patch is a dictionary with keys "X", "Y", "Z" and "level".
    # The outer list is over frames, and the inner list is over patches in that frame. 
    all_patches = []
    for frame_no in frame_nos:
        patches = []
        sol = Solution()
        sol.read(frame=frame_no, path=path,
                 read_aux=False,format=format)
        for state in sol.states:
            patch = {}
            xn = state.grid.x.nodes
            yn = state.grid.y.nodes
            XN, YN = np.meshgrid(xn, yn)
            patch['X'], patch['Y'] = XN.T, YN.T
            # patch['ETA'] = state.q[variable_index,:,:]+offset_z 
            xc = state.grid.x.centers
            yc = state.grid.y.centers
            XC, YC = np.meshgrid(xc, yc)
            patch['XC'], patch['YC'] = XC.T, YC.T
            # patch['ETA'] = state.q[variable_index,:,:]+offset_z 
            eta = state.q[-1,:,:]
            h = state.q[0,:,:]
            eta = np.where(h>0, eta, np.nan)
            bath = state.q[-1,:,:]-state.q[0,:,:] #Assuming eta = h + bath
            bath = np.where(h<=0, bath, np.nan) # 
            patch['ETA'] = eta
            patch['BATH'] = bath
            patch['level'] = state.patch.level
            patches.append(patch)
        if remove_overlaps:
            for p in patches:
                level = p["level"]
                #Go over patches in higher levels and set to nan the points of p that are covered by those patches.
                for p_higher in patches:
                    if p_higher["level"] > level:
                        x_min, x_max = p_higher["XC"].min(), p_higher["XC"].max()
                        y_min, y_max = p_higher["YC"].min(), p_higher["YC"].max()

                        mask = (
                            (p["XC"] >= x_min) &
                            (p["XC"] <= x_max) &
                            (p["YC"] >= y_min) &
                            (p["YC"] <= y_max)
                        )
                        p["ETA"][mask] = np.nan
                        p["BATH"][mask] = np.nan
        all_patches.append(patches)
    return all_patches

def fill_one_axis_surface_2D(axis,patches,plot_boundary_levels=[],
                             plot_cell_levels=[],cmap='Blues',norm=norm):
    #axis will be a pyplot axis object, and patches will be a list of patches,
    # e.g., an element (single frame) of the list returned by read_all_patches 
    vlims = np.array([-0.05, 0.05])
    #Sort patches by level, from lowest to highest
    patches = sorted(patches, key=lambda p: p["level"])
    boundary_colors = ['k','k', "#D209E8"]+['#D209E8']*10 # black for level 1, gold for level 2, and purple for level 3 and above
    for idx_patch,patch in  enumerate(patches):
        level = patch["level"]
        if level in plot_cell_levels:
            axis.pcolormesh(
                patch["X"], patch["Y"], patch["ETA"],
                # vmin=vlims[0],
                # vmax=vlims[1],
                cmap=cmap,
                norm=norm,
                shading="flat",
                edgecolors =  boundary_colors[patch["level"]-1],
                linewidth=0.5,
                zorder = level,
                # rasterized=True
            )
        else:
            # First topography, then surface elevation
            bath = patch["BATH"]
            bath = np.ma.masked_invalid(bath)

            for i in range(len(topo_levels) - 1):
                lower = topo_levels[i]
                upper = topo_levels[i+1]

                mask = (
                    (bath >= lower) &
                    (bath < upper)
                )

                if np.any(mask):
                    const_field = np.ones_like(bath)
                    const_field = np.ma.masked_where(~mask, const_field)

                    axis.pcolormesh(
                        patch["X"], patch["Y"], const_field,
                        cmap=mcolors.ListedColormap([topo_colors[i]]),
                        shading="flat",
                        edgecolors="none",
                        linewidth=0,
                        antialiased=False,
                        zorder=level-0.5
                    )

            axis.pcolormesh(
                patch["X"], patch["Y"], patch["ETA"],
                # vmin=vlims[0],
                # vmax=vlims[1],
                cmap=cmap,
                norm=norm,
                shading="flat",
                zorder = level,
                # rasterized=True
            )

        if level in plot_boundary_levels:
            dx = patch["X"][1,0] - patch["X"][0,0]
            dy = patch["Y"][0,1] - patch["Y"][0,0]
            xmin = patch["X"].min()
            xmax = patch["X"].max()
            ymin = patch["Y"].min()
            ymax = patch["Y"].max()

            rect = Rectangle(
                (xmin, ymin),
                xmax - xmin,
                ymax - ymin,
                fill=False,
                edgecolor=boundary_colors[level-1],
                linewidth=1,
                zorder=10
            )
            axis.add_patch(rect)
    axis.set_aspect('equal')
    
        
    

if __name__ == "__main__":
    #####################################################################
    frame_nos = [6]
    print("Plotting 2D surface...")
    print("Reading in data and processing patches...")
    systems = [# 'HSGN',
            'mbHSGN',
            'SGN',
            'SWE']
    paths = {
        # 'HSGN': path_HSGN,
        'mbHSGN': path_mbHSGN,
        'SGN': path_SGN,
        'SWE': path_SWE
    }
    all_frames = {}
    for system in systems:
        path = paths[system]
        all_frames[system] = read_all_patches(frame_nos)

    fig,axes = plt.subplots(2,2, figsize=(15,11))
    print("Plotting...")
    axes = axes.flatten()
    for idx_system, system in enumerate(systems):
        patches = all_frames[system][0] # only one frame in this case
        fill_one_axis_surface_2D(axes[idx_system], patches,cmap="RdBu_r",norm=norm)
        axes[idx_system].set_title(f"{system}", fontsize=15)
        axes[idx_system].set_xlabel('longitude',fontsize=15)
        #Set limits
        axes[idx_system].set_xlim(132,175)
        axes[idx_system].set_ylim(14,45)
        #Plot transect line
        axes[idx_system].plot([163.454,170.631], [31.1458,26.5282], color='k', linewidth=1.0,zorder=20)
        axes[idx_system].text(170.6,26.45,'transect',ha='center',va='top',fontsize=14,zorder=20)
        gauges = {
            '21413': (152.1167, 30.5153),
            '21418': (148.694, 38.711)
        }
        for gauge_id, (lon, lat) in gauges.items():
            axes[idx_system].plot(lon, lat, marker='o', color='k', markersize=5, zorder=20)
            axes[idx_system].text(lon, lat+0.5, gauge_id, ha='center', va='bottom', fontsize=14, zorder=20)
        #Add individual colorbar for this axis
        vlims = np.array([-0.05, 0.05])
        sm = plt.cm.ScalarMappable(cmap="RdBu_r",norm=norm)
        sm.set_array([])
        # cbar = fig.colorbar(sm, ax=axes[idx_system], orientation='vertical', fraction=0.03, pad=0.04)
        cbar = fig.colorbar(sm, ax=axes[idx_system], orientation='vertical',
                    fraction=0.03, pad=0.04,
                    ticks=[-1, -0.1, 0, 0.1, 1])
        cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{x:g}'))
        cbar.ax.minorticks_off()   # kill the small log subdivision ticks
        axes[idx_system].set_ylabel('latitude', fontsize=15)

    #Last axis is a 1d plot
    frameno = 6
    method = 'linear'  # for interpolation from 2D output
    # from Popinet (personal comm.): (163.454,31.1458) to (170.631,26.5282)
    xout = np.linspace(163.454, 170.631, 1000)
    yout = np.linspace(31.1458, 26.5282, 1000)

    swe_color = 'blue'
    HSGN_color = 'green'
    SGN_color = 'magenta'
    mbHSGN_color =  '#F05800'

    outdir = '_output_UCSB_SWE'#'_output_swe_%s' % event.lower()
    framesoln = Solution(frameno, path=outdir, file_format='binary')
    etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout, method=method)
    axes[-1].plot(xout,etaout,color=swe_color,linewidth=1.8,label='SWE')

    outdir = '_output_UCSB_mbHSGN'#'_output_sgn_%s' % event.lower()
    framesoln = Solution(frameno, path=outdir, file_format='binary')
    etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout, method=method)
    axes[-1].plot(xout,etaout,color=mbHSGN_color,linewidth=2.,label='mbHSGN',linestyle='-')

    outdir = '_output_UCSB_HSGN'#'_output_sgn_%s' % event.lower()
    framesoln = Solution(frameno, path=outdir, file_format='binary')
    etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout, method=method)
    axes[-1].plot(xout,etaout,color=HSGN_color,linewidth=1.8,label='HSGN',linestyle='--')

    outdir = '_output_UCSB_SGN'#'_output_sgn_%s' % event.lower()
    framesoln = Solution(frameno, path=outdir, file_format='binary')
    etaout = gridtools.grid_output_2d(framesoln, -1, xout, yout, method=method)
    axes[-1].plot(xout,etaout,color=SGN_color,linewidth=2,label='SGN',linestyle=':')
    axes[-1].set_aspect(4./5.)   # y unit is twice x unit

    axes[-1].set_xlim(163,171)
    axes[-1].set_ylim(-3,4)
    axes[-1].set_yticks(range(-2,4,1))
    axes[-1].set_ylabel(r'$\eta$ (m)',fontsize=15)
    axes[-1].grid(True)
    #xticks(range(163,171,1),[])
    axes[-1].set_xticks(range(163,172,1))
    axes[-1].set_xlabel('longitude',fontsize=15)
    axes[-1].legend(loc='upper left',framealpha=1,fontsize=15)
    thours = framesoln.t / 3600.
    axes[-1].set_title(f'Surface elevation on transect at time {thours:.2f} hours', fontsize=15)



        
    # plt.tight_layout()
    plt.savefig("Figures/surface_2D.png", dpi=300)
    plt.close()




    