import numpy as np
import matplotlib.pyplot as plt
from clawpack.pyclaw import Solution
from clawpack.visclaw import gridtools 
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors

norm = mcolors.Normalize(vmin=-0.05, vmax=0.05)
path = './_output_mbHSGN'  # run_automated_2d.sh writes mbHSGN (EDC, Bouss_Hyp_rel=-1) here

def read_all_patches(frame_nos, remove_overlaps=True,
                     variable_index=0, offset_z=0,
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
            patch['Z'] = state.q[variable_index,:,:]+offset_z 
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
                        p["Z"][mask] = np.nan
        all_patches.append(patches)
    return all_patches

def fill_one_axis_surface_2D(axis,patches,plot_boundary_levels=[2,3],
                             plot_cell_levels=[1],cmap="RdBu_r",norm=norm):
    #axis will be a pyplot axis object, and patches will be a list of patches,
    # e.g., an element (single frame) of the list returned by read_all_patches 
    vlims = np.array([-0.05, 0.05])
    #Sort patches by level, from lowest to highest
    patches = sorted(patches, key=lambda p: p["level"])
    boundary_colors = ['k','gold', "#D209E8"]
    for idx_patch,patch in  enumerate(patches):
        level = patch["level"]
        if level in plot_cell_levels:
            axis.pcolormesh(
                patch["X"], patch["Y"], patch["Z"],
                cmap=cmap,
                norm=norm,
                shading="flat",
                edgecolors =  boundary_colors[patch["level"]-1],
                linewidth=0.5,
                zorder = level
            )
        else:
            axis.pcolormesh(
                patch["X"], patch["Y"], patch["Z"],
                cmap=cmap,
                norm=norm,
                shading="flat",
                zorder = level
            )

        if level in plot_boundary_levels:
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
    frame_nos = [4, 8, 16]
    print("Plotting 2D surface...")
    print("Reading in data and processing patches...")
    all_frames = read_all_patches(frame_nos,variable_index=-1)
    fig,axes = plt.subplots(1,3, figsize=(14,4))
    print("Plotting...")
    for idx_frame, patches in enumerate(all_frames):
        fill_one_axis_surface_2D(axes[idx_frame], patches,cmap="RdBu_r",norm=norm)
        axes[idx_frame].set_title(f"t={frame_nos[idx_frame]}s")
        axes[idx_frame].set_xlabel('x (m)')
        if idx_frame>0:
            axes[idx_frame].set_yticklabels([])
    axes[0].set_ylabel('y (m)')
    #Add a common colorbar
    vlims = np.array([-0.075, 0.075])
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.015, pad=0.04)
    cbar.set_label(r'$\eta$ (m)', rotation=270, labelpad=15)

        
    # plt.tight_layout()
    plt.savefig("Figures/surface_2D.png", dpi=300)
    plt.close()
    # plt.show()
    #####################################################################
    print("Plotting 1D slice...")
    frame_nos = [16]
    fig,axes = plt.subplots(1,1, figsize=(12,3))
    
    sol_1d = Solution()
    sol_1d.read(frame=frame_nos[-1], path="1d_radial/_output_EDC_paper",
                      read_aux=False,file_format='ascii')
    eta_1d = sol_1d.state.q[-1,:]
    xc_1d = sol_1d.state.grid.c_centers[0]

    sol_2d = Solution()
    sol_2d.read(frame=frame_nos[-1], path=path,
                      read_aux=False,file_format='ascii')
    
    xc_2d = np.linspace(0.,5000.,5000)
    eta_2d = gridtools.grid_output_2d(sol_2d, -1, xc_2d, xc_2d)
    xc_2d = np.sqrt(2)*xc_2d # convert to radial distance along x=y

    axes.plot(xc_1d, eta_1d, label="1D HSGN with forcing source", linewidth=2, color='k')
    axes.plot(xc_2d, eta_2d, label=f"2D HSGN across $x=y$", linewidth=2,color='green')
    axes.set_title(f"t={frame_nos[-1]}s")

    axes.set_xlabel('radial distance (m)')
    axes.set_ylabel(r'$\eta (m)$')
    plt.legend()
    plt.tight_layout()
    plt.savefig("Figures/surface_1D_slice.pdf", dpi=300)




    