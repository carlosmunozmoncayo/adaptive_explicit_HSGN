import os
import numpy as np
import matplotlib.pyplot as plt
from clawpack.pyclaw import Solution
from clawpack.visclaw import gridtools 
from matplotlib.patches import Rectangle
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p

path = './_output'

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

def fill_one_axis_surface_2D(axis,patches,plot_boundary_levels=[3],
                             plot_cell_levels=[],
                             boundary_colors = ['k','gold', "#D209E8"]+10*["#FF0000"],
                             skip_levels = [], cmap='Blues'):
    #axis will be a pyplot axis object, and patches will be a list of patches,
    # e.g., an element (single frame) of the list returned by read_all_patches 
    vlims = np.array([0, 0.25])
    #Sort patches by level, from lowest to highest
    patches = sorted(patches, key=lambda p: p["level"])
    for idx_patch,patch in  enumerate(patches):
        level = patch["level"]
        if level in skip_levels:
            continue
        if level in plot_cell_levels:
            axis.pcolormesh(
                patch["X"], patch["Y"], patch["Z"],
                vmin=vlims[0],
                vmax=vlims[1],
                cmap=cmap,
                shading="flat",
                edgecolors =  boundary_colors[patch["level"]-1],
                linewidth=0.5,
                zorder = level
            )
        else:
            axis.pcolormesh(
                patch["X"], patch["Y"], patch["Z"],
                vmin=vlims[0],
                vmax=vlims[1],
                cmap=cmap,
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
                linewidth=0.5,
                zorder=10
            )
            axis.add_patch(rect)
    axis.set_aspect('equal')

def add_arrows(axis, arrow_starts, arrow_ends, color='black'):
    for start, end in zip(arrow_starts, arrow_ends):
        axis.annotate(
            '', xy=end, xytext=start,
            arrowprops=dict(facecolor=color, edgecolor=color, shrink=0.05, width=1, headwidth=8),
            zorder=20
        )

def add_diagonal(axis, color='black'):
    axis.plot([-400, 400], [-400, 400], color=color, linewidth=1, linestyle='--', zorder=15)

def plot_for_animation(frame_nos, plot_boundary_levels=[3], plot_cell_levels=[], cmap='Blues'):
    all_patches = read_all_patches(frame_nos,variable_index=-1)
    for idx_frame, patches in enumerate(all_patches):
        print(f"Plotting frame {idx_frame+1} of {len(all_patches)}...")
        fig, ax = plt.subplots(figsize=(6,6))
        fill_one_axis_surface_2D(ax, patches,
                                plot_boundary_levels=plot_boundary_levels,
                                plot_cell_levels=plot_cell_levels,
                                cmap=cmap)
        ax.set_title(f"t={frame_nos[idx_frame]}s")
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        plt.savefig(f'Figures/sequences/frame_{idx_frame:04d}.png', dpi=300)
        plt.close()
    

if __name__ == "__main__":
    #####################################################################
    #####################################################################
    #Plot surface in 2D, single plot
    frame_nos = [0,40,80]#list(range(0,81,9))#[4+4*i for i in range(9)]  # Generate frame numbers from 0 to 100 in steps of 10
    print("Plotting 2D surface in single figure...")
    print("Reading in data and processing patches...")
    all_frames = read_all_patches(frame_nos,variable_index=-1)
    fig,axes = plt.subplots(1,3, figsize=(14,4))
    axes = axes.flatten()
    print("Plotting...")
    for idx_frame, patches in enumerate(all_frames):
        boundary_colors = 10*['k']
        boundary_colors[4] = '#E69F00'
        boundary_colors[5] = "#DE00E6"
        fill_one_axis_surface_2D(axes[idx_frame], patches,
                                plot_boundary_levels=[5,6],
                                plot_cell_levels=[],
                                skip_levels = [1,2,3,4],
                                boundary_colors = boundary_colors,
                                cmap='Blues')
        axes[idx_frame].set_title(f"t={frame_nos[idx_frame]}s")
        axes[idx_frame].set_xlabel('x (m)')
        if idx_frame>0:
            axes[idx_frame].set_yticklabels([])
    axes[0].set_ylabel('y (m)')
    #Add arrows to first plot to indicate direction of wave propagation
    arrow_starts = [(-325, -325), (-150, 0), (150, 0), (325, 325)]
    arrow_ends = [(-395, -395), (-80,70), (80, -70), (395, 395)]
    add_arrows(axes[0], arrow_starts, arrow_ends, color='black')
    add_diagonal(axes[0], color='black')
    #Add a common colorbar
    vlims = np.array([0, 0.25])
    sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=vlims[0], vmax=vlims[1]))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.012, pad=0.04)
    cbar.set_label(r'$\eta$ (m)', rotation=270, labelpad=15)
    plt.savefig("Figures/surface_2D.png", dpi=300)
    plt.close()
    print("Plotting 1D slice in single figure...")
    fig,all_axes = plt.subplots(1,3, figsize=(14,3))
    for idx_frame, frame_no in enumerate(frame_nos):
        axes = all_axes[idx_frame]
        #Read and plot 1D SGN solution 
        outdir1 = "1d/_output"
        mapc2p1, mx_edge, xp_edge = make_mapc2p(os.path.join(outdir1,'celledges.data'))
        sol_1d = Solution()
        sol_1d.read(frame=frame_no, path=outdir1,
                          read_aux=False,file_format='ascii')
        xc_1d = mapc2p1(sol_1d.state.grid.c_centers[0])
        eta_1d = sol_1d.state.q[-1,:]
        axes.plot(xc_1d, eta_1d, label="1D SGN "+r"($\Delta x \approx 0.01 m$)", linewidth=1, color='m')


        #Read and plot 2D solution along x=y
        sol_2d = Solution()
        sol_2d.read(frame=frame_no, path='_output',
                        read_aux=False,file_format='ascii')
        
        xc_2d = np.linspace(-400,400,5000)
        eta_2d = gridtools.grid_output_2d(sol_2d, -1, xc_2d, xc_2d)
        xc_2d = np.sqrt(2)*xc_2d # convert to radial distance along x=y
        axes.plot(xc_2d, eta_2d, label=f"2D HSGN "+r"($\alpha =3$, $\Delta x \approx 0.35m$)", linewidth=1,color='green')
        axes.set_title(f"t={frame_no}s")

        axes.set_xlabel(r'$\sqrt{2} x(m)$')
        
    # plt.legend()
    all_axes[0].set_ylabel(r'$\eta (m)$')
    all_axes[1].legend(loc='center')
    plt.tight_layout()
    plt.savefig("Figures/surface_1D_slice.pdf")
    plt.close()
    #####################################################################
    #####################################################################
    # #Plot a sequence of frames, with only the patch boundaries, to show the grid evolution.
    # print("Plotting sequence of frames...")
    # frame_nos = list(range(0,92,2))  # Generate frame numbers from 0 to 90 in steps of 2
    # plot_for_animation(frame_nos, plot_boundary_levels=[6], plot_cell_levels=[])
    #####################################################################
    #####################################################################
    # print("Plotting 1D slice...")
    # frame_nos = [30]
    # fig,axes = plt.subplots(1,1, figsize=(12,3))

    # sol_2d = Solution()
    # sol_2d.read(frame=frame_nos[-1], path='_output',
    #                   read_aux=False,file_format='ascii')
    
    # xc_2d = np.linspace(-400,400,5000)
    # eta_2d = gridtools.grid_output_2d(sol_2d, -1, xc_2d, xc_2d)
    # xc_2d = np.sqrt(2)*xc_2d # convert to radial distance along x=y

    # # axes.plot(xc_1d, eta_1d, label="1D HSGN with forcing source", linewidth=2, color='k')
    # axes.plot(xc_2d, eta_2d, label=f"2D HSGN across $x=y$", linewidth=2,color='green')
    # axes.set_title(f"t={frame_nos[-1]}s")

    # axes.set_xlabel('radial distance (m)')
    # axes.set_ylabel(r'$\eta (m)$')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig("Figures/surface_1D_slice.pdf", dpi=300)
    #####################################################################
    #####################################################################
    # print("Plotting sqequence of 1D slice...")
    # frame_nos = [2*i for i in range(46)]  # Generate frame numbers from 0 to 90 in steps of 2
    # for idx_frame, frame_no in enumerate(frame_nos):
    #     print(f"Plotting frame {idx_frame+1} of {len(frame_nos)}...")
    #     fig,axes = plt.subplots(1,1, figsize=(15,8))

    #     #Read and plot 1D SGN solution 
    #     outdir1 = "1d/_output"
    #     mapc2p1, mx_edge, xp_edge = make_mapc2p(os.path.join(outdir1,'celledges.data'))
    #     sol_1d = Solution()
    #     sol_1d.read(frame=frame_no, path=outdir1,
    #                       read_aux=False,file_format='ascii')
    #     xc_1d = mapc2p1(sol_1d.state.grid.c_centers[0])
    #     eta_1d = sol_1d.state.q[-1,:]
    #     axes.plot(xc_1d, eta_1d, label="1D SGN "+r"($\Delta x \approx 0.01 m$)", linewidth=2, color='m')


    #     #Read and plot 2D solution along x=y
    #     sol_2d = Solution()
    #     sol_2d.read(frame=frame_no, path='_output',
    #                     read_aux=False,file_format='ascii')
        
    #     xc_2d = np.linspace(-400,400,5000)
    #     eta_2d = gridtools.grid_output_2d(sol_2d, -1, xc_2d, xc_2d)
    #     xc_2d = np.sqrt(2)*xc_2d # convert to radial distance along x=y
    #     axes.plot(xc_2d, eta_2d, label=f"2D HSGN "+r"($\alpha =3$, $\Delta x \approx 0.25$)", linewidth=2,color='green')
    #     axes.set_title(f"t={frame_no}s")

    #     axes.set_xlabel('x=y(m)')
    #     axes.set_ylabel(r'$\eta (m)$')
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.savefig(f'Figures/slices/frame_{idx_frame:04d}.png', dpi=300)
    #####################################################################
    #####################################################################



    