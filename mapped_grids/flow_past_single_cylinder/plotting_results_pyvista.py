import numpy as np
#Import plotting libraries, still deciding which one to use
from matplotlib import pyplot as plt
import plotly.graph_objects as go
from matplotlib.colors import LightSource
import pyvista as pv

from clawpack.pyclaw import Solution
from map_module_f2py import map_module
mapc2p = map_module.mapc2p
vectorized_mapc2p = map_module.vectorized_mapc2p

def edges_from_centers(c):
    # uniform spacing within a patch
    d = c[1] - c[0]
    return np.concatenate(([c[0] - d/2], c + d/2))   # length n+1

def centers_to_corners(Z):
    """(nx,ny) cell-centered -> (nx+1,ny+1) corner values,
    bilinear interior + linear extrapolation at the boundary."""
    nx, ny = Z.shape
    Zp = np.empty((nx+2, ny+2))
    Zp[1:-1, 1:-1] = Z
    Zp[0,    1:-1] = 2*Z[0, :]  - Z[1, :]
    Zp[-1,   1:-1] = 2*Z[-1,:]  - Z[-2,:]
    Zp[1:-1, 0]    = 2*Z[:, 0]  - Z[:, 1]
    Zp[1:-1, -1]   = 2*Z[:, -1] - Z[:, -2]
    Zp[0,0]   = 2*Zp[1,0]   - Zp[2,0]
    Zp[-1,0]  = 2*Zp[-2,0]  - Zp[-3,0]
    Zp[0,-1]  = 2*Zp[1,-1]  - Zp[2,-1]
    Zp[-1,-1] = 2*Zp[-2,-1] - Zp[-3,-1]
    return 0.25*(Zp[:-1,:-1] + Zp[1:,:-1] + Zp[:-1,1:] + Zp[1:,1:])


def read_patches(frame_no, remove_overlaps=True, remove_cylinder=True,
                 variable_index=0, offset_z=-1.0):
    patches = []
    sol = Solution()
    sol.read(frame = frame_no, read_aux=True)

    for state in sol.states:
        xc = state.grid.x.centers
        yc = state.grid.y.centers
        xe = edges_from_centers(xc)
        ye = edges_from_centers(yc)
        XE, YE = np.meshgrid(xe, ye, indexing='ij')          # (nx+1, ny+1)
        X, Y = vectorized_mapc2p(XE, YE)

        qc = state.q[variable_index, :, :] + offset_z         # (nx, ny)
        Z = centers_to_corners(qc)                            # (nx+1, ny+1)

        patches.append({
            'X': X, 'Y': Y, 'Z': Z,
            'level': state.patch.level,
            'xc': xc, 'yc': yc,                               # computational centers
            'xlim': (xe[0], xe[-1]),                          # computational extent
            'ylim': (ye[0], ye[-1]),
            'keep': np.ones((len(xc), len(yc)), dtype=bool),  # per-cell mask
        })

    if remove_overlaps:
        for p in patches:
            XC, YC = np.meshgrid(p['xc'], p['yc'], indexing='ij')   # (nx, ny) centers
            for q in patches:
                if q['level'] > p['level']:
                    x0, x1 = q['xlim']; y0, y1 = q['ylim']
                    inside = (XC > x0) & (XC < x1) & (YC > y0) & (YC < y1)
                    p['keep'] &= ~inside

    if remove_cylinder:
        cx, cy, r = 20.0, 0.0, 1.59
        for p in patches:
            XC, YC = np.meshgrid(p['xc'], p['yc'], indexing='ij')
            Xc, Yc = vectorized_mapc2p(XC, YC)                       # physical centers
            inside = (Xc - cx)**2 + (Yc - cy)**2 <= r**2
            p['keep'] &= ~inside

    return patches


def plot_amr_surface(
    patches,
    elev=30, azim=45,
    z_exag=1.0,
    level_colors=None,
    screenshot=None,
    show=True,
    window_size=(2500, 1750)
):
    """
    patches: list of dicts with X,Y,Z,level
    elev, azim: camera angles in degrees
    z_exag: vertical exaggeration factor
    level_colors: dict {level: (r,g,b)} or list by level index
    screenshot: filename to save static image
    """

    # Default colors (tab10-like)
    if level_colors is None:
        level_colors = {
            1: (31, 119, 180),
            2: (255, 127, 14),
            3: (44, 160, 44),
        }


    plotter = pv.Plotter(off_screen=not show,window_size=window_size)
    ###############################################



    ###############################################
    #Plot patches
    grids = []
    for p in patches:
        X, Y, Z = p["X"], p["Y"], p["Z"] * z_exag
        level = p["level"]

        grid = pv.StructuredGrid(X, Y, Z)
        grid.cell_data["keep"] = p["keep"].ravel(order="F").astype(np.uint8)
        grid = grid.threshold(0.5, scalars="keep")   # keep == 1, drop covered cells
        grids.append(grid)

        color = tuple(c / 255 for c in level_colors[level])
        plotter.add_mesh(grid, color=color, show_edges=False, smooth_shading=False)
        ###############################################
    #Add back line at the outer boundary of the surface
    # Merge patches
    combined = grids[0]
    for g in grids[1:]:
        combined = combined.merge(g)

    # Extract outer boundary
    edges = combined.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=False,
    )
    

    # Draw outline
    plotter.add_mesh(edges, color="black", line_width=2)
    ###############################################
    
    #Add cylinder
    z_min, z_max = -0.2*z_exag, 0.5*z_exag
    height = z_max - z_min
    z_center = 0.5 * (z_min + z_max)

    cyl = pv.Cylinder(
        center=(20.0, 0.0, z_center),
        direction=(0, 0, 1),
        radius=1.59,
        height=height,
        resolution=200,
    )

    plotter.add_mesh(cyl, color="gray", opacity=1.0,smooth_shading=True)
    ###############################################

    # ---- Camera ----
    plotter.view_isometric()   # good default
    plotter.camera.elevation = elev
    plotter.camera.azimuth = azim

    # Optional: equal horizontal scale, exaggerated vertical
    plotter.set_scale(zscale=z_exag)

    plotter.show_bounds(
        grid='back',
        location='outer',
        xtitle="x [m]",
        ytitle="y [m]",
        ztitle="η [m]",       # surface elevation symbol
        font_size=int(28*window_size[0]/1000),      # scales tick numbers
        # font_size=font_size_scaled,      # scales tick numbers
        # title_font_size=font_size_scaled + 6,   # scales axis titles
        color="black",
        ticks="outside",
    )



    # ---- Render ----
    if screenshot:
        plotter.show(screenshot=screenshot)
    elif show:
        plotter.show()

    return plotter


# frame_list = range(0, 25, 1)  # frames 0, 1, 2, ..., 24
frame_list = range(0, 25, 1)
for frame in frame_list:
    print(f"Plotting frame {frame}...")
    patches = read_patches(frame)
    #List frames = [0,13,24]

    colors = {
        # 1: (0, 0, 255),     # blue
        # 2: (255, 0, 0),     # red
        # 3: (0, 255, 0),     # green
    1: (179, 157, 219), # lavender
    2: (129, 216, 208), # turquoise
    3: (251, 180, 174) # coral
    # 1: (128, 128, 0), # olive
    # 2: (0, 200, 200), # cyan
    # 3: (200, 0, 200) # magenta
    }

    plot_amr_surface(
        patches,
        elev=20,
        azim=-180,
        z_exag=10,          # huge vertical exaggeration
        level_colors=colors,
        show=False,
        screenshot=f"Figures/surface{frame}.png",   # comment out to open window
    )