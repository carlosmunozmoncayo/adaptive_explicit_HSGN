import numpy as np
import matplotlib.pyplot as plt
from map_module_f2py import map_module
mapc2p = map_module.mapc2p
vectorized_mapc2p = map_module.vectorized_mapc2p


plot_horizontal = False
grid_data = np.loadtxt("single_coarse_grid.txt")
output_dir = "Figures"


print("Reading in data...")
# load data (i, j, aux1,...,aux7)
data = np.loadtxt("single_coarse_grid.txt")
# Split into pieces
xc, yc = data[:,0], data[:,1]      # computational coords
xp, yp = data[:,2], data[:,3]      # physical coords
q      = data[:,4:9]                 # solution variables
aux    = data[:,9:]                # everything after is aux

xc_temp = np.unique(xc)
dx = xc_temp[1] - xc_temp[0]
xc_edges = xc_temp - 0.5*dx
xc_edges = np.append(xc_edges, xc_edges[-1] + dx)  # add last edge
yc_temp = np.unique(yc)
dy = yc_temp[1] - yc_temp[0]
yc_edges = yc_temp - 0.5*dy
yc_edges = np.append(yc_edges, yc_edges[-1] + dy)  # add last edge

XC, YC = np.meshgrid(xc_edges, yc_edges, indexing='ij')
XP, YP = vectorized_mapc2p(XC, YC)

yp_edges = np.zeros_like(yc_edges)
xp_edges = np.zeros_like(xc_edges)
print("Mapping edges to physical space...")
for i in range(len(xc_edges)):
    for j in range(len(yc_edges)):
        xp_edges[i], yp_edges[j] = mapc2p(xc_edges[i], yc_edges[j])

print("Plotting...")
if plot_horizontal:
    # Plotting the grid
    plt.figure(figsize=(10, 10))
    #Plot a red circle of radius 1.5875 centered at (20,0)
    circle = plt.Circle((20, 0), 1.5875, color='gray', fill=True, lw=2)
    plt.gca().add_artist(circle)
    for i in range(XP.shape[0]):
        plt.plot(XP[i, :], YP[i, :], 'k-',lw=1)  # Horizontal lines
    for j in range(XP.shape[1]):
        plt.plot(XP[:, j], YP[:, j], 'k-',lw=1)  # Vertical lines
    plt.title('Mapped Grid in Physical Domain')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect('equal')
    plt.xlim(12,25)
    plt.ylim(-5,5)

    cylinder_center_x = 20.0
    cylinder_center_y = 0.0
    cylinder_radius = 1.5875
    #Everywhere at most level 3
    rtilde_list = [4.5, 2.61, 1.66, 1.03]
    for i,rtilde in enumerate(rtilde_list):
        x = cylinder_center_x - rtilde*cylinder_radius
    #Also, we have some gauges at 150 degrees from the origin. So top-right quadrant if the cylinder is at the origin.
    #For these we have r such that r/R = 2.92, 2.29, 1.35, 1.03

    gauge_colors = ['red', 'blue', 'green', 'orange']
    rtilde_list = [4.5, 2.61, 1.66, 1.03]
    for i,rtilde in enumerate(rtilde_list):
        x = cylinder_center_x - rtilde*cylinder_radius
        y = cylinder_center_y
        plt.scatter(x, y, c='red', s=70, edgecolors='black', zorder=5)
        #Place a marker at the gauge location with a different color for each gauge, with label next to it
        # plt.scatter(x, y, c=gauge_colors[i], s=50, label=f'Gauge {i+1}', edgecolors='black', zorder=5)
        plt.text(
            x-0.25, y - 0.5,  f'WG {i+1}',   # Offset slightly
            fontsize=10,
            bbox=dict(
                facecolor='white',      # Background color
                edgecolor='black',      # Border color
                boxstyle='round,pad=0.3'
            )
        )


    def gauge_location(rtilde,R):
        c = np.cos(30*np.pi/180.)
        s = np.sin(30*np.pi/180.)
        r = rtilde*R
        x = cylinder_center_x + r*c
        y = cylinder_center_y + r*s
        return x,y
    rtilde_list = [2.92, 2.29, 1.35, 1.15]
    for i,rtilde in enumerate(rtilde_list):
        x,y = gauge_location(rtilde,cylinder_radius)
        #scatter plot the gauge location
        plt.scatter(x, y, c='red', s=70, edgecolors='black', zorder=5)
        if i==2:
            plt.text(
                x-0.25, y + 0.35,  f'WG {i+5}',   # Offset slightly
                fontsize=10,
                bbox=dict(
                    facecolor='white',      # Background color
                    edgecolor='black',      # Border color
                    boxstyle='round,pad=0.3'
                )
            )
        else:
            plt.text(
                x-0.25, y - 0.5,  f'WG {i+5}',   # Offset slightly
                fontsize=10,
                bbox=dict(
                    facecolor='white',      # Background color
                    edgecolor='black',      # Border color
                    boxstyle='round,pad=0.3'
                )
            )

    plt.savefig(f'{output_dir}/mapped_grid_with_gauges.pdf', bbox_inches='tight')

else:
    #The same as above, but rotated, so that the figure is larger across the vertical direcation (e.g, x axis is vertical and y axis is horizontal)
    # Plotting the grid
    plt.figure(figsize=(10, 10))
    #Plot a red circle of radius 1.5875 centered at (20,0)
    circle = plt.Circle((0, 20), 1.5875, color='gray', fill=True, lw=2)
    plt.gca().add_artist(circle)
    for i in range(XP.shape[0]):
        plt.plot( YP[i, :],XP[i, :], 'k-',lw=1)  # Horizontal lines
    for j in range(XP.shape[1]):
        plt.plot(YP[:, j],XP[:, j], 'k-',lw=1)  # Vertical lines
    plt.title('Mapped Grid in Physical Domain', fontsize=16)
    plt.xlabel('y')
    plt.ylabel('x')
    plt.gca().set_aspect('equal')
    plt.ylim(12.5,24.2)
    plt.xlim(-4,4)


    cylinder_center_x = 20.0
    cylinder_center_y = 0.0
    cylinder_radius = 1.5875
    #Everywhere at most level 3
    rtilde_list = [4.5, 2.61, 1.66, 1.03]
    for i,rtilde in enumerate(rtilde_list):
        x = cylinder_center_x - rtilde*cylinder_radius
    #Also, we have some gauges at 150 degrees from the origin. So top-right quadrant if the cylinder is at the origin.
    #For these we have r such that r/R = 2.92, 2.29, 1.35, 1.03

    gauge_colors = ['red', 'blue', 'green', 'orange']
    rtilde_list = [4.5, 2.61, 1.66, 1.03]
    for i,rtilde in enumerate(rtilde_list):
        x = cylinder_center_x - rtilde*cylinder_radius
        y = cylinder_center_y
        plt.scatter(y, x, c='red', s=70, edgecolors='black', zorder=5)
        #Place a marker at the gauge location with a different color for each gauge, with label next to it
        # plt.scatter(x, y, c=gauge_colors[i], s=50, label=f'Gauge {i+1}', edgecolors='black', zorder=5)
        if i==0:
            plt.text(
                y-0.25, x + 0.35,  f'WG {i+1}',   # Offset slightly
                fontsize=12,
                bbox=dict(
                    facecolor='white',      # Background color
                    edgecolor='black',      # Border color
                    boxstyle='round,pad=0.3'
                )
            )
        else:
            plt.text(
                y-0.25, x - 0.5,  f'WG {i+1}',   # Offset slightly
                fontsize=12,
                bbox=dict(
                    facecolor='white',      # Background color
                    edgecolor='black',      # Border color
                    boxstyle='round,pad=0.3'
                )
            )


    def gauge_location(rtilde,R):
        c = np.cos(30*np.pi/180.)
        s = np.sin(30*np.pi/180.)
        r = rtilde*R
        x = cylinder_center_x + r*c
        y = cylinder_center_y + r*s
        return x,y
    rtilde_list = [2.92, 2.29, 1.35, 1.15]
    for i,rtilde in enumerate(rtilde_list):
        x,y = gauge_location(rtilde,cylinder_radius)
        #scatter plot the gauge location
        plt.scatter(y, x, c='red', s=70, edgecolors='black', zorder=5)
        if i==2:
            plt.text(
                y-0.25, x + 0.35,  f'WG {i+5}',   # Offset slightly
                fontsize=12,
                bbox=dict(
                    facecolor='white',      # Background color
                    edgecolor='black',      # Border color
                    boxstyle='round,pad=0.3'
                )
            )
        else:
            plt.text(
                y-0.25, x - 0.5,  f'WG {i+5}',   # Offset slightly
                fontsize=12,
                bbox=dict(
                    facecolor='white',      # Background color
                    edgecolor='black',      # Border color
                    boxstyle='round,pad=0.3'
                )
            )

    #Axis labels and ticks fontsize
    plt.xlabel(r'$y(m)$', fontsize=14)
    plt.ylabel(r'$x(m)$', fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=12)

    plt.savefig(f'{output_dir}/mapped_grid_with_gauges.pdf', bbox_inches='tight')
