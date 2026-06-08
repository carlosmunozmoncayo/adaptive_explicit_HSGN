
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
"""

from __future__ import absolute_import
import numpy as np
from mapc2p import mapc2p
from clawpack.visclaw import gridtools  # would work with v5.7.0
from clawpack.visclaw import gaugetools
from clawpack.visclaw import particle_tools
from clawpack.visclaw import legend_tools


#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    plotdata.clearfigures() # clear any old figures,axes,items data

    print('Reading all gauges...')
    gauge_solutions = particle_tools.read_gauges(gaugenos='all', 
                                                 outdir=plotdata.outdir)

    # gaugenos_lagrangian = [k for k in gauge_solutions.keys() \
    #             if gauge_solutions[k].gtype=='lagrangian']
    # gaugenos_stationary = [k for k in gauge_solutions.keys() \
    #             if gauge_solutions[k].gtype=='stationary']
    # def add_particles(current_data):
    #     t = current_data.t

    #     # plot recent path:
    #     t_path_length = 2.   # length of path trailing particle
    #     kwargs_plot_path = {'linewidth':1, 'color':'k'}
    #     particle_tools.plot_paths(gauge_solutions, 
    #                               t1=t-t_path_length, t2=t, 
    #                               gaugenos=gaugenos_lagrangian, 
    #                               kwargs_plot=kwargs_plot_path)

    #     # plot current location:
    #     kwargs_plot_point = {'marker':'o','markersize':3,'color':'k'}
    #     particle_tools.plot_particles(gauge_solutions, t, 
    #                                   gaugenos=gaugenos_lagrangian, 
    #                                   kwargs_plot=kwargs_plot_point)  

    # plot any stationary gauges:
    # gaugetools.plot_gauge_locations(plotdata, \
    #         gaugenos=gaugenos_stationary, format_string='kx', add_labels=False)

    
    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [17,23]
    plotaxes.ylimits = [-3,3]
    plotaxes.title = 'Water depth'
    plotaxes.scaled = True
    # plotaxes.axescmd = 'axes([.15,.5,.7,.45])'


    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = 'turbo'#colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.9
    plotitem.pcolor_cmax = 1.2
    plotitem.add_colorbar = True

    plotitem.amr_celledges_show = [0,0,0]  
    plotitem.amr_patchedges_show = [1,1,1]
    plotitem.amr_patchedges_color = ['green','k']
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

   
    def plot_gauges_and_circle(current_data):
        gaugetools.plot_gauge_locations(current_data.plotdata,
             gaugenos='all', format_string='ko', add_labels=True)
        from pylab import plot
        # Plot cylinder boundary as a circle
        cylinder_center_x = 20.0
        cylinder_center_y = 0.0
        cylinder_radius = 1.5875
        theta = np.linspace(0, 2*np.pi, 100)
        x_circle = cylinder_center_x + cylinder_radius * np.cos(theta)
        y_circle = cylinder_center_y + cylinder_radius * np.sin(theta)
        plot(x_circle, y_circle, 'w-', linewidth=1)

    plotaxes.afteraxes = plot_gauges_and_circle

    # plotaxes.afteraxes = plot_xsec


    # # Figure for contour plot
    # plotfigure = plotdata.new_plotfigure(name='contour', figno=1)

    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = [0,1]
    # plotaxes.ylimits = [0,1]
    # plotaxes.title = 'Solution'
    # plotaxes.scaled = True
    # plotaxes.afteraxes = addgauges

    # # Set up for item on these axes:
    # plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    # plotitem.plot_var = 0
    # plotitem.contour_nlevels = 20
    # plotitem.contour_min = 0.01
    # plotitem.contour_max = 0.99
    # plotitem.amr_contour_colors = ['r','g','b']  # color on each level
    # plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    # plotitem.celledges_show = 0
    # plotitem.patchedges_show = 0


    # # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = [-10,10]
    # plotaxes.ylimits = [-1,1]
    plotaxes.title = 'Grid patches'
    plotaxes.scaled = False
    plotitem.add_colorbar = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#3B4CC0', "#4ADE80", "#EAB308"]
    plotitem.amr_celledges_show = [1,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
    plotitem.amr_patchedges_color = ['#AAB8FF', '#248F40', '#8C6B00']
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    # #-----------------------------------------
    # # Figures for gauges
    # #-----------------------------------------
    # plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
    #                 type='each_gauge')
    # plotfigure.clf_each_gauge = True

    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'q'

    # # Plot q as blue curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 0
    # plotitem.plotstyle = 'b-'

    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.print_gaugenos = 'all'#[1,2,3,4,5,6]        # list of gauges to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata

def plot_mass(current_data):
    return np.sqrt(current_data.q[0,:,:]**2 + current_data.q[3,:,:]**2)

def plot_streamlines(current_data):
    from clawpack.pyclaw import Solution
    from pylab import quiver
    pd = current_data.plotdata
    frameno = current_data.frameno
    framesoln = Solution(frameno, path=pd.outdir, file_format=pd.format)
    xout = np.linspace(-14.,30.,100)
    yout = np.linspace(-6.,6.,100)
    XOUT, YOUT = np.meshgrid(xout,yout, indexing='ij')
    psi_real = gridtools.grid_output_2d(framesoln, 0, XOUT, YOUT)
    psi_imag = gridtools.grid_output_2d(framesoln, 3, XOUT, YOUT)
    u_real = gridtools.grid_output_2d(framesoln, 1, XOUT, YOUT) #x-velocity
    v_real = gridtools.grid_output_2d(framesoln, 2, XOUT, YOUT) #y-velocity
    u_imag = gridtools.grid_output_2d(framesoln, 4, XOUT, YOUT) #x-velocity
    v_imag = gridtools.grid_output_2d(framesoln, 5, XOUT, YOUT) #y-velocity
    #Computing superfluid velocity field
    # reconstruct complex fields
    psi = psi_real + 1j * psi_imag
    u = u_real + 1j * u_imag  
    v = v_real + 1j * v_imag   

    # compute density
    rho = np.abs(psi)**2

    # compute current components Jx, Jy
    Jx = np.imag(np.conj(psi) * u)
    Jy = np.imag(np.conj(psi) * v)

    # superfluid velocity (avoid division by small rho)
    eps = 1e-12
    vx = np.where(rho > eps, Jx / rho, 0.0)
    vy = np.where(rho > eps, Jy / rho, 0.0)

    ######################################
    # Plotting velocity field using quiver
    # quiver(XOUT, YOUT, vx, vy, color='k', scale=50)

    ######################################
    #Using streamplot instead of quiver 
    # from pylab import streamplot
    # # streamplot(XOUT, YOUT, vx, vy, color='k', density=1.5, linewidth=0.5, arrowsize=1)
    # streamplot(
    # XOUT.T, YOUT.T, vx.T, vy.T,
    # density=1.0, color=np.log10(rho+1e-12),linewidth=1, arrowsize=1)

    ######################################
    #Normalized quiver only around vortices'
    rho0 = 0.1
    mask = rho < rho0           # boolean mask for where to plot arrows
    mask = np.logical_and(mask, XOUT>22)

    # normalize (so all arrows same length)
    mag = np.sqrt(vx**2 + vy**2)
    vx_dir = np.zeros_like(vx)
    vy_dir = np.zeros_like(vy)
    vx_dir[mask] = vx[mask] / mag[mask]
    vy_dir[mask] = vy[mask] / mag[mask]

    # note the .T transposes (since your arrays use 'ij' indexing)
    quiver(
        XOUT.T[mask.T],
        YOUT.T[mask.T],
        vx_dir.T[mask.T],
        vy_dir.T[mask.T],
        color='k',width=0.003
    )



def addgauges(current_data):
    from clawpack.visclaw import gaugetools
    gaugetools.plot_gauge_locations(current_data.plotdata, \
         gaugenos='all', format_string='ko', add_labels=True)
