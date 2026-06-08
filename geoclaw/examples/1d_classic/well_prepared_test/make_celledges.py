"""
Make piecewise linear topography for wave tank.
"""

from pylab import *
from clawpack.geoclaw import nonuniform_grid_tools

# xlower = -70. 
# xupper = 70

xlower = -1e4
xupper = 1e4

#Flat topo:
xzpairs = [(-1e4, -100.),       # left edge
           ( 1e4, -100.)]      # right edge
        #    ( 23.15, 0.)]      # right edge (physical)

topo_fcn = nonuniform_grid_tools.make_pwlin_topo_fcn(xzpairs)

mx = 30000#14000    #flag_grid_mx PLEASE DO NOT DELETE THIS COMMENT
hmin = 2.0  # use uniform grid in shallower water

nonuniform_grid_tools.make_celledges_cfl(xlower, xupper, mx, topo_fcn,
        hmin, fname='celledges.data', plot_topo=True)

