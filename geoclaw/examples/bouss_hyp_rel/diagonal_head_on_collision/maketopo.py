
"""
Module to create topo and qinit data files for this example.
"""

from __future__ import absolute_import
from clawpack.geoclaw.topotools import Topography
from numpy import *
import numpy as np

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 2000
    nypoints = 2000
    xlower = -1000
    xupper = 1000
    ylower = -1000
    yupper = 1000
    outfile= "flat.tt3"     

    topography = Topography(topo_func=topo)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=3, Z_format="%22.15e")

def topo(x,y):
    base_topo = -1.0
    z = base_topo*np.ones(x.shape)
    return z


if __name__=='__main__':
    maketopo()
