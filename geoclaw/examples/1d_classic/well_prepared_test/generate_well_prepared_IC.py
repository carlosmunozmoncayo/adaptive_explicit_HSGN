import numpy as np
import matplotlib.pyplot as plt
from clawpack.geoclaw.nonuniform_grid_tools import make_mapc2p
from clawpack.pyclaw import Solution
#Import gaussian filter for smoothing
from scipy.ndimage import gaussian_filter1d

mapc2p = make_mapc2p('celledges.data')

def get_w(sol):
    xc = mapc2p[0](sol.grid.x.centers)
    dx =xc[1] - xc[0]
    eta = sol.q[-1,:]
    h = sol.q[0,:]
    bath = eta - h
    hu = sol.q[1,:]
    u = np.zeros_like(hu)
    mask = h > 1e-3
    u[mask] = hu[mask]/h[mask]
    divu = np.zeros_like(u)
    divu[1:-1] = (u[2:] - u[:-2])/(2*dx)

    divbath = np.zeros_like(bath)
    divbath[1:-1] = (bath[2:] - bath[:-2])/(2*dx)
    w = -0.5*h*divu + u*divbath
    return w

def get_p(sol0,sol1,sol2):
    xc = mapc2p[0](sol0.grid.x.centers)
    dx =xc[1] - xc[0]
    h0 = sol0.q[0,:]
    h1 = sol1.q[0,:]
    h2 = sol2.q[0,:]
    w0 = get_w(sol0)
    w1 = get_w(sol1)
    w2 = get_w(sol2)
    dt = sol1.t - sol0.t
    d_hw_dt = (w1*h1 - w0*h0)/dt
    d_hw_dt_2nd_order = (4*w1*h1 - w2*h2 - 3*w0*h0)/(2*dt)
    hu0 = sol0.q[1,:]
    div_huw = np.zeros_like(h0)
    div_huw[1:-1] = (hu0[2:]*w0[2:] - hu0[:-2]*w0[:-2])/(2*dx)
    p0 = (2./3.)*(d_hw_dt_2nd_order + div_huw)
    return p0

sol0 = Solution()
sol0.read(0, path='_output_WP_data', file_format='ascii', file_prefix='fort')
xc = mapc2p[0](sol0.grid.x.centers)
sol1 = Solution()
sol1.read(1, path='_output_WP_data', file_format='ascii', file_prefix='fort')
sol2 = Solution()
sol2.read(2, path='_output_WP_data', file_format='ascii', file_prefix='fort')
p = get_p(sol0, sol1, sol2)
w = get_w(sol0)
p_smooth = gaussian_filter1d(p, sigma=3)

#Save p0 and w0 to text files in the following format:
# 3         !mx
# 0.0 0.5 1.0 !x values
# 0.0 0.5 1.0 !w 
# 0.0 0.5 1.0 !p 

with open('w0_p0.txt', 'w') as f:
    f.write(f"{len(xc)}\n")
    f.write(" ".join(map(str, xc)) + "\n")
    f.write(" ".join(map(str, w)) + "\n")
    f.write(" ".join(map(str, p_smooth)) + "\n")