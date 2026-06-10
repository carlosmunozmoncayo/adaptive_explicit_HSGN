#!/usr/bin/env python
# encoding: utf-8

r"""
Dispersive water waves modelling
This non-hydrostatic hyperbolic system of eqs. was proposed by Escalante et al. (2019)
==================

Solve the one-dimensional equations:

.. math::
    q_t+F(q)_x+B(q)q_x=S(q),\\
    q=(h, hu, hw, hp, h0).

Here h is the depth, u is the horizontal velocity, p is the non hydrostatic pressure,
w is the vertical velocity, and h0(x) is the given depth of the ocean at rest.
This script simulates the propagation of a solitary wave solution for the hyperbolic
relaxation system over a frictionless bottom.
"""

from __future__ import absolute_import
import sys
import os
import numpy as np
from clawpack import pyclaw

#Importing Riemann solver and source terms
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.path.join(parent_dir, 'src'))
import rs_bouss_hyprel_fwave_splitting_1D as rs
from dispersive_tfluct import tfluct1


MIN_h = 1E10
MAX_u = -1E10 
MIN_dt = 1E10

# h=0
# hu=1
# hw=2
# hp=3
sech = lambda x: 1.0/np.cosh(x)
sinh = lambda x: np.sinh(x)
cosh = lambda x: np.cosh(x)
tanh = lambda x: np.tanh(x)

def setup(outdir='./_output',solver_type='sharpclaw',
          bouss_equations='BBBD',
          mx=4000, gamma=3./2., order=2,
          transition="no", filter_type="BTE",
          decayparam=(2,20),transition_after_time=-1,
          xlower=-40., xupper=40.,
          tfinal=1.,nout=50, relaxation_zone=(100,130),
          petsc=False,qinit=None, c_sq=1000,g=9.81):

    if petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if bouss_equations == 'EDC':
        print("Setting up neq, nwaves, and sources EDC system")
        num_eqn=4
        num_waves=4
        from source_terms import classic_source_step_EDC as classic_source_step
        from source_terms import sharpclaw_source_step_EDC as sharpclaw_source_step

    elif bouss_equations == 'BBBD':
        print("Setting up neq, nwaves, and sources BBBD system")
        num_eqn=6
        num_waves=6
        from source_terms import classic_source_step_BBBD as classic_source_step
        from source_terms import sharpclaw_source_step_BBBD as sharpclaw_source_step
        
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.order = order
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.fwave=True
        solver.step_source=classic_source_step
        solver.source_split = 2
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)
        solver.fwave=True
        solver.weno_order=5
        solver.lim_type=2
        solver.tfluct = tfluct1
        solver.dq_src = sharpclaw_source_step
        solver.call_before_step_each_stage = False

    solver.num_waves = num_waves
    solver.num_eqn = num_eqn
    
    solver.dt_initial = 1e-6
    solver.cfl_desired = 0.9
    solver.cfl_max = 1.0

    solver.max_steps=int(1e8)

    #Define spatial domain
    if qinit is not None:
        mx = qinit[0].shape[-1]
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    xc = x.centers
    dx = x.delta
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn,num_aux=3)
    #First compoment of aux will correspond to bathymetry
    #Second component will correspond to parameter c^2 (changing in space)
    #Third component will correspond to a decay parameter for damping
    #of hw and hp (EDC system)

    a = 0.2
    h0 = 1.0
    L = xupper-xlower #Physical domain length
    alpha = np.sqrt(mx/L/h0/g)
    h0=1.
    # c_sq=1000#100.*(alpha*np.sqrt(g*h0))**2 #Such that c^2~2/dx and we are in the resolved regime
    print(f"c_sq={c_sq}")
    print(f"1/dx={1./dx}")

    state.problem_data['grav'] = g
    state.problem_data['gama'] = gamma
    state.problem_data['eta_sq'] = 1.0
    state.problem_data['c_sq'] = c_sq
    state.problem_data['dry_tolerance'] = 1e-4
    state.problem_data['sea_level'] = 1.
    if bouss_equations == 'EDC':
        state.problem_data['bouss_equations'] = 1
    elif bouss_equations == 'BBBD':
        state.problem_data['bouss_equations'] = 2
    else:
        state.problem_data['bouss_equations'] = 0
    
    xc = state.grid.x.centers

    #Auxiliary vector will contain bathymetry
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap


    #Setting the transition method from the hyperbolic relaxation BTE
    #to shallow water equations
    filx, decayx= transition_HypRelBTE_SWE(transition=transition,filter_type=filter_type,
                                           relaxation_zone=relaxation_zone,xc=xc,decayparam=decayparam)
                                               

    #Setting initial conditions
    if qinit is not None:
        print("Setting initial conditions from input")
        h,u,w,sg,p,pb,z= qinit
    else:
        print("Generating initial conditions")
        h,u,w,sg,p,pb,z= set_IC_steady_SGN_BBBD(x=xc,h0=1,a=0.2,g=9.81)
    state.q[0,:] = h
    state.q[1,:] = h*u
    if bouss_equations == 'EDC':
        # print("Setting initial conditions for EDC system")
        # state.q[2,:] = 0.0
        # state.q[3,:] = 0.0
        ValueError("Only BBBD system is implemented")     
    if bouss_equations == 'BBBD':
        # print("Setting initial conditions for BBBD system")
        state.q[2,:] = h*w
        state.q[3,:] = h*sg
        state.q[4,:] = h*p
        state.q[5,:] = h*pb
        ICs_for_boundary = np.copy(state.q)

    state.aux[0,:] = z
    state.aux[1,:] = c_sq
    state.aux[2,:] = 0.#decayx


    if transition in ["filter", "no"]: 
        def b4step(self, state, filx=filx):
            #Get variables
            csq_space = state.aux[1]
            xc = state.grid.x.centers
            q = state.q
            #Damp non-hydrostatic variables
            if transition_after_time < state.t:
                q[2:,:] = q[2:,:]*filx
                state.aux[1,:] = c_sq*filx
            #Relaxation absorbing boundary conditions
            xc = state.grid.x.centers
            L = xc[-1] - xc[0]
            upperx = xc[-1] - 0.1*L
            lowerx = xc[0] + 0.1*L
            damp_gamma = np.ones_like(xc)
            damp_gamma = np.where(xc>upperx,(xc[-1]-xc)/(xc[-1]-upperx),damp_gamma)
            damp_gamma = np.where(xc<lowerx,(xc-xc[0])/(lowerx-xc[0]),damp_gamma)
            for i in range(state.q.shape[0]):
                state.q[i,:] = state.q[i,:]*np.sqrt(damp_gamma)+(1.-np.sqrt(damp_gamma))*ICs_for_boundary[i,:]
        # solver.before_step = b4step

    #solver.call_before_step_each_stage=True
    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = tfinal
    claw.num_output_times = nout
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.write_aux_init = True
    return claw

def set_IC_steady_SGN_BBBD(x,h0=1,a=0.2,g=9.81):
    r = (1.0/h0)*np.sqrt(3.0*a/(1.0+a))
    q = (1.0/h0)*np.sqrt(0.5*(1.0+a)*g*h0**3.)
    h = h0*(1.0+a/(cosh(r*x))**2.0)
    u = q/h
    w = (-2*a*(-1 + h0)*r*q*tanh(r*x))/(h0*(1 + 2*a + cosh(2*r*x)))
    sg = (-4*a*r*q*tanh(r*x))/(1 + 2*a + cosh(2*r*x))
    p = (a*(-3 + 4*h0)*r**2*(-3 - 4*a - 2*cosh(2*r*x) + cosh(4*r*x))*q**2*
        sech(r*x)**2)/(6.*h0*(1 + 2*a + cosh(2*r*x))**2)
    pb = (a*(-1 + h0)*r**2*(-3 - 4*a - 2*cosh(2*r*x) + cosh(4*r*x))*q**2*
        sech(r*x)**2)/(h0*(1 + 2*a + cosh(2*r*x))**2)
    z = -0.5*a/(cosh(r*x))**2-h0
    return h,u,w,sg,p,pb,z
    
def transition_HypRelBTE_SWE(transition,filter_type,xc,relaxation_zone, decayparam):
    
    def bump_filter(x,a,b):
        #Generates a smooth function =1 to the left of a and =0 to the right of b
        #Assumes b>a>=0
        f = lambda x: np.where(x>1.e-16, np.exp(-1/x), 0.)
        g = lambda x: f(x)/(f(x)+f(1-x))
        return -g((x-a)/(b-a))+1

    def affine_filter(x,a,b):
        #Defining relaxation zone
        #Generates an affine function equal to 1 for x<=a and equal to 0 for x>b
        Lrelax = b-a
        di = x-b #distance of points from relaxation zone's origin 
        mi = np.abs(di)/Lrelax #slope of filter in relaxation zone
        return (mi*(mi<=1)+(mi>1))*(di<=0)

    def sigma_decay_rate(x,x0,D,a,b):
        return (x>=x0)*(x<x0+D)*b*np.abs((x-x0)/D)**a+(x>=x0+D)*b

    filx=np.ones_like(xc)
    decayx=np.zeros_like(xc)

    if transition in ["filter", "no"]:
        #Filter type
        if filter_type == "jump": 
            filter_type = lambda x,a,b: (x<=a)
        elif filter_type == "affine":
            filter_type = affine_filter
        elif filter_type == "bump":
            filter_type = bump_filter
        elif filter_type == "BTE":
            filter_type = lambda x,a,b: np.ones_like(x)
        elif filter_type == "SWE":
            filter_type = lambda x,a,b: np.zeros_like(x)
        else:
            print(f"ERROR: No available filter: {filter_type}")
            filter_type = None 

        #Computing filter over domain
        filx = filter_type(x=xc,a=relaxation_zone[0],b=relaxation_zone[1])

    elif transition == "decay":
        #Computing decay rate
        #Decay rate starts increasing at x=a
        (a,b)=decayparam
        D = relaxation_zone[1]-relaxation_zone[0]
        decayx = sigma_decay_rate(x=xc,x0=relaxation_zone[0],D=D,a=a,b=b)

    else: print(f"ERROR: No available transition: {transition}")

    return filx, decayx

    
#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """
    #xc=plotdata.frames[0].grid.x.centers
    #Spatial limits
    #xmin=xc[0]#-20
    #xmax=xc[-1]#300
    # Plot variables
    def bathy(current_data):
        return current_data.aux[0, :]

    def eta(current_data):
        return current_data.q[0, :] + bathy(current_data)
    
    rgb_converter = lambda triple: [float(rgb) / 255.0 for rgb in triple]

    #plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = [xmin,xmax]
    # plotaxes.ylimits = [0.8,1.25]
    plotaxes.title = 'h'
    plotaxes.axescmd = 'subplot(221)'
    

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta
    plotitem.plot_var2= bathy
    plotitem.color = rgb_converter((67,183,219))
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bathy
    plotitem.color = 'k'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta
    plotitem.color = 'k'
    plotitem.kwargs = {'linewidth':1}

    # Figure for momentum[1]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(222)'
    # plotaxes.xlimits = [xmin,xmax]
    # plotaxes.ylimits = [-0.1,0.7]
    plotaxes.title = 'u-Momentum'
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':1}

    # Figure for w tracer[2]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(223)'
    # plotaxes.xlimits = [xmin,xmax]
    plotaxes.ylimits = [-0.15,0.15]
    plotaxes.title = 'hw'
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 2
    plotitem.plotstyle = '-'
    plotitem.color = 'g'
    plotitem.kwargs = {'linewidth':1}


    # Figure for p tracer[2]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(224)'
    # plotaxes.xlimits = [xmin,xmax]
    plotaxes.ylimits = [-1.e-12,1.e-12]
    plotaxes.title = 'hp'
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 3
    plotitem.plotstyle = '-'
    plotitem.color = 'r'
    plotitem.kwargs = {'linewidth':1}

    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)

  

