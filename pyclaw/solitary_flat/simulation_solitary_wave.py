#!/usr/bin/env python
# encoding: utf-8

r"""
Solitary wave test over a flat topography.
"""

from __future__ import absolute_import
import sys
import os
import numpy as np

#Importing Riemann solver and source terms
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.path.join(parent_dir, 'src'))
import rs_bouss_hyprel_fwave_splitting_1D as rs
from dispersive_tfluct import tfluct1
import bouss_hyprel_sharpclaw1


MIN_h = 1E10
MAX_u = -1E10 
MIN_dt = 1E10

h=0
hu=1
hw=2
hp=3

def setup(outdir='./_output',solver_type='sharpclaw',riemann_solver='dispersive',
          bouss_equations='EDC',
          mx=2**12, gamma=1.5,  bath='flat', order=2,
          transition="no", filter_type="BTE",
          decayparam=(2,20),transition_after_time=-1,
          IC='solitary',
          tfinal=3.,nout=10, relaxation_zone=(100,130),
          g = 9.81, alpha=3.0, H=1.,
          qinit=None,use_char_decomp=False,waveslope=False, use_petsc=False):
    
    if use_petsc:
        import clawpack.petclaw as pyclaw
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

    elif bouss_equations == 'DR':
        print("Setting up neq, nwaves, and sources DR system")
        num_eqn=5
        num_waves=5
        from source_terms import classic_source_step_DR as classic_source_step
        from source_terms import sharpclaw_source_step_DR as sharpclaw_source_step

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.order = order
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.source_split = 2
        solver.fwave=True
        solver.step_source=classic_source_step
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)
        solver.fwave=True
        solver.tfluct = tfluct1
        solver.weno_order=5
        solver.lim_type=2
        if use_char_decomp:
            solver.fmod = bouss_hyprel_sharpclaw1
            solver.lim_type = 2             # WENO reconstruction
            solver.char_decomp = 2          # characteristic-wise reconstruction
        elif waveslope:
            solver.char_decomp=1 #wave-slope reconstruction
        
        solver.dq_src = sharpclaw_source_step
        solver.call_before_step_each_stage = True

    solver.num_waves = num_waves
    solver.num_eqn = num_eqn
    
    solver.dt_initial = 1e-6
    solver.cfl_desired = 0.9
    solver.cfl_max = 1.0

    #############################
    #Initial conditions
    #Check if qinit is an array
    if qinit is None:
        #Load solitary data
        if bouss_equations == 'EDC':
            data_solitary = np.loadtxt('solitary_EDC_convergence.csv')
            #data_solitary = data_solitary[::4]
            xinit,h,u,w,p = data_solitary.T
            qinit = np.array([xinit,h,u,w,p])

        elif bouss_equations == 'BBBD':
            data_solitary = np.loadtxt('solitary_BBBD_convergence.csv')
            # data_solitary = data_solitary[::12]
            xinit,h,u,w,sg,p,pb = data_solitary.T
            qinit = np.array([xinit,h,u,w,sg,p,pb])

        else:
            raise ValueError(f"Unknown bouss_equations: {bouss_equations}. ")
    else:
        if bouss_equations == 'EDC':
            xinit,h,u,w,p = qinit
        elif bouss_equations == 'BBBD':
            xinit,h,u,w,sg,p,pb = qinit
        else:
            raise ValueError(f"Unknown bouss_equations: {bouss_equations}. ")
    #############################

    #Undersample to match mx if necessary
    if len(xinit) != mx:
        mxfine = len(xinit)
        ratio = int(mxfine/mx)
        print(f"mx= {mx} (2^{int(np.log2(mx))}), fine grid points per coarse cell = {ratio}")
        qavg = np.zeros((qinit.shape[0],mx))
        for j in range(mx):
            qavg[:,j] = np.mean(qinit[:,j*ratio:(j+1)*ratio],axis=1)
        xinit = qavg[0,:]
        h = qavg[1,:]
        u = qavg[2,:]
        if bouss_equations == 'EDC':
            w = qavg[3,:]
            p = qavg[4,:]
        elif bouss_equations == 'BBBD':
            w = qavg[3,:]
            sg = qavg[4,:]
            p = qavg[5,:]
            pb = qavg[6,:]
        else:
            print("Not handling DR system yet")
    else:
        print(f"mx= {mx} (2^{int(np.log2(mx))}), no undersampling needed")
        qavg[:,:] = qinit[:,:]
    
    dx = xinit[1]-xinit[0]
    print(f"dx = {dx}")
    xlower,xupper = xinit[0]-0.5*dx,xinit[-1]+0.5*dx
    mx = len(xinit)
    ##############################

    #Define spatial domain
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    xc = x.centers
    #Save xc to file
    np.savetxt('xc.csv',xc)
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn,num_aux=3)
    #First compoment of aux will correspond to bathymetry
    #Second component will correspond to parameter c^2 (changing in space)
    #Third component will correspond to a decay parameter for damping
    #of hw and hp (EDC system)

    #Problem parameters
    c_sq=(alpha*np.sqrt(g*H))**2
    state.problem_data['grav'] = g
    state.problem_data['gama'] = gamma
    state.problem_data['eta_sq'] = 1.0
    state.problem_data['c_sq'] = c_sq
    state.problem_data['dry_tolerance'] = 1e-4
    state.problem_data['sea_level'] = 0.
    state.problem_data["alphaldr"] = 1.153
    if bouss_equations == 'EDC':
        state.problem_data['bouss_equations'] = 1
    elif bouss_equations == 'BBBD':
        state.problem_data['bouss_equations'] = 2
    elif bouss_equations == 'DR':
        state.problem_data['bouss_equations'] = 3
    else:
        state.problem_data['bouss_equations'] = 0
    
    xc = state.grid.x.centers

    #Auxiliary vector will contain bathymetry
    solver.aux_bc_lower[0] = pyclaw.BC.periodic
    solver.aux_bc_upper[0] = pyclaw.BC.periodic

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic


    #Setting the transition method from the hyperbolic relaxation BTE
    #to shallow water equations
    filx, decayx= transition_HypRelBTE_SWE(transition=transition,filter_type=filter_type,
                                           relaxation_zone=relaxation_zone,xc=xc,decayparam=decayparam)
                                           

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
        solver.before_step = b4step
    

    #Setting initial conditions
    # state.q[:,:] = 0.
    if use_petsc:
        bound1 = np.where(xinit>=xc[0])[0][0]
        bound2 = np.where(xinit<=xc[-1])[0][-1]+1
        state.q[0,:] = h[bound1:bound2]
        state.q[1,:] = state.q[0,:]*u[bound1:bound2]
        if bouss_equations == 'EDC':
            print("Setting initial conditions for EDC system")
            state.q[2,:] = state.q[0,:]*w[bound1:bound2]
            state.q[3,:] = state.q[0,:]*p[bound1:bound2]
        if bouss_equations == 'BBBD':
            print("Setting initial conditions for BBBD system")
            state.q[2,:] = state.q[0,:]*w[bound1:bound2]
            state.q[3,:] = state.q[0,:]*sg[bound1:bound2]
            state.q[4,:] = state.q[0,:]*p[bound1:bound2]
            state.q[5,:] = state.q[0,:]*pb[bound1:bound2]
        state.aux[0,:] = 0.
        state.aux[1,:] = c_sq
        state.aux[2,:] = 0.#decayx
    else:
        state.q[0,:] = h
        state.q[1,:] = h*u
        if bouss_equations == 'EDC':
            print("Setting initial conditions for EDC system")
            state.q[2,:] = h*w
            state.q[3,:] = h*p
        if bouss_equations == 'BBBD':
            print("Setting initial conditions for BBBD system")
            state.q[2,:] = h*w
            state.q[3,:] = h*sg
            state.q[4,:] = h*p
            state.q[5,:] = h*pb
        state.aux[0,:] = 0.
        state.aux[1,:] = c_sq
        state.aux[2,:] = 0.#decayx

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
    plotaxes.ylimits = [0.8,1.25]
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
    plotaxes.ylimits = [-0.25,0.25]
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

  

