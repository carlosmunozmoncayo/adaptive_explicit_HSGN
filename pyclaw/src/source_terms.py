import numpy as np
from fortran_source_terms import source_terms_module
src_edc_1d = source_terms_module.src_edc_1d
src_bbbd_1d = source_terms_module.src_bbbd_1d
src_edc_1d_radial = source_terms_module.src_edc_1d_radial
src_bbbd_1d_radial = source_terms_module.src_bbbd_1d_radial
src_dr_1d = source_terms_module.src_dr_1d

# Bathymetry treatment for EDC (mbHSGN) and BBBD (HSGN):
# the Riemann solver only contains the hydrostatic bathymetry term g*hbar*db and one
# non-conservative product; every other bathymetry term is added here (b_x = cell_centred_bx).
# Classic solver (rs_bouss_hyprel_fwave_splitting_1D, non-conservative product -c^2*ubar*d(eta)):
#   EDC : -gamma p b_x (momentum), +c^2 u b_x (p-equation)
#   BBBD: -pb b_x (momentum), -c^2 u b_x (p-equation), +6 c^2 u b_x (pb-equation)
#   integrated implicitly together with the relaxation terms (one SDIRK(2,2) step of a local
#   linear system per cell, stable for any dt and any bathymetry jump).
# SharpClaw (rs_bouss_hyprel_sharpclaw_1D, non-conservative product -c^2*ubar*d(h), because SharpClaw
# reconstructs h): added explicitly to dq_src, like the relaxation terms there:
#   EDC : -gamma p b_x (momentum), +2 c^2 u b_x (p-equation)
#   BBBD: -pb b_x (momentum), +6 c^2 u b_x (pb-equation)


def cell_centred_bx(solver, state):
    """Cell-centred bottom slope (b_{i+1} - b_{i-1}) / (2 dx), with ghost values consistent
    with the aux BCs (periodic -> wrap; anything else -> constant extrapolation, which is also
    what a wall reflection of b gives). Serial only (uses state.aux without ghost cells)."""
    from clawpack import pyclaw
    b = state.aux[0]
    if solver.aux_bc_lower[0] == pyclaw.BC.periodic:
        bg = np.concatenate([b[-1:], b, b[:1]])
    else:
        bg = np.concatenate([b[:1], b, b[-1:]])
    return 0.5 * (bg[2:] - bg[:-2]) / state.grid.delta[0]

def classic_source_step_EDC(solver,state,dt):
    """SDIRK(2,2) on the local linear system for (u, w, p) that contains the relaxation
    terms AND the non-hydrostatic bathymetry terms -gamma p b_x, +c^2 u b_x."""
    q = state.q
    #Get parameters
    gamma = state.problem_data['gama']
    c_sq = state.problem_data['c_sq']
    csq_space = state.aux[1]
    bx = cell_centred_bx(solver, state)
    temp = src_edc_1d(q=q, bx=bx, dt=dt, gamma=gamma, csq=c_sq)
    state.q = np.where(csq_space>1., temp, state.q )

def classic_source_step_EDC_radial(solver,state,dt):
    q = state.q
    #Get parameters
    gamma = state.problem_data['gama']
    c_sq = state.problem_data['c_sq']
    r = np.abs(state.grid.x.centers)
    bx = cell_centred_bx(solver, state)
    state.q = src_edc_1d_radial(q=q, r=r, bx=bx, dt=dt, gamma=gamma, csq=c_sq)
    
def sharpclaw_source_step_EDC(solver,state,dt):
    q = state.q
    bath = state.aux[0]
    csq_space = state.aux[1]
    xc = state.grid.x.centers

    #Get parameters
    eta_sq = state.problem_data['eta_sq'] 
    g = state.problem_data['grav']
    gamma = state.problem_data['gama']
    c_sq = state.problem_data['c_sq']

    decayx = state.aux[2]

    #Get variables
    h = q[0,:]
    u = q[1,:]/h
    w = q[2,:]/h
    p = q[3,:]/h
    
    bx = cell_centred_bx(solver, state)

    dq = np.zeros_like(q)
    dq[1,:] = np.where(csq_space>1.e-6,-dt*gamma*p*bx,0.)
    dq[2,:] = np.where(csq_space>1.e-6,dt*gamma*p-dt*decayx*q[2,:],0.)
    dq[3,:] = np.where(csq_space>1.e-6,-2*dt*w*csq_space+2*dt*csq_space*u*bx-dt*decayx*q[3,:],0.)

    return dq

def classic_source_step_BBBD(solver,state,dt):
    """SDIRK(2,2) on the local linear system for (u, w, sigma, p, pb) that contains the relaxation
    terms AND the bathymetry terms that are not in the Riemann solver: -pb b_x (momentum),
    -c^2 u b_x (p-equation) and +6 c^2 u b_x (pb-equation).
    Constant c^2 = problem_data['c_sq'], applied in every cell (skipped for h < 1e-3)."""
    q = state.q
    #Get parameters
    c_sq = state.problem_data['c_sq']
    bx = cell_centred_bx(solver, state)
    state.q = src_bbbd_1d(q=q, bx=bx, dt=dt, csq=c_sq)

def classic_source_step_BBBD_radial(solver,state,dt):
    q = state.q
    #Get parameters
    c_sq = state.problem_data['c_sq']
    r = np.abs(state.grid.x.centers)
    bx = cell_centred_bx(solver, state)
    state.q = src_bbbd_1d_radial(q=q, r=r, bx=bx, dt=dt, csq=c_sq)

def sharpclaw_source_step_BBBD(solver,state,dt):
    q = state.q
    bath = state.aux[0]
    csq_space = state.aux[1]
    decayx = state.aux[2]
    xc = state.grid.x.centers

    #Get parameters
    eta_sq = state.problem_data['eta_sq'] 
    g = state.problem_data['grav']
    gamma = state.problem_data['gama']
    c_sq = state.problem_data['c_sq']

    #Get variables
    h = q[0,:]
    u = q[1,:]/h
    w = q[2,:]/h
    sg = q[3,:]/h
    p = q[4,:]/h
    pb = q[5,:]/h

    bx = cell_centred_bx(solver, state)
    
    dq = np.zeros_like(q)
    dq[1,:] = dt*np.where(csq_space>1.e-6,-pb*bx,0.)
    dq[2,:] = dt*np.where(csq_space>1.e-6,pb,0.)
    dq[3,:] = dt*np.where(csq_space>1.e-6,-6*pb+12*p,0.)
    dq[4,:] = dt*np.where(csq_space>1.e-6,-csq_space*sg,0.)
    dq[5,:] = dt*np.where(csq_space>1.e-6,-6*csq_space*(w-0.5*sg)+6*csq_space*u*bx,0.)

    return dq

def classic_source_step_DR(solver,state,dt):
    q = state.q
    #Get parameters
    c_sq = state.problem_data['c_sq']
    state.q = src_dr_1d(q=q ,dt=dt,csq=c_sq)

def sharpclaw_source_step_DR(solver,state,dt):
    q = state.q
    bath = state.aux[0]
    csq_space = state.aux[1]
    xc = state.grid.x.centers

    decayx = state.aux[2]

    #Get variables
    h = q[0,:]
    u = q[1,:]/h
    w = q[2,:]/h
    p = q[3,:]/h
    
    dq = np.zeros_like(q)
    dq[2,:] = np.where(csq_space>1.e-6,dt*1.5*p-dt*decayx*q[2,:],0.)
    dq[3,:] = np.where(csq_space>1.e-6,-2*dt*w*csq_space-dt*decayx*q[3,:],0.)

    return dq

if __name__ == '__main__':
    print("Not meant to be run as a script")