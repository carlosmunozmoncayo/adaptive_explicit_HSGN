import numpy as np
from fortran_source_terms import source_terms_module
src_edc_1d = source_terms_module.src_edc_1d
src_bbbd_1d = source_terms_module.src_bbbd_1d
src_edc_1d_radial = source_terms_module.src_edc_1d_radial
src_bbbd_1d_radial = source_terms_module.src_bbbd_1d_radial
src_dr_1d = source_terms_module.src_dr_1d

def classic_source_step_EDC(solver,state,dt):
    q = state.q
    #Get parameters
    gamma = state.problem_data['gama']
    c_sq = state.problem_data['c_sq']
    csq_space = state.aux[1]
    # state.q = src_edc_1d(q=q ,dt=dt,gamma=gamma,csq=c_sq)
    temp = src_edc_1d(q=q ,dt=dt,gamma=gamma,csq=c_sq)
    state.q = np.where(csq_space>1., temp, state.q )

def classic_source_step_EDC_radial(solver,state,dt):
    q = state.q
    #Get parameters
    gamma = state.problem_data['gama']
    c_sq = state.problem_data['c_sq']
    r = np.abs(state.grid.x.centers)
    state.q = src_edc_1d_radial(q=q ,dt=dt,gamma=gamma,csq=c_sq,r=r)
    
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
    
    dq = np.zeros_like(q)
    dq[2,:] = np.where(csq_space>1.,dt*gamma*p-dt*decayx*q[2,:],0.)
    dq[3,:] = np.where(csq_space>1.,-2*dt*w*csq_space-dt*decayx*q[3,:],0.)

    return dq

def classic_source_step_BBBD(solver,state,dt):
    q = state.q
    #Get parameters
    c_sq = state.problem_data['c_sq']
    state.q = src_bbbd_1d(q=q ,dt=dt, csq=c_sq)

def classic_source_step_BBBD_radial(solver,state,dt):
    q = state.q
    #Get parameters
    c_sq = state.problem_data['c_sq']
    r = np.abs(state.grid.x.centers)
    state.q = src_bbbd_1d_radial(q=q ,dt=dt, csq=c_sq, r=r)

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
    # u = q[1,:]/h
    w = q[2,:]/h
    sg = q[3,:]/h
    p = q[4,:]/h
    pb = q[5,:]/h
    
    dq = np.zeros_like(q)
    dq[2,:] = dt*np.where(csq_space>-1.,pb,0.)
    dq[3,:] = dt*np.where(csq_space>-1.,-6*pb+12*p,0.)
    dq[4,:] = dt*np.where(csq_space>-1.,-csq_space*sg,0.)
    dq[5,:] = dt*np.where(csq_space>-1.,-6*csq_space*(w-0.5*sg),0.)

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
    dq[2,:] = np.where(csq_space>1.,dt*1.5*p-dt*decayx*q[2,:],0.)
    dq[3,:] = np.where(csq_space>1.,-2*dt*w*csq_space-dt*decayx*q[3,:],0.)

    return dq

if __name__ == '__main__':
    print("Not meant to be run as a script")