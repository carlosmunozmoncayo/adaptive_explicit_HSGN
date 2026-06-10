import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import tabulate
import os
from clawpack import pyclaw
from simulation_nontrivial_steady_SGN import setup

sech = lambda x: 1.0/np.cosh(x)
sinh = lambda x: np.sinh(x)
cosh = lambda x: np.cosh(x)
tanh = lambda x: np.tanh(x)

H0 = 1.0
a = 0.2
g = 9.81
csq = 1e5
tfinal = 0.1

def set_IC_steady_SGN_BBBD(x,h0=H0,a=a,g=g):
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


### Get fine grid solution
dom_len = 80.
mxfine = int(dom_len*2**20)
xfine = np.linspace(-dom_len/2,dom_len/2,mxfine)
hfine,ufine,wfine,sgfine,pfine,pbfine,bathfine = set_IC_steady_SGN_BBBD(x=xfine)
qfine = np.array([hfine,ufine,wfine,sgfine,pfine,pbfine,bathfine])


xupper = xfine[-1]
xlower = xfine[0]


### Get cell-averaged initial conditions
mx_list = [int(dom_len*2**i) for i in range(2,10)]
q_ICs = []
x_ICs = []
for mxcoarse in mx_list:
    ratio = int(mxfine/mxcoarse)
    print(f"mx= {mxcoarse} (2^{int(np.log2(mxcoarse/dom_len))}), fine grid points per coarse cell = {ratio}")
    qavg = np.zeros((7,mxcoarse))
    xavgs = np.zeros(mxcoarse)
    for j in range(mxcoarse):
        xavgs[j] = np.mean(xfine[j*ratio:(j+1)*ratio])
        qavg[:,j] = np.mean(qfine[:,j*ratio:(j+1)*ratio],axis=1)
    q_ICs.append(qavg)        
    x_ICs.append(xavgs)

### Run convergence tests Classic
for i in range(len(mx_list)):
    mx = mx_list[i]
    qinit = q_ICs[i]
    outdir = f"./_output_convergence_classic/mx_{mx}"
    claw = setup(mx=mx,qinit=qinit,outdir=outdir,c_sq=csq,solver_type='classic',
                xupper=xupper,xlower=xlower,tfinal=tfinal)
    claw.run()
    #Deallocate memory
    claw.solver.__del__()

### Run convergence tests SharpClaw
for i in range(len(mx_list)):
    mx = mx_list[i]
    qinit = q_ICs[i]
    outdir = f"./_output_convergence/mx_{mx}"
    claw = setup(mx=mx,qinit=qinit,outdir=outdir,c_sq=csq,
                xupper=xupper,xlower=xlower,tfinal=tfinal)
    claw.run()
    #Deallocate memory
    claw.solver.__del__()
