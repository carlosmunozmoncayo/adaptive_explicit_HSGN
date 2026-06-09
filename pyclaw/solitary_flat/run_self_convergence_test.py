import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from simulation_solitary_wave import setup


####################################
#EDC
####################################
g = 9.81
alpha = 3.
gamma = 1.5
H= 1.
c_sq = (alpha*np.sqrt(g*H))**2
A = 0.2
V = np.sqrt(g*(H+A))
l = np.sqrt(H*(A+H))


### Get fine grid solution
data_solitaryEDC = np.loadtxt('solitary_EDC_convergence.csv')
xfine,hfine,ufine,wfine,pfine, = data_solitaryEDC.T
qfine = np.array([xfine,hfine,ufine,wfine,pfine])
dom_len = xfine[-1] - xfine[0]
mxfine = len(xfine) #2**20


### Get cell-averaged initial conditions
mx_list = [2**i for i in range(6,12)] #For convergence tests
mx_list.append(2**15) #For reference solution
q_ICs = []
for mxcoarse in mx_list:
    ratio = int(mxfine/mxcoarse)
    print(f"mx= {mxcoarse} (2^{int(np.log2(mxcoarse))}), fine grid points per coarse cell = {ratio}")
    qavg = np.zeros((5,mxcoarse))
    xavgs = np.zeros(mxcoarse)
    for j in range(mxcoarse):
        qavg[:,j] = np.mean(qfine[:,j*ratio:(j+1)*ratio],axis=1)
    q_ICs.append(qavg)        

print(mx_list)

tfinal = 2.0

### Run convergence tests Classic
for i in range(len(mx_list)):
    mx = mx_list[i]
    qinit = q_ICs[i]
    outdir = f"./_output_self_convergence_classic_EDC/mx_{mx}"
    claw = setup(mx=mx,outdir=outdir,gamma=gamma,g=g,tfinal=tfinal,bouss_equations='EDC',
                 alpha=alpha,H=H,qinit=None,solver_type='classic')
    claw.run()
    #Deallocate memory
    claw.solver.__del__()


####################################
#BBBD
####################################
g = 9.81
alpha = 3.
gamma = 1.5
H= 1.
c_sq = (alpha*np.sqrt(g*H))**2
A = 0.2
V = np.sqrt(g*(H+A))
l = np.sqrt(H*(A+H))


### Get fine grid solution
data_solitaryBBBD = np.loadtxt('solitary_BBBD_convergence.csv')
xfine,hfine,ufine,wfine,sgfine,pfine,pbfine = data_solitaryBBBD.T
qfine = np.array([xfine,hfine,ufine,wfine,sgfine,pfine,pbfine])
dom_len = xfine[-1] - xfine[0]
mxfine = len(xfine)

q_ICs = []
for mxcoarse in mx_list:
    ratio = int(mxfine/mxcoarse)
    print(f"mx= {mxcoarse} (2^{int(np.log2(mxcoarse))}), fine grid points per coarse cell = {ratio}")
    qavg = np.zeros((7,mxcoarse))
    xavgs = np.zeros(mxcoarse)
    for j in range(mxcoarse):
        # xavgs[j] = np.mean(xfine[j*ratio:(j+1)*ratio])
        qavg[:,j] = np.mean(qfine[:,j*ratio:(j+1)*ratio],axis=1)
    q_ICs.append(qavg)        

print(mx_list)

tfinal = 2.0

### Run convergence tests Classic
for i in range(len(mx_list)):
    mx = mx_list[i]
    qinit = q_ICs[i]
    outdir = f"./_output_self_convergence_classic_BBBD/mx_{mx}"
    claw = setup(mx=mx,outdir=outdir,gamma=gamma,g=g,tfinal=tfinal,bouss_equations='BBBD',
                 alpha=alpha,H=H,qinit=None,solver_type='classic')
    claw.run()
    #Deallocate memory
    claw.solver.__del__()
