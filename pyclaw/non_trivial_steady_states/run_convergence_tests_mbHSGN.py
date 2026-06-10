import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from simulation_nontrivial_steady_mbHSGN import setup

H0 = 1.0
q0 = 1.8
g = 9.81
gamma = 1.5
c_sq = 50.

def get_aux_f(a,b,c):
    f = lambda x: 2*c*(x - a)*np.exp(-b*(x - a)**2)
    fx = lambda x: 2*c*np.exp(-b*(a - x)**2)*(1 - 2*b*(a - x)**2)
    fxx = lambda x: -4*b*c*np.exp(-b*(a - x)**2)*(-3 + 2*b*(a - x)**2)*(a - x)
    return f, fx, fxx
f,fx,fxx = get_aux_f(5.,3.4,1.5)

def steady_EDC(x,q=q0,g=g,f=f,fx=fx,fxx=fxx,h0=H0,z0=0.,gamma=gamma,c_sq=c_sq):
    y0=np.array([h0,z0])
    def func(xlocal,y):
        h = y[0]
        z = y[1]
        rhs = np.zeros(2)
        w = f(xlocal)
        p = (q/gamma)*fx(xlocal)
        px = (q/gamma)*fxx(xlocal)
        rhs[0] = (-h*px-(g*h+gamma*p)*(h*0.5*px/c_sq+w*h/q))/(0.5*g*h-q**2/(h**2)+(1-0.5*gamma)*p)
        rhs[1] = -0.5*rhs[0] + h*w/q + h*0.5*px/c_sq
        return rhs
    x_span = np.array([x[0],x[-1]])
    sol = solve_ivp(func,x_span,y0,method='DOP853',t_eval=x,rtol=1.e-14,atol=1.e-14)
    h = sol.y[0]
    bath = sol.y[1]
    w = f(x)
    p = (q/gamma)*fx(x)
    u = q/h
    print(sol.status)
    return h,u,p,w,bath


### Get fine grid solution
dom_len = 40.
mxfine = int(dom_len*2**20)
xfine = np.linspace(-dom_len/2+5.,dom_len/2+5.,mxfine)
hfine,ufine,pfine,wfine,bathfine = steady_EDC(xfine)
plt.plot(hfine)
plt.show()
for x in [hfine, ufine, wfine, pfine, bathfine]:
    print(x.shape)
qfine = np.array([hfine,ufine,wfine,pfine,bathfine])

xupper = xfine[-1]
xlower = xfine[0]

### Get cell-averaged initial conditions
mx_list = [int(dom_len*2**i) for i in range(5,12)]
q_ICs = []
x_ICs = []
for mxcoarse in mx_list:
    ratio = int(mxfine/mxcoarse)
    print(f"mx= {mxcoarse} (2^{int(np.log2(mxcoarse/dom_len))}), fine grid points per coarse cell = {ratio}")
    qavg = np.zeros((5,mxcoarse))
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
    claw = setup(mx=mx,qinit=qinit,outdir=outdir,gamma=gamma,
                 xupper=xupper,xlower=xlower,c_sq=c_sq,solver_type='classic')
    claw.run()
    #Deallocate memory
    claw.solver.__del__()