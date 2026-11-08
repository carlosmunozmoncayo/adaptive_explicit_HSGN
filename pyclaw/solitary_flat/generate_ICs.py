### Generate initial conditions for HSGN

import numpy as np
from scipy.integrate import solve_ivp

g=9.81
alpha=3.
H=1.
c_sq=(alpha*np.sqrt(g*H))**2
A=0.2
V=np.sqrt(g*(H+A))

mx = 2**20

def sourcetermBBBD(t,y):
    h=y[0]
    u=y[1]/y[0]
    w=y[2]/y[0]
    sigma=y[3]/y[0]
    p=y[4]/y[0]
    pb=y[5]/y[0]
    
    uvsq = (u - V)**2
    soundsq = c_sq + g*h + p 
    denom = (u - V)*(soundsq - uvsq)
    out = np.array([c_sq*sigma,
                    c_sq*V*sigma,
                    c_sq*w*sigma + (soundsq - uvsq)*pb,
                    12*p*(soundsq - uvsq) + c_sq*sigma**2 - 6*(soundsq - uvsq)*pb,
                    -c_sq*(g*h - uvsq)*sigma,
                    c_sq*(-3*(soundsq - uvsq)*(2*w - sigma) + sigma*pb)
                    ])
    return out/denom

q_eq = np.array([1., 0., 0., 0., 0., 0.])

#start along the unstable eigenvector
def numjac(f, q, h=1e-7):
    n = len(q); J = np.zeros((n, n))
    for i in range(n):
        d = np.zeros(n); d[i] = h
        J[:, i] = (f(0, q + d) - f(0, q - d)) / (2 * h)
    return J

w, Vmat = np.linalg.eig(numjac(sourcetermBBBD, q_eq))
i = np.argmax(w.real)  # the unstable eigenvalue
v = np.real(Vmat[:, i]); v /= np.linalg.norm(v)
if v[0] < 0: v = -v                   # orient so h rises

eps = 1e-12
q0 = q_eq + eps * v

# stop exactly at the crest:
def crest(t, y): return y[3]
crest.terminal = True
crest.direction = 0     # flip to +1/-1 if it triggers at t=0

sol = solve_ivp(sourcetermBBBD, [0, 500], q0, method='DOP853',
                rtol=1e-12, atol=1e-15, events=crest, dense_output=True)
tc = sol.t_events[0][0]

half_t = np.linspace(0, tc, mx // 2 + 1)
half_y = sol.sol(half_t)              
P = np.array([1, 1, -1, -1, 1, 1])[:, None]

left_t  = half_t - tc              
left_y  = half_y
right_t = (tc - half_t[::-1])[1:] 
right_y = (P * half_y[:, ::-1])[:, 1:]

full_t = np.concatenate([left_t, right_t])
full_y = np.concatenate([left_y, right_y], axis=1)

np.savetxt('solitary_BBBD_convergence.csv', np.transpose([full_t, full_y[0], full_y[1]/full_y[0], full_y[2]/full_y[0], full_y[3]/full_y[0], full_y[4]/full_y[0], full_y[5]/full_y[0]]))

### Generate initial conditions for mbHSGN
g, H, A = 9.81, 1.0, 0.2
alpha = 3         
c2 = (alpha*np.sqrt(g*H))**2
V  = np.sqrt(g*(H+A))
mx = 2**20

def hpprime(h): return -g*h+(H*V)**2/h**2-(0.5*g*H**2-0.5*g*h**2+(V*(h-H)/h)*H*V)/h

def pfunc(h): return g*(H**2-h**2)/(2*h) + H*V**2*(h-H)/h**2
def denom_func(h): return c2-hpprime(h)

def rhs(s, y):
    h, w = y
    return [-2*c2*w*h/(H*V*denom_func(h)), -1.5*pfunc(h)/(H*V)]

DH  = c2 - g*A 
assert DH > 0, "c^2 <= gA: equilibrium is a center — no solitary wave exists"
a   = -2*c2/(V*DH)
b   = -1.5*(g*A/H)/(H*V)
lam = np.sqrt(a*b)   # unstable eigenvalue > 0
v   = np.array([-a, -lam]); v /= np.linalg.norm(v)   # h rises (v[0]>0), w<0

eps = 1e-8
y0  = np.array([H, 0.0]) + eps*v

crest = lambda s, y: y[1]         # w = 0 at the crest
crest.terminal  = True
crest.direction = +1              # w dips negative, then climbs back through 0

sol = solve_ivp(rhs, [0, 40/lam], y0, method='DOP853',
                rtol=1e-12, atol=1e-15, events=crest, dense_output=True)
sc  = sol.t_events[0][0]          # "time" from tail to crest

# left half, crest centered at xi=0:  xi = s - sc  in [-sc, 0]
s_half = np.linspace(0, sc, mx//2 + 1)
y_half = sol.sol(s_half)
xi_left, h_left, w_left = s_half - sc, y_half[0], y_half[1]

# mirror: h even, w odd about the crest
xi = np.concatenate([xi_left, -xi_left[-2::-1]])
h  = np.concatenate([h_left,   h_left[-2::-1]])
w  = np.concatenate([w_left,  -w_left[-2::-1]])
p = pfunc(h)
u = V*(h-H)/h

np.savetxt('solitary_EDC_convergence.csv', np.transpose([xi, h, u, w, p]))