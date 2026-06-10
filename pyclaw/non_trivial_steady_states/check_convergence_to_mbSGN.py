import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

H0 = 1.0
q0 = 1.8
g = 9.81
gamma = 3./2.
alpha = 3.
c_sq = alpha**2 * g * H0

def get_steady_EDC(xc, H0=H0, q=q0, g=g, gamma=3./2., c_sq=c_sq,limit=False):
    def get_aux_f(a,b,c):
        f = lambda x: 2*c*(x - a)*np.exp(-b*(x - a)**2)
        fx = lambda x: 2*c*np.exp(-b*(a - x)**2)*(1 - 2*b*(a - x)**2)
        fxx = lambda x: -4*b*c*np.exp(-b*(a - x)**2)*(-3 + 2*b*(a - x)**2)*(a - x)
        return f, fx, fxx
    f,fx,fxx = get_aux_f(5.,3.4,1.5)

    def steady_EDC(x,q=q0,g=g,f=f,fx=fx,fxx=fxx,h0=H0,z0=0.,gamma=gamma,c_sq=c_sq):
        y0=np.array([h0,z0])
        if limit:
            def func(xlocal,y):
                h = y[0]
                z = y[1]
                rhs = np.zeros(2)
                w = f(xlocal)
                p = (q/gamma)*fx(xlocal)
                px = (q/gamma)*fxx(xlocal)
                rhs[0] = (-h*px-(g*h+gamma*p)*(w*h/q))/(0.5*g*h-q**2/(h**2)+(1-0.5*gamma)*p)
                rhs[1] = -0.5*rhs[0] + h*w/q
                return rhs
        else:
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
    return steady_EDC(xc)


if __name__ == "__main__":
    alpha_list = [2.,3.,5.,100.,1000000]
    c_sq_list = [(alpha*np.sqrt(g*H0))**2 for alpha in alpha_list]#100,200,300,400,500,1000,10000000]#[2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
    hlist = []
    hdotlist = []
    bathlist = []

    x= np.linspace(0, 10, 10000)
    q0= 10
    min_len = 10000000000
    for c_sq in c_sq_list:
        if c_sq > 100000:
            h,u,p,w,bath= get_steady_EDC(xc=x,q=q0,limit=True)
        else:
            h,u,p,w,bath = get_steady_EDC(xc=x,q=q0,c_sq=c_sq)
        hdot = np.gradient(h,x)
        hlist.append(h)
        hdotlist.append(hdot)
        bathlist.append(bath)

    Linf_error_list = []
    Linf_error_hdot_list = []
    L1_error_list = []
    dx = x[1] - x[0]
    for i in range(len(hlist)):
        Linf_error_list.append(np.linalg.norm(hlist[i]-hlist[-1], ord=np.inf))
        Linf_error_hdot_list.append(np.linalg.norm(hdotlist[i]-hdotlist[-1], ord=np.inf))
        L1_error_list.append(dx*np.linalg.norm(hlist[i]-hlist[-1], ord=1))

    #Plot both solutions and rate of convergence
    #Color list for left plot (increasingly lighter)
    color_list_left = ["#307adb", "#f87605", "#05750e", "#ff0000"]
    fig,ax = plt.subplots(1,2, figsize=(15,4))
    ax[0].set_title("Solutions of mbHSGN", fontsize=14)
    ax[0].set_xlabel(r"$x\; (m)$",fontsize=12)
    ax[0].set_ylabel(r"Water depth $h\;(m)$",fontsize=12)
    ax[1].set_title("Rate of convergence to mbSGN", fontsize=14)
    ax[1].set_xlabel(r"$c^2 = \alpha^2 g H\; (m^2/s^2)$",fontsize=12)
    ax[1].set_ylabel(r"$\|h - {h}_{mbSGN}\|$",fontsize=12)
    for i in range(len(hlist)):
        if alpha_list[i]> 100000:
            ax[0].plot(x,hlist[i],label=r"mbSGN", linestyle='--', color='black',lw=2)
        else:
            ax[0].plot(x,hlist[i],label=r"$\alpha=$"+f"{(alpha_list[i])}",lw=2, color=color_list_left[i])
    ax[1].loglog(c_sq_list[:-1], Linf_error_list[:-1], label=r"$L^\infty$ error h", marker='o',lw=2,color='blue')
    ax[1].loglog(c_sq_list[:-1], Linf_error_hdot_list[:-1], label=r"$L^\infty$ error $\partial_x h$", marker='x',markersize=8,lw=2,color='m')
    #Dashed grey line for linear convergence
    ax[1].loglog(c_sq_list[:-1], 0.5/np.array(c_sq_list[:-1]), linestyle='--', color='grey', label=r'$\mathcal{O}(1/c^2)$',lw=2)
    #Dashed grey line for quadratic
    ax[0].set_xlim(3.5,6.5)
    ax[0].legend(title=r"$c=\alpha \sqrt{g H}$",fontsize=12,title_fontsize=12,loc='upper right')
    #Force legend to right side
    ax[1].legend(fontsize=12,title_fontsize=12,loc='upper right')
    plt.tight_layout()
    #Save as pdf
    plt.savefig("compare_steady_states.pdf", bbox_inches='tight')