import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import tabulate
import os
from clawpack import pyclaw

sech = lambda x: 1.0/np.cosh(x)
sinh = lambda x: np.sinh(x)
cosh = lambda x: np.cosh(x)
tanh = lambda x: np.tanh(x)


get_n_frames = lambda path: len([name for name in os.listdir(path) if (os.path.isfile(os.path.join(path, name)) 
                                                                       and name.startswith('fort.q'))])

dom_len = 20.
dom_len = 4*dom_len
mx_list = [int(dom_len*2**i) for i in range(2,8)]#range(2,10)]


L1_distance = lambda q1, q2: np.sum(np.abs(q1-q2))

def errorfunc(sol,refsol,var=0):
    h = sol.q[var,:]
    h0 = refsol.q[var,:]
    dx = sol.grid.delta[0]
    return np.max(np.abs(h-h0))

def rate_conv(errold, errnew, dxold, dxnew):
    return (np.log(errnew)-np.log(errold))/(np.log(dxnew)-np.log(dxold))

def get_solutions(mx_list,path='_output_convergence/mx_'):
    path_head = path
    list_sols = []
    list_refsols = []
    for i in range(len(mx_list)):
        path = path_head + str(mx_list[i])+"/"
        n_frames = get_n_frames(path)
        sol = pyclaw.Solution(n_frames-1, file_format='ascii', file_prefix='fort', path=path)
        list_sols.append(sol)
        refsol = pyclaw.Solution(0, file_format='ascii', file_prefix='fort', path=path)
        list_refsols.append(refsol)
    return list_sols, list_refsols

def get_errors(mx_list,var=0,path='_output_convergence/mx_'):
    list_sols, list_refsols = get_solutions(mx_list, path=path)
    list_errors = []
    for i in range(len(list_sols)):
        list_errors.append(errorfunc(list_sols[i], list_refsols[i], var=var))
    return list_errors

def get_rate_conv(mx_list, var=0,path='_output_convergence/mx_'):
    list_sols, _ = get_solutions(mx_list, path=path)
    list_rate_conv = []
    error_list = get_errors(mx_list, var=var,path=path)
    for i in range(1, len(list_sols)):
        dxold = list_sols[i-1].grid.delta[0]
        dxnew = list_sols[i].grid.delta[0]
        list_rate_conv.append(rate_conv(error_list[i-1], error_list[i], dxold, dxnew))    
    return list_rate_conv


errors = get_errors(mx_list,path='_output_convergence_classic/mx_',var=0)
rates = get_rate_conv(mx_list,path='_output_convergence_classic/mx_',var=0)
#Make a table with errors and rates of convergence (scientific notation for errors)
table = []
#First row empty for rate of convergence
table.append([ mx_list[0], f"{errors[0]:.2e}", "-"])
for i in range(1,len(mx_list)):
    table.append([ mx_list[i], f"{errors[i]:.2e}", f"{rates[i-1]:.2f}"])
print(tabulate.tabulate(table, headers=['mx', 'Error', 'Rate'], tablefmt='github'))