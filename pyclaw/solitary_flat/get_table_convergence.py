import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import tabulate
import os
from clawpack import pyclaw

### Self-convergence analysis

get_n_frames = lambda path: len([name for name in os.listdir(path) if (os.path.isfile(os.path.join(path, name)) 
                                                                       and name.startswith('fort.q'))])
mx_list = [2**i for i in range(6,11)]

mx_fine = 2**15 


def errorfunc(sol,refsol,var=0):
    h = sol[var,:]
    h0 = refsol[var,:]
    return np.max(np.abs(h-h0))

def rate_conv(errold, errnew, dxold, dxnew):
    return (np.log(errnew)-np.log(errold))/(np.log(dxnew)-np.log(dxold))

def get_solutions(mx_list,path='_output_self_convergence_classic_EDC/mx_'):
    path_head = path
    n_frames = get_n_frames(path_head + str(mx_fine) + "/")
    finesol = pyclaw.Solution(n_frames-1, file_format='ascii', file_prefix='fort', path=path_head + str(mx_fine) + "/").q
    list_sols = []
    dx_list = []
    list_refsols = []
    for i in range(len(mx_list)):
        path = path_head + str(mx_list[i])+"/"
        n_frames = get_n_frames(path)
        sol = pyclaw.Solution(n_frames-1, file_format='ascii', file_prefix='fort', path=path)
        dx_list.append(sol.grid.delta[0])
        sol = sol.q
        list_sols.append(sol)
        ratio = int(mx_fine/mx_list[i])
        refsol = np.zeros((sol.shape[0],mx_list[i]))
        for j in range(mx_list[i]):
            #Cell averaging fine solution
            refsol[:,j] = np.mean(finesol[:,j*ratio:(j+1)*ratio],axis=1)
        list_refsols.append(refsol)
    return list_sols, list_refsols, dx_list


def get_errors(mx_list,var=0,path='_output_self_convergence_classic_EDC/mx_'):
    list_sols, list_refsols, dx_list = get_solutions(mx_list, path=path)
    list_errors = []
    for i in range(len(list_sols)):
        list_errors.append(errorfunc(list_sols[i], list_refsols[i], var=var))
    return list_errors

def get_rate_conv(mx_list, var=0,path='_output_self_convergence_classic_EDC/mx_'):
    list_sols, _, dx_list = get_solutions(mx_list, path=path)
    list_rate_conv = []
    error_list = get_errors(mx_list, var=var, path=path)
    for i in range(1, len(list_sols)):
        dxold = dx_list[i-1]
        dxnew = dx_list[i]
        list_rate_conv.append(rate_conv(error_list[i-1], error_list[i], dxold, dxnew))    
    return list_rate_conv

if __name__ == "__main__":
    path_EDC = '_output_self_convergence_classic_EDC/mx_'
    path_BBBD = '_output_self_convergence_classic_BBBD/mx_'

    var = 0
    table = []
    table.append(["mx", "Error", "Rate", "mx", "Error", "Rate"])
    list_errors = []
    list_rates = []
    list_paths = ['_output_self_convergence_classic_EDC/mx_',
                '_output_self_convergence_classic_BBBD/mx_']
    for path in list_paths:
        errors = get_errors(mx_list,path=path,var=var)
        rates = get_rate_conv(mx_list,path=path,var=var)
        list_errors.append(errors)
        list_rates.append(rates)

    #First row empty for rate of convergence
    table.append([ mx_list[0], f"{list_errors[0][0]:.2e}", "-", mx_list[0], f"{list_errors[1][0]:.2e}", "-"])
    for i in range(1,len(mx_list)):
        table.append([ mx_list[i], f"{list_errors[0][i]:.2e}", f"{list_rates[0][i-1]:.2f}", mx_list[i], f"{list_errors[1][i]:.2e}", f"{list_rates[1][i-1]:.2f}"])

    table = tabulate.tabulate(table, headers=["", "mbHSGN","","", "HSGN",""], tablefmt='github')
    print(table)


