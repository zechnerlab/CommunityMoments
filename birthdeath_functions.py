import pandas as pd
import numpy as np
import scipy.stats as st
from scipy.integrate import solve_ivp
import tellurium as te

# Returns the gaussian approximation of the non-linear portion of the system.
def gaussian_closure(u, k_mean, k_std):
    b = k_mean[0]
    bb = k_std[0]**2 + k_mean[0]**2
    b0b1 = k_mean[0]*k_mean[0]
    bd = k_mean[0]*k_mean[1]
    d = k_mean[1]
    dd = k_std[1]**2 + k_mean[1]**2
    d0d1= k_mean[1]*k_mean[1]    
    return np.array([0,0,
                     -2*(d*u[2]+u[0]*u[6]+u[0]*u[6]-2*d*u[0]*u[0]),
                     -1*(d*u[3]+u[0]*u[11]+u[1]*u[6]-2*d*u[0]*u[1]),
                     -1*(d*u[4]+u[0]*u[8]+u[0]*u[6]-2*d*u[0]*u[0])-1*(d*u[4]+u[0]*u[6]+u[0]*u[8]-2*d*u[0]*u[0]),
                     -1*(d*u[5]+b*u[6]+u[0]*bd-2*d*b*u[0]),
                     -1*(d*u[6]+d*u[6]+u[0]*dd-2*d*d*u[0]),
                     -1*(d*u[7]+b*u[6]+u[0]*bd-2*d*b*u[0]),
                     -1*(d*u[6]+d*u[8]+u[0]*d0d1-2*d*d*u[0]),
                     0,0,0])

# Compiles the closed and approximated (non-closed) parts of the moment dynamics.
def approx_moments(u, t, A, a, k_mean, k_std):
    return np.dot(A, u) + a + gaussian_closure(u, k_mean, k_std)

# Returns the initial conditions of the system.
def initialize(i_mean, i_std, k_mean, k_std):
    return np.array([i_mean[0],
                     i_mean[1],
                     i_std[0]**2+i_mean[0]**2,
                     i_mean[0]*i_mean[1],
                     i_mean[0]*i_mean[0],
                     i_mean[0]*k_mean[0],
                     i_mean[0]*k_mean[1],
                     i_mean[0]*k_mean[0],
                     i_mean[0]*k_mean[1],
                     i_std[1]**2+i_mean[1]**2,
                     i_mean[1]*k_mean[0],
                     i_mean[1]*k_mean[1]])

# Returns the linear (closed) portion of the system.
def closed_moments(N, i_mean, i_std, k_mean, k_std):
    kout = k_mean[2]
    kin = k_mean[3]
    A = np.array([[-kout,kin,0,0,0,0,-1,0,0,0,0,0],
                  [N*kout,-N*kin,0,0,0,0,0,0,0,0,0,0],
                  [kout,kin,-2*kout,2*kin,0,2,1,0,0,0,0,0],
                  [-kout,-kin,kout,-N*kin-kout,(N-1)*kout,0,0,0,0,kin,1,0],
                  [0,0,0,2*kin,-2*kout,0,0,2,0,0,0,0],
                  [0,0,0,0,0,-kout,0,0,0,0,kin,0],
                  [0,0,0,0,0,0,-kout,0,0,0,0,kin],
                  [0,0,0,0,0,0,0,-kout,0,0,kin,0],
                  [0,0,0,0,0,0,0,0,-kout,0,0,kin],
                  [N*kout,N*kin,0,2*N*kout,0,0,0,0,0,-2*N*kin,0,0],
                  [0,0,0,0,0,kout,0,(N-1)*kout,0,0,-N*kin,0],
                  [0,0,0,0,0,0,kout,0,(N-1)*kout,0,0,-N*kin]])
    b = k_mean[0]
    bb = k_std[0]**2+k_mean[0]**2
    b0b1 = k_mean[0]*k_mean[0]
    bd = k_mean[0]*k_mean[1]
    d = k_mean[1]
    dd = k_std[1]**2+k_mean[1]**2
    d0d1= k_mean[1]*k_mean[1]
    a = np.array([b,0,b,0,0,bb,bd,b0b1,bd,0,0,0])
    return A, a

# Numerically solves the moment dynamics of the system.
def solve_moments(N, i_mean, i_std, k_mean, k_std):
    t_span = (0,1000)
    t_steps = np.linspace(0, 1000, 100)
    A, a = closed_moments(N, i_mean, i_std, k_mean, k_std)
    initial = initialize(i_mean, i_std, k_mean, k_std)
    u = solve_ivp(lambda t,y: approx_moments(y, t, A, a, k_mean, k_std), t_span, initial, t_eval = t_steps, method = 'Radau')
    return u

# Writes the Antimony file for the birth-death process.
def writeBirthDeath(initials, rates, N):
    ant = open('birthdeath.ant','w')
    ant.write("model birthdeath\n")
    rxn_id = 0
    for i in range(0,N):
        ant.write("R"+str(rxn_id+1)+": -> P"+str(i)+"; kb"+str(i)+" \n")
        ant.write("R"+str(rxn_id+2)+": P"+str(i)+" -> ; kd"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+3)+": P"+str(i)+" -> Ss; kout"+str(i)+"*P"+str(i)+" \n")
        ant.write("R"+str(rxn_id+4)+": Ss -> P"+str(i)+"; kin"+str(i)+"*Ss \n\n")
        rxn_id = rxn_id + 4

    for i in range(N):
        ant.write("kb"+str(i)+" = "+str(rates[i][0])+";\n")
        ant.write("kd"+str(i)+" = "+str(rates[i][1])+";\n")
        ant.write("kout"+str(i)+" = "+str(rates[i][2])+";\n")
        ant.write("kin"+str(i)+" = "+str(rates[i][3])+";\n")
        ant.write("P"+str(i)+" = "+str(initials[i][0])+";\n")
        
    ant.write("Ss = "+str(initials[0][1])+";\n")
    ant.write("end")
    ant.close()
    return

# Runs multiple stochastic simulations using Tellurium.
def runSSA(N, k, runs, std):
    # P, Ss.
    i_mean = np.array([20, 0]) 
    i_std = np.array([5, 0])
    # kb, kd, kout, kin.
    k_mean = np.array([1, 0.01, k, k])
    k_std = np.array([std, 0.001, 0, 0])
    s_stack = np.zeros(N+2)
    for i in range(runs):
        # Errors may be thrown out here if a negative value is drawn.
        initials = np.random.normal(i_mean,i_std, size = (N,len(i_mean)))
        rates = np.random.normal(k_mean,k_std, size = (N,len(k_mean)))
        # Rewrite Antimony file for N population cells.
        writeBirthDeath(initials.astype(int), rates, N)
        path = 'birthdeath.ant'
        r = te.loada(path)
        r.integrator = 'gillespie'
        selections = ['time'] + r.getFloatingSpeciesIds()
        r.integrator.variable_step_size = False
        r.resetToOrigin()
        s = r.simulate(0, 1000, 100, selections=selections)
        s_stack = np.vstack((s_stack, s))
    s_stack = s_stack[1:]
    s_stack_df = pd.DataFrame(s_stack, columns = s.colnames)
    s_stack_df['N'] = np.repeat(N,runs*100)
    s_stack_df['k'] = np.repeat(k,runs*100)
    s_stack_df['std'] = np.repeat(std,runs*100)
    return s_stack_df

# Bootstrapping to get confidence intervals of SSA runs.
def bootstrapping(data, typ, runs, pts, label):
    # typ bootstrap of mean, var, cv, or pv.
    # runs is for the number bootstraps.
    # pts is how many time points there are.
    # labels is the list of column names needed for the bootstrap calculation.
    raw_arr = np.array(data[label[0]]).reshape(runs,pts)
    idx = np.random.randint(runs, size = runs)
    boot_arr = raw_arr[idx,:]
    if typ == 'mean':
        boots = np.mean(boot_arr,axis = 0)
    elif typ == 'var':
        boots = np.std(boot_arr,axis = 0)**2 + np.mean(boot_arr,axis = 0)**2
    elif typ == 'cv':
        boots = np.std(boot_arr,axis = 0)/np.mean(boot_arr,axis = 0)
    elif typ == 'pv':
        d1 = np.array(data[label[1]]).reshape(runs,pts)[idx,:]
        d2 = np.array(data[label[2]]).reshape(runs,pts)[idx,:]
        boots = (np.mean(boot_arr,axis = 0)/(d1*d2))**0.5
    for i in range(runs-1):
        raw_arr = np.array(data[label[0]]).reshape(runs,pts)
        idx = np.random.randint(runs, size = runs)
        boot_arr = raw_arr[idx,:]
        if typ == 'mean':
            boots = np.vstack((boots, np.mean(boot_arr,axis = 0)))
        elif typ == 'var':
            boots = np.vstack((boots, np.std(boot_arr,axis = 0)**2 + np.mean(boot_arr,axis = 0)**2))
        elif typ == 'cv':
            boots = np.vstack((boots, np.std(boot_arr,axis = 0)/np.mean(boot_arr,axis = 0)))
        elif typ == 'pv':
            d1 = np.array(data[label[1]]).reshape(runs,pts)[idx,:]
            d2 = np.array(data[label[2]]).reshape(runs,pts)[idx,:]
            boots = np.vstack((boots, (np.mean(boot_arr,axis = 0)/(d1*d2))**0.5))
    ci = st.t.interval(0.95, len(boots), loc=np.mean(boots, axis = 0), scale=np.std(boots,axis = 0))
    ci_df = pd.DataFrame(np.array((ci[0],ci[1])).transpose())
    ci_df.columns = ['lower','upper']
    return ci_df