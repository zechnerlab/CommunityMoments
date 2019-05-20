import numpy as np
import pandas as pd
import tellurium as te
import itertools
import sympy
from sympy import symbols
from sympy.parsing.sympy_parser import parse_expr
from scipy.integrate import odeint

def getRateLaws(path):
    ant_file = open(path).readlines()
    rate_laws = []
    for line in ant_file:
        if line[0] == 'R':
            rate_laws.append(line.split(';')[1].strip())
    return rate_laws

def getReactionInfo(path):
    r = te.loada(path)
    keys = r.getFloatingSpeciesIds() + r.getGlobalParameterIds()
    rxn_stoich = np.zeros(shape = [len(r.getGlobalParameterIds()), r.getNumReactions()])
    net_stoich = np.concatenate((r.getFullStoichiometryMatrix(), rxn_stoich), axis = 0)
    return keys, net_stoich

def getMoments(net_stoich, rate_laws, keys):
    (M, R) = np.shape(net_stoich) # M = species + rate constant, R = reactions
    order = 2 # moments up to 2nd order only
    odes = [list(itertools.combinations_with_replacement(range(1,M+1),i)) for i in range(1,order+1)]
    odes_list = list(itertools.chain.from_iterable(odes))
    # collect moment terms
    moment_idx = []
    for item in odes_list:
        product = 1
        for factor in item:
            product = sympy.simplify(product * parse_expr(keys[int(factor)-1]))
        moment_idx.append(product)
    # collect moment equations
    moments = np.array([])
    for ode in odes_list:
        powers = [ode.count(i+1) for i in range(M)]
        summ = 0
        for rxn in range(R):
            rate_law = parse_expr(rate_laws[rxn])
            prod1 = 1
            prod2 = 1
            for idx,x in enumerate(keys):
                y = parse_expr(x)
                prod1 = prod1 * (y+net_stoich[:,rxn][idx])**powers[idx]
                prod2 = prod2 * y**powers[idx]
            summ = summ + (rate_law * (prod1 - prod2))
        summ = sympy.simplify(summ)
        summ = summ.expand(basic = True)
        moments = np.append(moments, summ)
    return M, moments, moment_idx

def getCommunityEffect(N, moments, moment_idx, ss, kout0, kin0, kout1, kin1):
    moment_idx_str = [str(i) for i in moment_idx]
    pop_moments = np.array([])
    for idx, m in enumerate(moments):
        replace = 0
        if moment_idx[idx].has(ss):
            if '1' in moment_idx_str[idx]:
                to_rep = [kout0, kin0]
            else:
                to_rep = [kout1, kin1]
            for arg in m.args:
                if arg.has(to_rep[0]) or arg.has(to_rep[1]):
                    replace = replace + arg*(N-1)
                else:
                    replace = replace + arg
            pop_moments = np.append(pop_moments, replace)
        else:
            pop_moments = np.append(pop_moments, m)
    return pop_moments   

def getMatrixForm(equation, idx, odes_rev, odes_str):
    moments_row = pd.DataFrame([np.repeat(0,len(odes_rev))], columns=odes_str)
    moments_row.index = [odes_str[::-1][idx]]
    for expr in equation.args:
        for i in range(len(odes_rev)):
            coeff = expr.coeff(odes_rev[i])
            if coeff != 0:
                moments_row[[odes_str[i]]] = moments_row[[odes_str[i]]] + coeff
                break
    return moments_row

def getMomentMatrix(M, pop_moments, keys):
    order = 3 # close moments at 3rd order 
    odes = [list(itertools.combinations_with_replacement(range(1,M+1),i)) for i in range(1,order+1)]
    odes_list = list(itertools.chain.from_iterable(odes))
    odes_idx = []
    for item in odes_list:
        product = 1
        for factor in item:
            product = sympy.simplify(product * parse_expr(keys[int(factor)-1]))
        odes_idx.append(product)
    odes_rev = odes_idx[::-1]
    odes_str = [str(i) for i in odes_rev]
    moment_mat = pd.DataFrame(columns=odes_str)
    for idx, equation in enumerate(pop_moments):
        moment_mat = moment_mat.append(getMatrixForm(equation, idx, odes_rev, odes_str))
    return moment_mat

def getClosureKey(B_cols):
    cols = np.array([parse_expr(i) for i in B_cols])
    mu_bk = []
    for i in cols:
        f = []
        for j in sympy.factor_list(i)[1]:
            if j[1] == 1:
                f.append(j[0])
            elif j[1] == 2:
                f.append(j[0])
                f.append(j[0])
            elif j[1] == 3:
                f.append(j[0])
                f.append(j[0])
                f.append(j[0])
            elif j[1] == 4:
                f.append(j[0])
                f.append(j[0])
                f.append(j[0])
                f.append(j[0])
        if len(f) == 3:
            mu_bk.append([str(sympy.simplify(f[0]*f[1])),
                         str(sympy.simplify(f[0]*f[2])),
                         str(sympy.simplify(f[1]*f[2])),
                         str(sympy.simplify(f[0])),
                         str(sympy.simplify(f[1])),
                         str(sympy.simplify(f[2]))])
        elif len(f) == 4:
            mu_bk.append([str(sympy.simplify(f[0]*f[1]*f[2])),
                         str(sympy.simplify(f[0]*f[1]*f[3])),
                         str(sympy.simplify(f[1]*f[2]*f[3])),
                         str(sympy.simplify(f[0]*f[2]*f[3])),
                         str(sympy.simplify(f[0]*f[1])),
                         str(sympy.simplify(f[0]*f[2])),
                         str(sympy.simplify(f[0]*f[3])),
                         str(sympy.simplify(f[1]*f[2])),
                         str(sympy.simplify(f[1]*f[3])),
                         str(sympy.simplify(f[2]*f[3])),
                         str(sympy.simplify(f[0])),
                         str(sympy.simplify(f[1])),
                         str(sympy.simplify(f[2])),
                         str(sympy.simplify(f[3]))])
    return mu_bk

def getInitialConditions(mu, mu_dict):
    mu_0 = np.zeros(len(mu))
    for idx,i in enumerate(mu):
        factors = sympy.factor_list(i)
        mu_c = 1
        for j in factors[1]:
            mu_c = mu_c * mu_dict[str(j)]
        mu_0[idx] = mu_c
    return mu_0

def getLognormalClosure(u, mu_bk, mu):
    mu_b = np.array([])
    for i in mu_bk:
        if len(i) == 6:
            num = (u[mu.index(i[0])]*u[mu.index(i[1])]*u[mu.index(i[2])])
            den = (u[mu.index(i[3])]*u[mu.index(i[4])]*u[mu.index(i[5])])
        elif len(i) == 14:
            num = (u[mu.index(i[0])]*u[mu.index(i[1])]*u[mu.index(i[2])]*u[mu.index(i[3])]*
                   u[mu.index(i[10])]*u[mu.index(i[11])]*u[mu.index(i[12])]*u[mu.index(i[13])])
            den = (u[mu.index(i[4])]*u[mu.index(i[5])]*u[mu.index(i[6])]*
                   u[mu.index(i[7])]*u[mu.index(i[8])]*u[mu.index(i[9])])
        mu_b = np.append(mu_b, num/den)
    return mu_b

def getGaussianClosure(u, mu_bk, mu):
    mu_b = np.array([])
    for i in mu_bk:
        if len(i) == 6:
            term_1 = u[mu.index(i[3])]*u[mu.index(i[2])]
            term_2 = u[mu.index(i[4])]*u[mu.index(i[1])]
            term_3 = u[mu.index(i[5])]*u[mu.index(i[0])]
            term_4 = 2*u[mu.index(i[3])]*u[mu.index(i[4])]*u[mu.index(i[5])]
        mu_b = np.append(mu_b, term_1+term_2+term_3-term_4)
    return mu_b

def derivativeLognormal(u, t, A_arr, B_arr, mu_bk, mu):
    return np.dot(A_arr, u) + np.dot(B_arr, getLognormalClosure(u, mu_bk, mu))

def derivativeGaussian(u, t, A_arr, B_arr, mu_bk, mu):
    return np.dot(A_arr, u) + np.dot(B_arr, getGaussianClosure(u, mu_bk, mu))