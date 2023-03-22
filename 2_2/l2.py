import numpy as np
import scipy as sc
import pandas as pd
from math import exp, sqrt, sin
from matplotlib import pyplot as plt


def model(x0, u0, xa, xb, u1, k, f, q):
    
    ka = k[0]
    kb = k[1]
    fa = f[0]
    fb = f[1]
    qa = q[0]
    qb = q[1]
    
    mu_a = fa/qa
    mu_b = fb/qb
    
    l_a = sqrt(qa/ka)
    l_b = sqrt(qb/kb)
    
    A11 = exp(-l_a * x0) - exp(l_a * x0)
    A12 = exp(l_b * (2 - x0)) - exp(l_b * x0)
    A21 = ka * l_a * (exp(l_a * x0) + exp (-l_a * x0))
    A22 = kb * l_b * (exp(l_b * (2 - x0)) + exp(l_b * x0))
    
    B1 = mu_b - mu_a + (mu_a - u0) * exp(l_a * x0) - (mu_b - u1) * exp(l_b * (1 - x0))
    B2 = ka * l_a * (u0 - mu_a) * exp(l_a * x0) + kb * l_b * (u1 - mu_b) * exp(l_b * (1 - x0))
    
    c1 = (((u0 - mu_a) * A11 - B1) * A22 - ((u0 - mu_a) * A21 - B2) * A12) / (A11 * A22 - A12 * A21)
    c2 = (B1 * A22 - B2 * A12) / (A11 * A22 - A12 * A21)
    c3 = (B2 * A11 - B1 * A21) / (A11 * A22 - A12 * A21)
    c4 = (u1 - mu_b) * exp(l_b) - c3 * exp(2 * l_b)
    
    u = [c1 * np.exp(l_a * xa) + c2 * np.exp(-l_a * xa) + mu_a, c3 * np.exp(l_b * xb) + c4 * np.exp(-l_b * xb) + mu_b]
    
    return u

def array_division(x, x0):

    xa = []
    xb = []
    for i in x:
        if (i <= x0):
            xa.append(i)
        else:
            xb.append(i)
            
    return [xa, xb]

def norm(Y, Y_):
    Y = np.array(Y)
    Y_ = np.array(Y_)
    return (np.abs(Y - Y_))

def print_results(Y, Y_, a, b, sections):
    y = []
    y_= []
    x = np.linspace(a, b, 11)
    p = int(sections/10)
    if (sections == 10):
        y = Y
        y_= Y_
    else:
        for i in range(0, len(Y), p):
            y.append(Y[i])
            y_.append(Y_[i])

    d = norm(y, y_)
    df = pd.DataFrame({'x': x, 'Y_mod': y, 'Y': y_, 'delta': d})
    print(df)
    print('\nmax delta: ', d.max())
    df.to_csv('out.csv', sep = ',', index=False)


def vstr_progonka(nodes, x0, xa, xb, x, h, u0, u1):
    ua = np.zeros(len(xa))
    ub = np.zeros(len(xb))
    ua[0] = u0
    ub[len(xb)-1] = u1
    
    ind_a = len(xa) - 1
    ind_b = ind_a + 1
    
    ka = np.full(len(xa), 1)
    kb = np.exp(np.sin(xb))
    qa = np.full(len(xa), 1) 
    qb = np.full(len(xb), 2)
    fa = np.exp(xa)
    fb = np.exp(xb)
    
    aa = ka
    ab = kb
    ba = -2*ka - qa*h**2
    bb = -2*kb - qb*h**2
    ca = ka
    cb = kb
    da = -fa*h**2
    db = -fb*h**2

    alpha_a = np.zeros(ind_a)
    beta_a = np.zeros(ind_a)
    
    alpha_a[0] = -aa[1]/ba[1]
    beta_a[0]  = (da[1] - ca[1]*u0)/ba[1]
    
    a = nodes - 2 - ind_a  
    alpha_b = np.zeros(a)
    beta_b = np.zeros(a) 
    
    alpha_b[a-1] = -cb[a-1]/bb[a-1] 
    beta_b[a-1] = (db[a-1]-cb[a-1]*u1)/bb[a-1]
    

    for i in range(1, ind_a-1):
        alpha_a[i] = -aa[i+1]/(ba[i+1] + ca[i+1] * alpha_a[i-1])
        beta_a[i] = (da[i+1] - ca[i+1] * beta_a[i-1])/(ba[i+1] + ca[i+1] * alpha_a[i-1])
    for i in range(a-2, 0, -1):
        alpha_b[i] = -cb[i]/(bb[i] + ab[i] * alpha_b[i+1])
        beta_b[i] = (db[i] - ab[i] * beta_b[i+1])/(bb[i] + ab[i] * alpha_b[i+1])

    ua[ind_a] = ub[0] = (ka[ind_a]*beta_a[ind_a-2] + kb[0]*beta_b[1]) / (ka[ind_a]*(1 - alpha_a[ind_a-2]) + kb[0]*(1 - alpha_b[1]))

    for i in range(1, a):
        ub[i] = ub[i-1] * alpha_b[i] + beta_b[i]

    for i in range(ind_a - 1, 0, -1):
        ua[i] = ua[i + 1] * alpha_a[i-1] + beta_a[i-1]

    u = [ua, ub]
    
    return u

nodes = int(input("number of sections: ")) + 1
    
x = np.linspace(0, 1, nodes)
h = 1/(nodes - 1)

u0 = 1
u1 = 0
x0 = 1/(sqrt(2))

x_ = array_division(x, x0)
xa = np.array(x_[0])
xb = np.array(x_[1])

k = [1, exp(sin(x0))]
q = [1, 2]
f = [exp(x0), exp(x0)]

u_mod = np.array(model(x0, u0, xa, xb, u1, k, f, q), dtype=object)
u = np.concatenate((u_mod[0], u_mod[1]), axis=0)

u_ = np.array(vstr_progonka(nodes, x0, xa, xb, x, h, u0, u1), dtype=object)
u_num = np.concatenate((u_[0], u_[1]), axis=0)

print_results(u, u_num, 0, 1, nodes - 1)

plt.plot(x, u, label='модельная задача')
plt.plot(x, u_num, label='численное решение')
plt.legend()
plt.show()
