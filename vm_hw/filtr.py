import numpy as np
import matplotlib.pyplot as plt

def rho(p, p_0, rho0, cf):
    return rho0*(1 + cf*(p-p_0))

def progonka(a, b, c, d, nodes_x):
    alpha = [-b[0]/a[0]]
    beta = [d[0]/a[0]]
    for i in range(1, nodes_x):
        alpha.append(-b[i]/(a[i] + c[i] * alpha[i-1]))
        beta.append((d[i] - c[i] * beta[i-1])/(a[i] + c[i] * alpha[i-1]))
    p_pr = np.zeros(nodes_x)
    L = nodes_x - 1
    
    p_pr[L] = (d[L] - c[L]*beta[L-1])/(a[L] + c[L]*alpha[L-1])
    for i in range(L-1, -1, -1):
        p_pr[i] = alpha[i] * p_pr[i+1] + beta[i]
    return p_pr

L = 500

p0 = 1e7
p_inj = 1.5e7
p_prod = 5e6

k = 1e-13
mu = 1e-2
phi = 0.2
cf = 1e-9
rho0 = 1000
p_0 = 1.2e7

h = 5
nodes_x = int(L/h)
x = np.linspace(0, L, nodes_x)
tau = 60*60
T = 24* 60 * 60 * 10
nodes_t = int(T/tau)
time = np.linspace(0, T, nodes_t)

press = np.zeros((nodes_t, nodes_x))
press[0, :] = np.full((1, nodes_x), p0)

for t in range(1, nodes_t):
    
    a = [1]
    b = [0]
    c = [0]
    d = [p_inj]
    p = np.copy(press[t-1])
    for i in range(1, nodes_x-1):

        if (p[i] >= p[i+1]):
            rho2 = rho(p[i], p_0, rho0, cf)
        else:
            rho2 = rho(p[i+1], p_0, rho0, cf)
            
        if (p[i-1] >= p[i]):
            rho1 = rho(p[i-1], p_0, rho0, cf)
        else:
            rho1 = rho(p[i], p_0, rho0, cf)
        
        c.append(k*rho1/(mu*h**2))
        b.append(k*rho2/(mu*h**2))
        a.append(-c[i] - b[i] - phi*cf*rho0/tau)
        d.append(-phi*cf*rho0*p[i]/tau)
        
    a.append(1)
    b.append(0)
    c.append(0)
    d.append(p_prod)
    press[t] = np.copy(progonka(a, b, c, d, nodes_x))
    
    
step = int(nodes_t/10)
for i in range(0, nodes_t, step):
    plt.plot(x, press[i])
    plt.title(str(i*tau/(60*60*24)) + " ( дня)")
    plt.xlabel("x")
    plt.ylabel("p, Па")
    plt.savefig("res/" + str(i) + ".png")
    plt.clf()


