import numpy as np
import math
from math import sin, cos, pi, sqrt, log2, log
from matplotlib import pyplot as plt

def make_n(r):
    data = [[1, 2]]
    k = 0
    for j in range(2, r+1):
        n = []
        for i in data[k]:
            n.append(i)
            n.append(pow(2, j) + 1 - i)
            i+=1
        k += 1
        data.append(n)
    return data[len(data)-1]

def make_tau(X, Y, hx, hy):
    min_mu = 4 * ((sin(pi/(2*X))/hx)**2 + (sin(pi/(2*Y))/hy)**2)
    #min_mu = 2*pi**2
    max_mu = 4 * ((cos(pi/(2*X))/hx)**2 + (cos(pi/(2*Y))/hy)**2)
    #max_mu = 4*(X**2 + Y**2)
    N_min = log(2/1e-3) / log((sqrt(max_mu) + sqrt(min_mu))/(sqrt(max_mu) - sqrt(min_mu)))
    print("N >= ", N_min)
    r = int(log2(N_min)) + 1
    #r = 5
    n = make_n(r)
    N = len(n)
    print("N: ", N, "\n")
    tauu = [2/((max_mu + min_mu) + (max_mu - min_mu)*cos(pi*(2*i-1)/(2*N))) for i in range(N)]
    tau = []
    for i in n:
        tau.append(tauu[i-1])
    return tau


def analit(x, y, X, Y):
    u = np.zeros((Y, X))
    for i in range(1, Y-1):
        for j in range(1, X-1):
            u[i][j] = cos(3*pi*x[j])*sin(4*pi*y[i])
    return u

def null_cond(x, y, X, Y, hx, hy):
    u = np.zeros((Y, X))
    """for i in range(1, Y-1):
        for j in range(1, X-1):
            u[i][j] = sin(4*pi*y[i])*cos(3*pi*x[j])
            #u[i][j] = 0"""
    return u

def numerical(x, y, X, Y, hx, hy, u, tau):
    for t in tau:
        u0 = np.copy(u)
        for i in range(1, Y-1):
            for j in range(1, X-1):
                u[i][j] = (u0[i][j] + t*((u0[i][j-1] - 2*u0[i][j] + u0[i][j+1])/(hx**2) 
                                       + (u0[i-1][j] - 2*u0[i][j] + u0[i+1][j])/(hy**2) 
                                       + 25*(pi**2)*cos(3*pi*x[j])*sin(4*pi*y[i])))
    return u


X = Y = 82

x = np.linspace(-0.5, 0.5, X)
y = np.linspace(0, 1, Y)
hx = 1/(X-1)
hy = 1/(Y-1)

tau = make_tau(X, Y, hx, hy)
#print(tau)
u = null_cond(x, y, X, Y, hx, hy)
U = numerical(x, y, X, Y, hx, hy, u, tau)

u_an = analit(x, y, X, Y)
u_num = U

"""plt.imshow(u_an)
plt.colorbar()
plt.show()
plt.imshow(u_num)
plt.colorbar()
plt.show()"""

np.set_printoptions(precision=4)

print("ANALITICAL:\n", u_an)
print("NUMERICAL:\n", u_num)
print("DELTA:\n", np.abs(u_an - u_num))
print("MAX DELTA: ", (np.abs(u_an-u_num)).max())
