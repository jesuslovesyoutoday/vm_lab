import numpy as np
import pandas as pd
from math import log, sin, cos
from matplotlib import pyplot as plt

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

def analit(x, t, T, X, tau, h):
    u = np.zeros((T, X))
    for i in range(T):
        for j in range(X):
            u[i][j] = log(1 + (x[j] - t[i])**2) + sin(x[j])
    return u

def numerical(x, t, T, X, tau, h):
    u = np.zeros((T, X))
    u[0] = np.log(1 + np.power(x, 2)) + np.sin(x)
    
    for i in range(len(u)):
        u[i][0] = log(1 + t[i]**2)
        u_x = -2*t[i]/(1+t[i]**2)
        u_x_1 = u_x + sin(1)
        u_xx = 2*(1 - t[i]**2)/(1+t[i]**2)**2
        u_xx_1 = u_xx - sin(1)
        u[i][1] = u[i][0] + h*u_x
        u[i][2] = u[i][1] + h*u_x_1
        #u[i][1] = u[i][0] + h*u_x + h**2*u_xx/2
        #u[i][2] = u[i][1] + h*u_x_1 + h**2*u_xx_1/2
    for i in range(1, T):
        n = i-1
        for j in range(3, X):
            u[i][j] = (u[n][j] + tau*(2*u[n][j-3] - 9*u[n][j-2] + 18*u[n][j-1] - 11*u[n][j])/(6*h)
                               + (tau**2)*(-u[n][j-3] + 4*u[n][j-2] - 5*u[n][j-1] + 2*u[n][j])/(2*h**2)
                               + (tau**3)*(-u[n][j-3] + 3*u[n][j-2] - 3*u[n][j-1] + u[n][j])/(6*h**3)
                               + tau*cos(x[j]) + 0.5*tau**2*sin(x[j]) - tau**3*cos(x[j])/6)
    return u
            
    
nodes = int(input()) + 1
h = 1/(nodes-1)
tau = h/10
nodes_t = int(1/tau)

x = np.linspace(0, 1, nodes)
t = np.linspace(0, 1, nodes_t)

u = analit(x, t, nodes_t, nodes, tau, h)
U = numerical(x, t, nodes_t, nodes, tau, h)

print_results(u[1], U[1], 0, 1, nodes)

plt.plot(x, u[0], label='аналит. t = 0')
plt.plot(x, U[0], label='числен. t = 0')
plt.plot(x, u[nodes_t-1], label='аналит. t = 0')
plt.plot(x, U[nodes_t-1], label='числен. t = 0')
plt.legend()
plt.grid()
plt.show()


            
