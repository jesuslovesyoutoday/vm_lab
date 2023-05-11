import numpy as np
import pandas as pd
from math import sqrt
import matplotlib.pyplot as plt

def norm(Y, Y_):
    Y = np.array(Y)
    Y_ = np.array(Y_)
    return (np.abs(Y - Y_))

def print_results(Y, Y_, a, b, X):
    y = np.zeros((6, 6))
    y_= np.zeros((6, 6))
    p = int(X/6)
    if (X == 6):
        y = Y
        y_= Y_
    else:
        k = 0
        for i in range(0, len(Y), p):
            l = 0
            for j in range(0, len(Y[0]), p):
                y[k][l] = Y[i][j]
                y_[k][l] = Y_[i][j]
                l += 1
            k += 1
            
    d = norm(y, y_)
    np.set_printoptions(precision=4)
    print("ANALITICAL:\n", y, '\n')
    print("NUMERICAL:\n", y_, "\n")
    print("DELTA:\n", d, '\n')
    print('MAX DELTA: ', d.max())

def progonka(a, b, c, d, nodes):
    alpha = [-a[0]/b[0]]
    beta = [d[0]/b[0]]
    
    for l in range(1, nodes):
        alpha.append(-a[l]/(b[l]+c[l]*alpha[l-1]))
        beta.append((d[l]-c[l]*beta[l-1])/(b[l]+c[l]*alpha[l-1]))
    u_p = np.zeros(nodes)
    L = nodes - 1

    #np.set_printoptions(precision=4)
    u_p[L] = (d[L] - c[L]*beta[L-1])/(b[L] + c[L]*alpha[L-1])
    for l in range(L-1, -1, -1):
        u_p[l] = alpha[l] * u_p[l+1] + beta[l] 
    return u_p
    

def analit(x, y, X, Y, t):
    u = np.zeros((Y, X))
    for i in range(Y):
        for j in range(X):
            u[i][j] = np.power((1 + x[j] + y[i]), 4)/(21 - 20*t)**2
    return u

def numerical(x, y, hx, hy, X, Y, t, tau, e):
    u = analit(x, y, X, Y, 0) # НУ по времени
    
    for n in range(1, len(t)):

        u[:, 0]   =  np.power((1 + y), 4)/(21 - 20*t[n])**2 #
        u[:, X-1] =  np.power((2 + y), 4)/(21 - 20*t[n])**2 # <- ГУ по пр-ву
        u[0, :]   =  np.power((1 + x), 4)/(21 - 20*t[n])**2 #
        u[Y-1, :] =  np.power((2 + x), 4)/(21 - 20*t[n])**2 #

        u0 = np.copy(u) # на n-ом шаге
        prodoljaem = True
        u_0 = np.copy(u)  # на полушаге
        u_00 = np.copy(u) # на k-ой итерации 

        while (prodoljaem):

            for m in range(1, Y-1):
                a = np.zeros(X)
                b = np.zeros(X)
                c = np.zeros(X)
                d = np.zeros(X)

                for l in range(1, X-1):
                    a[l] = -(sqrt(u_0[m][l+1]) + sqrt(u_0[m][l])) * tau/(2*hx**2)
                    c[l] = -(sqrt(u_0[m][l-1]) + sqrt(u_0[m][l])) * tau/(2*hx**2)
                    b[l] = 1 - a[l] - c[l]
                    d[l] = u0[m][l]
                a[0] = 0
                b[0] = 1
                c[0] = 0
                d[0] = u0[m][0]
                a[X-1] = 0
                b[X-1] = 1
                c[X-1] = 0
                d[X-1] = u0[m][X-1]
                
                u[m, 1:X-1] = progonka(a, b, c, d, X)[1:X-1]
                #u[m, :] = progonka(a, b, c, d, X)
                #u[m, 1:X-1] = progonka(a[1:X-1], b[1:X-1], c[1:X-1], d[1:X-1], X-2)
                u_0 = np.copy(u)
            #print("после 1")
            #plt.imshow(u)
            #plt.colorbar()
            #plt.show()

            for l in range(1, X-1):
                a = np.zeros(Y)
                b = np.zeros(Y)
                c = np.zeros(Y)
                d = np.zeros(Y)

                for m in range(1, Y-1):
                    a[m] = -(sqrt(u0[m+1][l]) + sqrt(u0[m][l])) * tau/(2*hy**2)
                    c[m] = -(sqrt(u0[m-1][l]) + sqrt(u0[m][l])) * tau/(2*hy**2)
                    b[m] = 1 - a[m] - c[m]
                    d[m] = u_0[m][l]
                a[0] = 0
                b[0] = 1
                c[0] = 0
                d[0] = u_0[0][l]
                a[X-1] = 0
                b[X-1] = 1
                c[X-1] = 0
                d[X-1] = u_0[Y-1][l]

                u[1:Y-1, l] = progonka(a, b, c, d, Y)[1:Y-1]
            #print("после 2")
            #plt.imshow(u)
            #plt.colorbar()
            #plt.show()
                #u[:, l] = progonka(a, b, c, d, Y)
                #u[1:Y-1, l] = progonka(a[1:Y-1], b[1:Y-1], c[1:Y-1], d[1:Y-1], Y-2)

            k = 0
            for m in range(1, Y-1):
                for l in range(1, X-1):
                    if (u[m][l] != 0):
                        if (abs((u[m][l] - u_00[m][l])/u[m][l]) > e):
                            k+=1
            if(k==0):
                prodoljaem = False
            u_0 = np.copy(u)
            u_00 = np.copy(u)
    return u



X = Y = int(input("Enter nodes: "))
N = int(input("Enter time steps: "))
tau = 1/N
sections = X - 1
hx = hy = 1/(sections)
T = 1
e = 0.001

x = np.linspace(0, 1, X)
y = np.linspace(0, 1, Y)
t = np.linspace(0, 1, N)

u_an = analit(x, y, X, Y, T)
u_num = numerical(x, y, hx, hy, X, Y, t, tau, e)

d = norm(u_an, u_num)
np.set_printoptions(precision=4)
print("ANALITICAL:\n", u_an, '\n')
print("NUMERICAL:\n", u_num, "\n")
print("DELTA:\n", d, '\n')
print('MAX DELTA: ', d.max())
    
plt.imshow(u_an)
plt.colorbar()
plt.show()
plt.imshow(u_num)
plt.colorbar()
plt.show()

"""d = [-5, -18, -40, -27]
a = [0, 1, 1, 1]
b = [2, 10, -5, 4]
c = [1, -5, 2, 0]

print(progonka(a, b, c, d))
d = [-5, -18, -40, -27]
b = [2, 10, -5, 4]
a = [1, -5, 2, 0]
c = [0, 1, 1, 1]

print(progonka(a, b, c, d))"""

