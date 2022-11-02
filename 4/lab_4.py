from matplotlib import pyplot as plt
from sympy import simplify
from sympy import Symbol
import pandas as pd
import numpy as np
import math

def solve(X, diff):
    x = Symbol('x')
    print(simplify(diff[0][0] + (x-X[0])*diff[1][0] + (x-X[0])*(x-X[1])*diff[2][0] + (x-X[0])*(x-X[1])*(x-X[2])*diff[3][0] + (x-X[0])*(x-X[1])*(x-X[2])*(x-X[3])*diff[4][0] + (x-X[0])*(x-X[1])*(x-X[2])*(x-X[3])*(x-X[4])*diff[5][0] + (x-X[0])*(x-X[1])*(x-X[2])*(x-X[3])*(x-X[4])*(x-X[5])*diff[6][0]), '\n')

a_0 = 0.983 
a_1 = -4.142
a_2 = 7.193
a_3 = -6.658
a_4 = 3.519
a_5 = -1.036
a_6 = 0.142

#------------------Разделенные разности-----------------#
def dev_dif(x, y):

    diff = []  
    diff.append(y)

    for i in range(1, len(data)):
        f = []
        y_ = diff[i-1]
        p = len(x) - len(y_)
        
        for j in range(len(y_) - 1):
            f.append((y_[j] - y_[j+1])/(x[j] - x[j+1+p]))
        diff.append(f)
    return diff
#----------------Производная интерполянта---------------#

def derivate(x0):
    return 6*x0**5*a_6 + 5*x0**4*a_5 + 4*x0**3*a_4 + 3*x0**2*a_3 + 2*x0*a_2 + a_1

#-----------------------Сплайн--------------------------#

def interpolate(x0, diff):

    m = len(diff[2])

    xx = [0]*(m+1)
    alpha = -1/4
    beta = []
    for i in range(m):
        beta.append(diff[2][i]/2)

    xx[m-1] = (diff[2][m-1] - beta[m-1]/2)/(alpha/2 + 2)
    xx[m] = 0
    i = m-2
    while (i != 0):
        xx[i] = xx[i+1] * alpha + beta[i+1]
        i -= 1

    d = [0, 0]

    for i in range(1, m):
        d.append((xx[i] - xx[i-1])/step)

    b = [0, 0]
    b.append((xx[1]*step/3 + diff[1][0]))

    for i in range(2, len(diff[1])):
        b.append(xx[i]*step/3 + xx[i-1]*step/6 + diff[1][i])

    ind = (np.abs(x - x0)).argmin()
    x1 = x[ind]
    dx = x0 - x1

    a0 = y[ind]
    a1 = b[ind]
    a2 = xx[ind]/2
    a3 = d[ind]/6

    xx = [i/2 for i in xx]
    d = [i/6 for i in d]

    print('a0', y)
    print('a1', b)
    print('a2', xx)
    print('a3', d)

    print('a0 = ', a0)
    print('a1 = ', a1)
    print('a2 = ', a2)
    print('a3 = ', a3)

    res = a0 + a1 * dx + a2 * dx**2 + a3 * dx**3
    print('x0 = ', x0, ' S(x0) = ', res)    

    return res

#-----------------------------------------------------#

def coeff(x, f):
    a0 = []
    a1 = []
    a2 = []
    a3 = []
    for i in range(len(x)-1):
        xi  = x[i]
        xi1 = x[i+1]
        fi  = f[i]
        fi1 = f[i+1]
        pi  = derivate(xi)
        pi1 = derivate(xi1) 
        
        a0.append((-pi1 * xi**2 *xi1 * (xi1 - xi) + fi1 * xi**2 * (3 * xi1 - xi) + fi * xi1**2 * (xi1 - 3 * xi) - pi * xi * xi1**2 * (xi1 - xi))/(xi1 - xi)**3)
        a1.append((pi1 * xi * (2 * xi1 + xi) * (xi1 - xi) - 6 * (fi1 - fi) * xi * xi1 + pi * xi1 * (xi1 + 2 * xi) * (xi1 - xi))/(xi1-xi)**3)
        a2.append((-pi1 * (xi1 - xi) * (xi1 + 2 * xi) + 3 * (fi1 - fi) * (xi1 + xi) - pi * (xi1 - xi) * (xi + 2 * xi1))/ (xi1 - xi)**3)
        a3.append((pi1 * (xi1 - xi) - 2 * (fi1 - fi) + pi * (xi1 - xi))/ (xi1 - xi)**3)        

    print('a0: ', a0)
    print('a1: ', a1)
    print('a2: ', a2)
    print('a3: ', a3, '\n')
    return[a0, a1, a2, a3]               

#------------------------------------------------------#

def res(x0, x, a0, a1, a2, a3):
    ind = (np.abs(x - x0)).argmin()
    dx = x0 - x[ind]
    return x0**3 * a3[ind] + x0**2 * a2[ind] + x0 * a1[ind] + a0[ind]

#------------------------------------------------------#

data = pd.read_csv('data.txt', sep = ' ')
x = (data.loc[:, "x"].to_numpy())
y = (data.loc[:, "y"].to_numpy())

diff = dev_dif(x, y)
solve(x, diff)
c = coeff(x, diff[0])

print("Введите точку x0:\n")

x0 = float(input())

print(res(x0, x, c[0], c[1], c[2], c[3]))




'''x_ = np.linspace(1, 2.7, 20)
y_ = []

for i in x_:
    y_.append(interpolate(i, diff))

plt.scatter(x, y)
plt.plot(x_, y_)
plt.show()'''
