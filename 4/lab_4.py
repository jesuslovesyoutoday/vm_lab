from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import math

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

    print('a0 = ', a0)
    print('a1 = ', a1)
    print('a2 = ', a2)
    print('a3 = ', a3)

    res = a0 + a1 * dx + a2 * dx**2 + a3 * dx**3
    print('x0 = ', x0, ' S(x0) = ', res)    

    return res

data = pd.read_csv('data.txt', sep = ' ')
x = (data.loc[:, "x"].to_numpy())
y = (data.loc[:, "y"].to_numpy())

step = x[1] - x[0]

print("Введите точку x0:\n")
x0 = float(input())

diff = dev_dif(x, y)
interpolate(x0, diff)

'''x_ = np.linspace(1, 2.7, 20)
y_ = []

for i in x_:
    y_.append(interpolate(i, diff))

plt.scatter(x, y)
plt.plot(x_, y_)
plt.show()'''
