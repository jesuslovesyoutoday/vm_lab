import csv
from math import e
from math import pi
from math import cos
import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return pow(10, x)

def first_norma(y, Y):
    y = np.array(y)
    dy = y-Y
    return (np.abs(dy)).max()

def dev_dif(x, y):

    diff = []  
    diff.append(y)

    for i in range(1, len(x)):
        f = []
        y_ = diff[i-1]
        p = len(x) - len(y_)
        
        for j in range(len(y_) - 1):
            f.append((y_[j] - y_[j+1])/(x[j] - x[j+1+p]))
        diff.append(f)
    return diff

def interpolate(diff, x, nodes):
    x_c = [1]
    Y = []
    y_ = []
    """for i in range(len(x)):
        res = diff[0][0]
        x_c = [1]
        for j in range(1, len(nodes)):
            x_c.append(x_c[j-1] * (x[i] - nodes[j-1]))
            res += x_c[j] * diff[j][0]
        Y.append(res)"""
    for i in range(len(nodes)):
        res = diff[0][0]
        x_c = [1]
        for j in range(1, len(nodes)):
            x_c.append(x_c[j-1] * (nodes[i] - nodes[j-1]))
            res += x_c[j] * diff[j][0]
        Y.append(res)
        y_.append(func(nodes[i]))
    return (Y,y_,nodes)

def uni_grid(N, x):
    
    nodes = np.linspace(-1, 1, N)
    y_n = [func(x) for x in nodes]
    #print('x_ravn: ', nodes)

    diff = dev_dif(nodes, y_n)
    return interpolate(diff, x, nodes)

def chebishev_grid(N, x):
    
    nodes = np.array([cos((pi + 2 * pi * i)/(2*N)) for i in range(N)])
    nodes = np.flip(nodes)
    #print('x_cheb: ', nodes)
    y_n = [func(x) for x in nodes]
    
    diff = dev_dif(nodes, y_n)
    return interpolate(diff, x, nodes)


Nmax = 100
   
x = np.linspace(-1, 1, 10000)
y = np.array([func(i) for i in x])

n  = []
n1 = []
N_ = np.arange(3, Nmax, 1)

for N in N_:
    print("Counting ", N, "nodes")
    Y, y_, nodes = uni_grid(N, x)
    n.append(first_norma(y_, Y))
    
    Y1, y_1, nodes_1 = chebishev_grid(N, x)
    n1.append(first_norma(y_1, Y1))
    
    """plt.plot(x, y, color = 'g')
    plt.scatter(x, Y, color = 'r')
    plt.scatter(x, Y1, color = 'b')
    plt.show()
    plt.grid()"""

with open('uni.csv', 'w') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerows(zip(N_,n))
       
with open('non_uni.csv', 'w') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerows(zip(N_,n1))

plt.scatter(N_, n1, label = 'chebishev_grid')
plt.scatter(N_, n, label = 'uni_grid')
plt.legend()
plt.grid()
plt.show()
