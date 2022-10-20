import csv
from math import pi
from math import cos
import numpy as np
import matplotlib.pyplot as plt

def first_norma(y, Y):
    dy = y-Y
    return np.abs(dy).max()

def newton_interpolate(diff, nodes, x):

    for i in range(1, len(nodes)):
        f = np.array([])
        y_ = diff[i-1]
        p = len(nodes) - len(y_)
        
        for j in range(len(y_) - 1):
            f = np.append(f, (y_[j] - y_[j+1])/(nodes[j] - nodes[j+1+p]))
        diff.append(f)
        
    Y = np.array([])
    for x0 in x:
        x_c = np.array([])
        x_c = np.append(x_c, 1)

        for i in range(1, len(diff)):
            x_c = np.append(x_c, x_c[i-1] * (x0 - x[i-1]))

        res = 0

        for i in range(len(diff)):
            res +=  diff[i][0] * x_c[i]
        Y = np.append(Y, res) 
        
    return Y

def uni_grid(N, x):
    
    nodes = np.arange(-1, 1, 2/N)
    y_n   = [pow(10, x) for x in nodes]
    diff = [y_n]
    
    return newton_interpolate(diff, nodes, x)
    
 
def chebishev_grid(N, x):
    
    nodes = np.array([cos((pi + 2 * pi * i)/(2*N)) for i in range(N)])
    y_n = [pow(10, x) for x in nodes]
    diff = [y_n]
    
    return newton_interpolate(diff, nodes, x)

   
x = np.arange(-1, 1, 2/10000)
y = np.array([pow(10, i) for i in x])

n  = []
n1 = []
N_ = np.arange(3, 100, 1)

for N in N_:
    print("Counting ", N, "nodes")
    
    Y = uni_grid(N, x)
    
    """print(Y)
    print(y)
    print(x)
    with open('test1.csv', 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(zip(x,Y))
    with open('test2.csv', 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(zip(x,y))
    break"""
    
    n.append(first_norma(y, Y))
    
    Y1 = chebishev_grid(N, x)
    n1.append(first_norma(y, Y1))

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
