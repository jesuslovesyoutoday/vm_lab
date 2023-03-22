import os
import csv
import numpy as np
from math import pi
import matplotlib.pyplot as plt

# TODO: animate results

def save_results(u, f_path, n):

    num = np.array([n])
    U = np.concatenate((num, u), axis=0)

    f = open(f_path, 'a')
    writer = csv.writer(f)
    writer.writerow(U)
    f.close()

def ugolok(u0, h, t, T, L, nodes, f_path, x):

    u_prev = u0
    save_results(u_prev, f_path, 0)
    for j in range(T):
    
        plt.plot(x, u_prev)
        plt.title("num: " + str(j))
        plt.grid()
        plt.show()
    
        u_next = np.zeros(nodes)
            
        for i in range(len(u0)-1):
            u_next[i+1] = u_prev[i]*t/h + u_prev[i+1]*(1-t/h)
                
        u_next[0] = u_next[len(u_next)-1]
        save_results(u_next, f_path, j+1)
        u_prev = u_next
        
    
f_path = "out.csv"

if (os.path.exists(f_path)):
    os.remove(f_path)

L = 20
T = 18
K = float(input("t/h: "))

h = 0.5
t = K*h
nodes = int(L/h + 1)

x = np.linspace(0, L, nodes)
u0 = np.sin(4*pi*x/L)

ugolok(u0, h, t, T, L, nodes, f_path, x)


