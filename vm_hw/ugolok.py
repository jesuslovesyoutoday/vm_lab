import os
import csv
import numpy as np
from math import pi
import matplotlib.pyplot as plt

def save_results(u, f_path, n):

    num = np.array([n])
    U = np.concatenate((num, u), axis=0)

    f = open(f_path, 'a')
    writer = csv.writer(f)
    writer.writerow(U)
    f.close()

def ugolok(u0, h, t, T, L, nodes, f_path, x, step_plot):

    u_prev = u0
    save_results(u_prev, f_path, 0)
    j = 0
    count = 0
    while j < T:
    
        if (count % step_plot) == 0:
            plt.plot(x, u_prev, label="num: " + str(j))
            plt.grid()
    
        u_next = np.zeros(nodes)
            
        for i in range(len(u0)-1):
            u_next[i+1] = u_prev[i]*t/h + u_prev[i+1]*(1-t/h)
                
        u_next[0] = u_next[len(u_next)-1]
        save_results(u_next, f_path, j+1)
        u_prev = u_next
        j += t
        count += 1
        
    plt.legend()
    plt.show()
        
    
f_path = "out_ugolok.csv"

if (os.path.exists(f_path)):
    os.remove(f_path)

L = 20
T = 18
K = float(input("t/h: "))
moments = [0, 5, 10, 15, 18]

h = 0.5
t = K*h
nodes = int(L/h + 1)
step_plot = int(int(T/t)/5)

x = np.linspace(0, L, nodes)
u0 = np.sin(4*pi*x/L)

ugolok(u0, h, t, T, L, nodes, f_path, x, step_plot)


