import pandas as pd
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def print_results(Y, Y_, a, b, sections):
    y = [[], []]
    y_= [[], []]
    x = np.linspace(a, b, 11)
    p = int(sections/10)
    if (sections == 10):
        y = Y
        y_= Y_
    else:
        for i in range(0, len(Y[0]), p):
            y[0].append(Y[0][i])
            y[1].append(Y[1][i])
            y_[0].append(Y_[0][i])
            y_[1].append(Y_[1][i]) 
    d = norm(y, y_)
    df = pd.DataFrame({'x': x, 'Y1': y[0], 'Y*1': y_[0], 'delta1': d[0], 'Y2': y[1], 'Y*2': y_[1], 'delta2': d[1]})
    print(df)
    df.to_csv('out.csv', sep = ',', index=False)

def norm(Y, Y_):
    Y = np.array(Y)
    Y_ = np.array(Y_)
    return (np.abs(Y[0] - Y_[0]).max(), np.abs(Y[1] - Y_[1]).max())

def func(Y, J):
    return np.matmul(J,Y.transpose())

def analit(x, A, B):
    h1 = [5, 2]
    h2 = [-5, 2]
    y1 = (A/10 + B/4)*np.exp(x)*h1[0] + (B/4 - A/10)*np.exp(-199*x)*h2[0]
    y2 = (A/10 + B/4)*np.exp(x)*h1[1] + (B/4 - A/10)*np.exp(-199*x)*h2[1]
    Y = np.array([y1, y2])
    return (Y)

def explicit_euler(x, NU, J, nodes, h):
    Y = np.array([NU])
    for i in range(nodes-1):
        Y = np.concatenate((Y,[Y[i] + h * np.matmul(J, Y[i].transpose())]), axis=0)
    return (Y)


def euler_pereschet(x, NU, J, nodes, h):
    Y = np.array([NU])
    for i in range(nodes-1):
        f = func(Y[i], J)
        y = Y[i] + h * f
        f_ = func(y, J)
        Y = np.concatenate((Y, [Y[i] + h * (f + f_)/2]), axis=0)
    return (Y)

def implicit_euler(x, NU, J, nodes, h):
    Y1 = [NU[0]]
    Y2 = [NU[1]]
    for i in range(nodes-1):
        #y1, y2 = sp.symbols('y1, y2')                                                                          #
        #res = list(sp.linsolve([y1 - Y1[i] + 99*h*y1 - 250*h*y2, y2 - Y2[i] - 40*h*y1 + 99*h*y2], (y1, y2)))   # <-- this is 
        #y1 = res[0][0]                                                                                         # really slow
        #y2 = res[0][1]                                                                                         #
        
        #y1 = (Y1[i]*(1+99*h) + Y2[i]*250*h)/(-199*h*h + 198*h + 1)
        #y2 = (y1*(99*h+1) - Y1[i])/(250*h)

        y1 = (250*h*Y2[i] + Y1[i]*(1+99*h))/((1+99*h)**2 - 1000*(h**2))
        y2 = (40*h*y1 + Y2[i])/(1+99*h) 
        Y1.append(y1)
        Y2.append(y2)        
    return ([Y1, Y2])


sections = int(input("sections = "))
nodes = sections+1
A = int(input("A = "))
B = int(input("B = "))
NU = np.array([A, B])
a = 0
b = 1
l = b - a
x = np.linspace(a, b, nodes)
h = float(l/(nodes))
print('h = ', h)
J = np.array([[-99, 250],[40, -99]])


Y  = analit(x, A, B)
Y_ = np.array(implicit_euler(x, NU, J, nodes, h))
#Y_ = explicit_euler(x, NU, J, nodes, h).transpose()

print_results(Y, Y_, a, b, sections)
plt.plot(x, Y[1], 'r', label = 'аналитическое')
plt.plot(x, Y_[1], 'b', label = 'численное')
plt.legend()
plt.grid()
plt.show()    
