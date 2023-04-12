import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#def start():

def determinant(a, b, c):
    f = [0, 1, b[0]]
    for i in range(1, len(a)):
        f.append(b[i]*f[i+1] - c[i]*a[i-1]*f[i])
    #print('d: ', f[len(f)-1])
    return(f[len(f)-1])

def dihotomy(x, h, l1, l2, e, found): #l1 - d < 0 l2 - d > 0
    if found:
        return True
    l = (l1 + l2)/2
    #print('l: ', l)
    if (abs(l1 - l2) < e):
        found = True
        print('lambda: ')
        print(l)
        return True
    a, b, c = coeff(x, l, h)
    d = determinant(a, b, c)
    if (d > 0):
        #print('[', l1, ', ', l, ']')
        dihotomy(x, h, l1, l, e, found)
    elif (d < 0):
        #print('[', l, ', ', l2, ']')
        dihotomy(x, h, l, l2, e, found)
    #elif (d == 0):
        

def coeff(x, l, h):

    a = [1]
    b = [-1]
    c = [0]
    
    for i in range(1, len(x)-1):
        a.append((2-x[i]) - h)
        b.append(l*(x[i]**2)*(h**2)/(1+x[i]) - 2*(2-x[i]) + h)
        c.append((2-x[i]))
    
    a.append(0)
    b.append(h)
    c.append(0)
    return (a, b, c)


nodes = 321

x = np.linspace(0, 1, nodes)
l0 = 0 
e = 1e-6
h = 1/(nodes-1)
a, b, c = coeff(x, l0, h)

#df = pd.DataFrame({'x': x, 'a': a, 'b': b, 'c': c})
#print(df)

"""l0 = 1 
a, b, c = coeff(x, l0, h)

df = pd.DataFrame({'x': x, 'a': a, 'b': b, 'c': c})
print(df)"""


lambdas = []
determinants = []

for i in range(4):
    a, b, c = coeff(x, l0, h)
    d0 = determinant(a, b, c)
    if (d0 < 0):
        l = l0
        d = determinant(a, b, c)
        while (d < 0):
            determinants.append(d)
            lambdas.append(l)
            l += 1
            a, b, c = coeff(x, l, h)
            d = determinant(a, b, c)
        #print('d<0 while l < ', l)
        dihotomy(x, h, l0, l, e, False)
    elif (d0 > 0):
        l = l0
        d = determinant(a, b, c)
        while (d > 0):
            determinants.append(d)
            lambdas.append(l)
            l += 1 
            a, b, c = coeff(x, l, h)
            d = determinant(a, b, c)
        #print('d>0 while l < ', l)
        dihotomy(x, h, l, l0, e, False)
    elif (d0==0):
        l = l0
        print('lambda: ', l0)
    l0 = l + 1



plt.plot(lambdas, determinants)
plt.grid()
plt.show()


