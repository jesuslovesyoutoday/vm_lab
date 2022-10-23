import pandas as pd
import numpy as np
import math

data = pd.read_csv('data.txt', sep = ' ')
x = (data.loc[:, "x"].to_numpy())
y = (data.loc[:, "y"].to_numpy())
data = data.to_numpy()
step = x[1] - x[0]


#-----------------------Интерполяция полиномом------------------------#


diff = []  
diff.append(y)

x0 = 0.95

for i in range(1, len(data)):
    f = []
    y_ = diff[i-1]
    p = len(x) - len(y_)
    
    for j in range(len(y_) - 1):
        f.append((y_[j] - y_[j+1])/(x[j] - x[j+1+p]))
    diff.append(f)

x_c = []
x_c.append(1)

for i in range(1, len(diff)):
    x_c.append(x_c[i-1] * (x0 - x[i-1]))

res = 0

for i in range(len(diff)):
    res +=  diff[i][0] * x_c[i]
print('Интерполяция методом Ньютона:', res)


#-----------------Интерполяция куюическим сплайном---------------------#


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

res1 = y[ind] + b[ind]*dx + (xx[ind] * dx**2)/2 + (d[ind] * dx**3)/6

print('Интерполяция кубическим сплайном', res1)


#---------------------------Погрешности--------------------------------#


N = len(data)
delta_m = 10**6 * step / (N * math.factorial(N))
L = 2**N/(math.exp(1)*(N-1)*math.log(N-1))

yy_1 = []
yy_2 = []

for i in range(m):
    x0 = x[i]
    x_c = [1]
    for j in range(1, len(diff)):
        x_c.append(x_c[j-1] * (x0 - x[j-1]))
    res = 0
    for k in range(len(diff)):
        res +=  diff[k][0] * x_c[k]
    yy_1.append(res)
    
    dx = step
    res1 = y[i] + b[i]*dx + (xx[i] * dx**2)/2 + (d[i] * dx**3)/6
    
    yy_2.append(res1)

y = y[:len(yy_1)]
print(yy_1)    
d1 = (y - yy_1).max()
d2 = (y - yy_2).max()

delta_v1 = d1 * L
delta_v2 = d2 * L

delta_1 = delta_m + delta_v1
delta_2 = delta_m + delta_v2

print('Ошибка интерполяции методом Ньютона:', delta_1)
print('Ошибка интерполяции сплайном:', delta_2)
