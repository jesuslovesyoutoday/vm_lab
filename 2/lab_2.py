import math
import numpy as np

#-----------------Начальные данные-----------------#

gamma_0 = 5/3
rho_0   = 1.694 * pow(10,-4)
P_0     = 1.013 * pow(10, 6)
U_0     = 0

gamma_3 = 7/5
C_3     = 3.6537 * pow(10, 4)
P_3     = 1.6768 * pow(10, 6)
U_3     = 0

rho_3 = gamma_3 * P_3 / C_3**2 

#----------------Поиск коэффициентов---------------#

alpha_0 = (gamma_0 + 1) / (gamma_0 - 1) 
n       = 2 * gamma_3 / (gamma_3 - 1)
n = int(n)
print('n =', n)
X = P_3/P_0

#Z = pow(P_1/P_3, 1/n)

mu = (U_3 - U_0) * math.sqrt((gamma_0 - 1) * rho_0 / (2 * P_0))
nu = (2 / (gamma_3 - 1)) * math.sqrt(gamma_3 * (gamma_0 - 1) * P_3 * rho_0 / (2 * P_0 * rho_3))

a0 = X**2
a1 = - alpha_0 * nu**2 * X
a2 = 2 * alpha_0 * nu * (mu + nu) * X
a3 = -(2 + ((mu + nu)**2) * alpha_0) * X
a4 =  - nu**2
a5 = 2 * nu * (nu + mu)
a6 = -((nu + mu)**2) + 1

print('a0 = ', a0 , ' * Z**(2*n)')
print('a1 = ', a1 , ' * Z**(n+2)')
print('a2 = ', a2 , '* Z** (n+1)')
print('a3 = ', a3 , '* Z**n')
print('a4 = ', a4 ,' * Z**2')
print('a5 = ', a5 , '* Z')
print('a6 = ', a6, '\n')

#---------------Локализация корней-----------------#

coef = np.array([a0, a1, a2, a3, a4, a5, a6])
A = np.max(coef[1:])
B = np.max(coef[0:len(coef) - 2])

print('A = ', A)
print('B = ', B)

x1 = abs(a6) / (abs(a6) + B)
x2 = 1 + A/(abs(a0))

print(x1, '< x <', x2)

#--------------------Поиск корней-------------------#

def f(x):
    f = (a0 * x**(2*n) + a1 * x**(n+2) + a2 * x**(n+1) +
        + a3 * x**n + a4 * x**2 + a5 * x + a6)
    return f 

def f_1(x):
    f_1 = (2*n * a0 * x**(2*n - 1) + (n + 2) * a1 * x**(n+1) +
         + (n + 1) * a2 * x**(n) + n * a3 * x**(n-1) + 2 * a4 * x +
         + a5)
    return f_1

def newton(x0):
    x = [x0]
    i = 0
    y = f(x[i])
    while (abs(y) > 10**(-14)):
        x_ = x[i] - f(x[i])/f_1(x[i])
        y = f(x_)
        x.append(x_)
        i += 1
    return x[len(x) - 1]

h = 0.001
xi = [x1]
roots = []

for i in range(int((x2 - x1)/h)):
    xi.append(xi[i] + h)

for i in range(len(xi)-1):
    if (f(xi[i])*f(xi[i+1]) < 0):
        roots.append(newton((xi[i+1] + xi[i])/2))

for i in range(len(roots)):
    roots[i] = int(roots[i] * 10000)/10000

print('\nroots: ', roots) 

#----------------Вычисление параметров------------------#
P_1 = []

for i in roots:
    P_1.append(int(P_3*(i**n)*10000)/10000)

print('\nЗначения P1: ', P_1)

P_2 = P_1
C_0 = math.sqrt(gamma_0 * P_0 / rho_0)

D_0 = []
U_1 = []
U_2 = []

for P1 in P_1:
    C_2 = C_3 * pow(P1/P_3, (gamma_3 - 1)/(2 * gamma_3))
    U2 = U_3 + 2 * (C_3 - C_2)/(gamma_3 - 1)
    U_2.append(U2)
    D_0.append(U_0 - (1/rho_0) * (P1 - P_0)/(U_0 - U2))
U_1 = U_2

for i in range(len(D_0)):
    D_0[i] = int(D_0[i] * 10000)/10000

print('\nЗначения D0: ', D_0)

#-------------Определяем физический случай--------------#

lambda_0 = []
lambda_1 = []

for i in range(len(D_0)):

    C_1 = C_3 * pow(P_2[i]/P_3, (gamma_3 - 1)/(2*gamma_3))
    #rho_1 = gamma_1 * P_1[i] / C_1**2
    rho_1 = rho_0 * (U_0 - D_0[i]) / (U_1[i] - D_0[i])
    a_kr = math.sqrt((P_1[i] - P_0)/(rho_1 - rho_0))
    
    lambda_0.append(abs(U_0 - D_0[i])/a_kr)
    lambda_1.append(abs(U_1[i] - D_0[i])/a_kr)
    
    print('\na_кр_', i+1,'= ', a_kr)
    print('lambda_0', i+1, '= ', lambda_0[i])
    print('lambda_1', i+1, '= ', lambda_1[i])
    print('l_1 * l_2', i+1, '= ', lambda_0[i]*lambda_1[i], '\n')


