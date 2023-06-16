import sympy as sp
import numpy as np
from scipy import constants
from scipy import integrate
from math import sqrt, sin
from matplotlib import pyplot as plt

def integrand(x, t, En, x1, sigma, a, h, m):
    return N * np.exp(-(x-x1)**2/(2*sigma**2)) * sqrt(2/a) * np.exp(1j*En*t/h)*(np.sin(sqrt(2*m*En/(h**2))*x))

def psi_0(x):
    return N * np.exp(-(x-x1)**2/(2*sigma**2))

a = 0.1
x1 = 0.05
sigma = 1/100
N = 1/100
m = constants.m_e
h = constants.hbar
n = 10

X = np.linspace(0, a, 100)
T = np.linspace(0, 10, 200)

"""plt.plot(X, np.exp(-(X-x1)**2/(2*sigma**2)))
plt.plot(X, np.power(np.abs(np.exp(-(X-x1)**2/(2*sigma**2))), 2))
plt.show()"""


psi = np.zeros(len(X))
rho = []

for t in T:
    psi = np.zeros(len(X))
    for i in range(1, n):
        En = (np.pi*i*h/a)**2 * 1/(2*m)
        #plt.plot(X, N*np.exp(-(X-x1)**2/(2*sigma**2)))
        #plt.plot(X, sqrt(2/a) * np.exp(1j*En*t/h)*(np.sin(sqrt(2*m*En/(h**2))*X)))
        #plt.show()
        cn, err = integrate.quad(integrand, -np.inf, np.inf, args=(t, En, x1, sigma, a, h, m))
        #plt.plot(X, psi)
        #plt.plot(X, cn * sqrt(2/a) * np.exp(-1j*En*t/h)*(np.sin(np.sqrt(2*m*En/(h**2))*X)))
        #plt.show()
        psi = psi + cn * sqrt(2/a) * np.exp(-1j*En*t/h)*(np.sin(np.sqrt(2*m*En/(h**2))*X))
    rho.append(np.power(np.abs(psi), 2))

for i in range(len(T)):

    plt.plot(X, rho[i])
    plt.xlabel("x, m")
    plt.ylabel("|Psi|^2")
    plt.savefig("pic/" + str(i) + ".png")
    plt.clf()






"""def null_cond(x):
    return N*np.exp(-np.power(x, 2)/(2*sigma**2))

def rho(psi):
    return (sp.Abs(psi))**2

def psi_n(X, T, en):
    return np.exp(-1j*en*T/h)*sqrt(2/a)*np.sin(np.sqrt(2*m*en/h**2)*X)

def psi(X, T, cn, en):
    psi = np.array([c * psi_n(X, T, en) for c in cn])
    return np.sum(psi) 


X = np.linspace(-a, a, 100)

#psi_x_0 = null_cond(x)

x, t, n = sp.symbols("x t n")
Psi_x_0 = N * sp.exp(-x**2/(2*sigma**2))

En = (sp.pi*n*h/a)**2 * 1/(2*m)
Psi_n = sqrt(2/a)*sp.exp(-sp.I*En*t/h)*(sp.sin(sp.sqrt(2*m*En/(h**2))*x))

expr = Psi_x_0 * sp.conjugate(Psi_n)
print(expr)
cn = sp.integrate(expr, x)
print(cn)
"""
