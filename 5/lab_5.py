import numpy as np
from math import pi, sqrt
from sympy.calculus.util import minimum, maximum
from sympy import diff, symbols, sin, cos, solveset, Min, Max, Interval

def func1(x):
    return x**2

def func(x):
    return sin(100*x)/(1+x)

# Step counting

def step_choose():
    x = symbols('x')                   
    f1 = diff(sin(100*x)/(1+x), x)      #
    f2 = diff(sin(100*x)/(1+x), x, 2)   # Derivates
    f4 = diff(sin(100*x)/(1+x), x, 4)   #   

    f1 = 100/(x + 1)                                                  #
    f2 = 2*(-5000*(-1)/(x + 1) + (-1)/(x + 1)**2)/(x + 1)             # Limit values
    f4 = 8*(12500000*1 - 15000*1/(x + 1)**2 + 3*1/(x + 1)**4)/(x + 1) #

    f = [f1, f2, f4]
    maximums = []
    l_b = 0
    u_b = pi/2
    e = 10e-6
    interval = Interval(l_b, u_b)
    inter = u_b - l_b

    for i in range(3):
        maximums.append(maximum(f[i], x, domain = interval))
    
    h1 = pow(e * 2 / (maximums[0]*inter), 1/2)    # Rectangle method (left)
    h2 = h1                                       # Rectangle method (right)
    h3 = pow(24 * e / (maximums[0]*inter), 1/3)   # Rectangle method (center)
    h4 = pow(12 * e / (maximums[1]*inter), 1/3)   # Trapezoidal method
    h5 = pow(2800 * e / (maximums[2]*inter), 1/4) # Rule 3/8
    h6 = pow(6480 * e/ (maximums[2]*inter), 1/5)  # Simpson method

    return [h1, h2, h3, h4, h5, h6]

def rect_left(h):
    x = 0
    I = 0
    while (x < pi/2):
        I += func(x) * h
        x += h
    return I

def rect_right(h):
    x = 0
    I = 0
    while (x < pi/2 - h):
        x += h
        I += func(x) * h
    return I

def rect_center(h):
    x = h/2
    I = 0
    while (x < pi/2):
        I += func(x) * h
        x += h
    return I

def trap(h):
    x1 = 0
    x2 = h
    I = 0
    while (x2 < pi/2):
        I += (func(x1) + func(x2)) * h / 2
        x1 += h
        x2 += h
    return I

def rule_3_8(h):
    x = 0
    x1_3 = h/3
    x2_3 = h * 2 / 3
    x1   = h
    I = 0
    while (x < pi/2):
        I += (func(x) + 3 * func(x1_3) + 3 * func(x2_3) + func(x1)) * h / 8
        x += h
        x1_3 += h
        x2_3 += h
        x1 += h
    return I
    
def simpson(h):
    x = 0
    x1_2 = h/2
    x1 = h
    I = 0
    while(x < pi/2):
        I += (func(x) + 4 * func(x1_2) + func(x1)) * h / 6
        x += h
        x1_2 += h
        x1 += h
    return I

def gauss():
    a = 0                                        #
    b = pi/2                                     #
    x1 = (pow(1/3, 1/2) * (b - a) + a + b)/2     # doesn't work properly
    x2 = (-pow(1/3, 1/2) * (b - a) + a + b)/2    #
    I = (func1(x1) + func1(x2)) * pi/4           #
    return I

def gauss_2(N):
    a = 0
    b = pi/2    
    h = (b - a)/N
    ai = 0
    bi = ai + h
    I = 0
    for i in range(N):
        x1 = (pow(1/3, 1/2) * (bi - ai) + ai + bi)/2
        x2 = (-pow(1/3, 1/2) * (bi - ai) + ai + bi)/2
        I += (bi-ai)*(func(x1) + func(x2))/2
        ai += h
        bi += h
    return I

def uzly(h):
    inter = pi/2
    print([inter/i for i in h])
    

h = step_choose()
print('шаги: ', h)
print('метод прямоугольников (левая)', rect_left(h[0]))
print('метод прямоугольников (правая)', rect_right(h[1]))
print('метод прямоугольников (центр)', rect_center(h[2]))
print('метод трапеций', trap(h[3]))
print('метод 3/8', rule_3_8(h[4]))
print('метод симпсона', simpson(h[5]))
print('метод гаусса', gauss_2(200))
uzly(h)


