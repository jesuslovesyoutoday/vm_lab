from sympy import symbols, exp, simplify, I, re, im, plot_parametric

k = symbols("k", positive=True)
h = 0.01
t = 0.002

l = 1 + t*(2*exp(-I*k*3) - 9*exp(-I*k*2) + 18*exp(-I*k) - 11)/(6*h) + t**2*(-exp(-I*k*3) + 4*exp(-I*k*2) - 5*exp(-I*k) + 2)/(2*h**2) - t**3*(-exp(-I*k*3) + 3*exp(-I*k*2) - 3*exp(-I*k) + 1)/(6*h**3)

x = re(l)
y = im(l)

print(x)
print(y)
plot_parametric((x, y), xlabel='Re(l)', ylabel='Im(l)')
