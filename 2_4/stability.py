import sympy as sm
from sympy import exp, simplify, Abs, re, im, plot_parametric
from IPython.display import display
from sys import exit

sm.init_printing(use_latex=True)

h, t, tau, lam, k, m, n = sm.symbols("h t tau lambda k m n", real=True)

u6x = 0
# u5x = 0
#h = 0.01
#tau = 0.001


un1 = lam**(n+1) * exp(1j*k*m*h)
umn = lam**n * exp(1j*k*m*h)
um1 = lam**n * exp(1j*k*(m+1)*h)
um2 = lam**n * exp(1j*k*(m+2)*h)
um3 = lam**n * exp(1j*k*(m+3)*h)

exp = umn - un1 + tau/(6*h)*(tau/2 + tau + 6)*(2*um3 - 9*um2 + 18*um1 - 11*umn) + \
    tau**2/(2*h**2) *(tau + tau + 6)*(tau+6) *(-um3 + 4*um2 - 5*um1 + 2*umn) - \
    tau**3/(6*h**3)*(tau+6)**3 * (-um3 + 3*um2 - 3*um1 + umn)
    
# exp = exp

# print(sm.simplify(exp))
res = sm.solvers.solve(exp, lam)[1]
# print(sm.solvers.solve(exp, lam))
# display(simplify(Abs(res.expand(complex=True))))
# display(Abs(res.expand(complex=True)))
# print(str(Abs(res.expand(complex=True))).replace('**', '^').replace('tau', 't'))
print('(' + 
      str(re(res.expand(complex=True))).replace('**', '^').replace('tau', 'T') + 
      ', ' + 
      str(im(res.expand(complex=True))).replace('**', '^').replace('tau', 'T') + 
      ')')

cur = 0.8

H = 1/20
T = cur*H

display(lam)
x = re(res.expand(complex=True))
x = x.subs({h:H, tau:T})
print("Re:", x)
#display(re(res.expand(complex=True)))
y = im(res.expand(complex=True))
y = y.subs({h:H, tau:T})
print("Im:", y)
#display(im(res.expand(complex=True)))

plot_parametric((x, y), xlabel='Re(l)', ylabel='Im(l)')
exit(0)

res2 = sm.solvers.solve(res, tau)[1]/h
display(simplify(res2))
print(simplify(Abs(res2)))
print(str(simplify(res2)).replace('I', 'i'))
