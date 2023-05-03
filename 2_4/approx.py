import sympy as sm
from IPython.display import display

sm.init_printing(use_latex=True)

# un1, umn, um1, um2, um3 = sm.symbols('u_m^{n+1} u_m^n u_{m+1} u_{m+3} u_{m+3}')
umn = sm.symbols('u_m^n')
t, h, tau, ut, utt, uttt = sm.symbols("t h tau u^{'}_t u^{''}_t u^{'''}_t")
# h, tau = sm.symbols("h tau")
ux, uxx, uxxx = sm.symbols("u^{'}_x u^{''}_x u^{'''}_x")
u4x, u5x, u6x = sm.symbols("u^{(4)}_x u^{(5)}_x u^{(6)}_x")

u6x = 0
u5x = 0
#u4x = 0

ut = -(tau+6)*ux
utt = (tau+6)*uxx + ux
uttt = (tau+6)*uxxx + 2*uxx


un1 = umn + tau*ut + tau**2/2*utt + tau**3/6*uttt
um1 = umn + h*ux + h**2/2*uxx + h**3/6*uxxx + \
    h**4/24*u4x + h**5/120*u5x + h**6/720*u6x
um2 = umn + (2*h)*ux + (2*h)**2/2*uxx + (2*h)**3/6*uxxx + \
    (2*h)**4/24*u4x + (2*h)**5/120*u5x + (2*h)**6/720*u6x
um3 = umn + (3*h)*ux + (3*h)**2/2*uxx + (3*h)**3/6*uxxx +\
    (3*h)**4/24*u4x + (3*h)**5/120*u5x + (3*h)**6/720*u6x

exp = umn - un1 + tau/(6*h)*(tau/2 + t + 6)*(2*um3 - 9*um2 + 18*um1 - 11*umn) + \
    tau**2/(2*h**2) *(tau + t + 6)*(t+6) *(-um3 + 4*um2 - 5*um1 + 2*umn) - \
    tau**3/(6*h**3)*(t+6)**3 * (-um3 + 3*um2 - 3*um1 + umn)
    
exp = exp/tau

print('r:')
display(sm.simplify(exp))
