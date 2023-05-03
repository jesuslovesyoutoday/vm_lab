from sympy import symbols, simplify

umn, t, tau, h, ut, utt, uttt, ux, uxx, uxxx, uxxxx = symbols("umn t tau h ut utt uttt ux uxx uxxx uxxxx")

un1 = umn + tau*ut + tau**2/2*utt + tau**3/6*uttt
um1 = umn + h*ux + h**2/2*uxx + h**3/6*uxxx + h**4/24*uxxxx
um2 = umn + (2*h)*ux + (2*h)**2/2*uxx + (2*h)**3/6*uxxx + (2*h)**4/24*uxxxx
um3 = umn + (3*h)*ux + (3*h)**2/2*uxx + (3*h)**3/6*uxxx + (3*h)**4/24*uxxxx

exp = umn - un1 + tau/(6*h)*(tau/2 + t + 6)*(2*um3 - 9*um2 + 18*um1 - 11*umn) + \
    tau**2/(2*h**2) *(tau + t + 6)*(t+6) *(-um3 + 4*um2 - 5*um1 + 2*umn) - \
    tau**3/(6*h**3)*(t+6)**3 * (-um3 + 3*um2 - 3*um1 + umn)
print(simplify(exp))
