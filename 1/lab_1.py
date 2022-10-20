import math
import numpy as np
import decimal
from decimal import Decimal

print('Введите N:')
N = int(input())

pi = Decimal(math.pi)
e  = Decimal(math.e)

y = np.zeros(N)
y[0] = Decimal(pi/(e + Decimal(17)))
y[1] = Decimal(y[0])

for i in range(2, N):
    y[i] =(pi - Decimal(5)  * Decimal(y[i-2]) 
              - Decimal(12) * Decimal(y[i-1]))/e

np.savetxt('res.txt', y)
