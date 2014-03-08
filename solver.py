import numpy as np
import math
from copy import deepcopy
        
def to_linear_system(n):
        '''Returns (A, b, x), from where y can be found by 
        solving Ay=b and x is discretized interval'''

        h = 0.5 / (n + 1.0)
        A = np.zeros((n, n))
        x = np.zeros(n)

        for i in range(1, n + 1):
            x[i - 1] = 0.5 + (i - 0.5) * h

        for i in range(1, n - 1): 
                xi = x[i]
                ai = 2.0 * xi**3 - h * xi**2
                bi = -4.0 * xi**3 - 2.0 * h**2 * xi
                ci = 2.0 * xi**3 + h * xi**2

                A[i][i - 1] = ai
                A[i][i] = bi
                A[i][i + 1] = ci

        x1 = 0.5 + 0.5 * h
        a1 = 2.0 * x1**3 - h * x1**2
        b1 = -4.0 * x1**3 - 2.0 * h**2 * x1
        c1 = 2.0 * x1**3 + h * x1**2
        print(a1, b1, c1)
        A[0][0] = b1 - a1
        A[0][1] = c1

        xn = 0.5 + (n - 0.5) * h
        an = 2.0 * xn**3 - h * xn**2
        bn = -4.0 * xn**3 - 2.0 * h**2 * xn
        cn = 2.0 * xn**3 + h * xn**2
        A[n - 1][n - 2] = an
        A[n - 1][n - 1] = bn + cn

        b = np.zeros(n)
        for i in range(1, n):
            b[i - 1] =  -4.0 * h**2
        b[0] = -4.0 * h**2 - 2 * (-2.0 * math.log(2.0)) * a1
        
        return (A, b, x)
               
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    nf = len(a) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy the array
    for it in range(1, nf):
        mc = ac[it]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]
     
    xc = ac
    xc[-1] = dc[-1]/bc[-1]
      
    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
       
    del bc, cc, dc # delete variables from memory
        
    return xc

def TDMASolve(a, b, c, d):
    n = len(d) # n is the numbers of rows, a and c has length n-1
    for i in range(n-1):
        d[i+1] -= d[i] * a[i] / b[i]
        b[i+1] -= c[i] * a[i] / b[i]
    for i in reversed(range(n-1)):
        d[i] -= d[i+1] * c[i] / b[i+1]
    return [d[i] / b[i] for i in range(n)] # return the solution


def solve_system(A, b):
    return TDMASolve(np.diag(A, -1).copy(), np.diag(A, 0).copy(), np.diag(A, 1).copy(), b)
