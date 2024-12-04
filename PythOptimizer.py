import sympy as sm
from sympy.abc import x
import tabulate as tab

# Parameters used
h = 0.001
max_iter = 200
tol = 1e-20 

# Inputs
funct = sm.sympify(input("Enter the function: "))
vrbl = sm.symbols(input("Enter the variable used: "))
guess = float(input("Initial Guess: "))

# 1st Derivative
term1 = funct.subs(vrbl, vrbl-2*h)
term2 = funct.subs(vrbl, vrbl-h)
term3 = funct.subs(vrbl, vrbl+h)
term4 = funct.subs(vrbl, vrbl+2*h)
derv1 = (term1 -8*term2 + 8*term3 - term4)/(12*h)

# 2nd Derivative
derv2 = (-term1 + 16*term2 - 30*funct + 16*term3 - term4)/(12*h**2)

# Newton-Raphson
i = 0
apprx_error = 1
x0 = guess - (funct.subs(vrbl, guess)/derv1.subs(vrbl, guess))

while i <= max_iter or apprx_error != tol:
    i += 1
    xnew = x0 - (funct.subs(vrbl, x0)/derv1.subs(vrbl, x0))
    fxnew = derv1.subs(vrbl, xnew)
    x0 = xnew
