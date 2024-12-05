import sympy as sm
from sympy.abc import x
import tabulate

# Parameters used
h = 0.001
max_iter = 200
tol = 1e-20 

# Inputs
funct = sm.sympify(input("Enter the function: "))
vrbl = sm.symbols(input("Enter the variable used: "))
guess = float(input("Initial Guess: "))
degree = sm.degree(funct)

# 1st Derivative
term1 = funct.subs(vrbl, vrbl-2*h)
term2 = funct.subs(vrbl, vrbl-h)
term3 = funct.subs(vrbl, vrbl+h)
term4 = funct.subs(vrbl, vrbl+2*h)
derv1 = (term1 -8*term2 + 8*term3 - term4)/(12*h)

# 2nd Derivative
derv2 = (-term1 + 16*term2 - 30*funct + 16*term3 - term4)/(12*h**2)

# Newton-Raphson on the 1st Derivative
i = 0
data = []
params = {"xprev": guess, "fxprev": 1, "xnew":1, "fxnew":1, "apprx_error":1}
params["xprev"] = guess - (derv1.subs(vrbl, guess)/derv2.subs(vrbl, guess))

if degree != 0 and degree != 1:
    while i <= max_iter and params["apprx_error"] > tol:
        i += 1
        params["xnew"] = params["xprev"] - (derv1.subs(vrbl, params["xprev"])/derv2.subs(vrbl, params["xprev"]))
        params["fxnew"] = derv1.subs(vrbl,  params["xnew"])
        if params["xnew"] == 0:
            params["apprx_error"] = 0
        else:
            params["apprx_error"] = abs((params["xnew"]-params["xprev"]))/abs(params["xnew"])
        params["xprev"] = params["xnew"]
        params["fxprev"] = derv1.subs(vrbl,  params["xprev"])

        data.append(params.copy().values())

    # Tests for Maxima or Minima
    if derv2.subs(vrbl, params["xnew"]) > 0:
        max_min = "Minima"
    elif derv2.subs(vrbl, params["xnew"]) < 0:
        max_min = "Maxima"

    # print the data in tabular form and final answer in 15 decimal digit precision
    print("\n" + tabulate.tabulate(data, params.keys(), "double_outline", showindex=True,  colalign = ("center", "center", "center" , "center", "center", "center")) + "\n")
    print(tabulate.tabulate([[params["xnew"], params["apprx_error"], max_min ]], ["Critical Value (x)", "Approximation Error", "Maxima/Minima"], "fancy_grid",floatfmt=[".15f"], colalign = ("center", "center", "center")), "\n")

else:
    print("The function has no critical value!")