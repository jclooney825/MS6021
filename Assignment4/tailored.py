import numpy as np 

'''

James Clooney 
MS6021
Assignment 5 
November 21, 2021 

A more advanced form of trapezoidal integration 
for functions of the form

            f(x) = x^alpha * g(x)

where g(x) is a function specified by the user on 
the interval [0,1]. 
  
    tailored(g, alpha, N):
    g: function handle of g(x) function
    alpha: power of x 
    N: number of subintervals 

    Example below: 
        f(x) = x^alpha (2 + 3x^3)

'''

# Number of intervals 
N = 1000

# Alpha value 
alpha = 1.2

# Function handle 
g = lambda x: 2 + 3*x**3

# Tailored method
def tailored(g, alpha, N):

    # Bounds 
    a = 0 
    b = 1   

    # x values from bounds a to b with N subintervals 
    x = np.linspace(a, b, N)

    # x_(i+1) and x_(i) values 
    x_i = x[0:N-1]
    x_iplus1 = x[1:N]
   
    # Coefficient of interpolated line: g^I(x) = Ax + B
    A = (g(x_iplus1) - g(x_i))/(x_iplus1 - x_i)
    B =  g(x_i) - A*x_i

    # Integrals over region x_(i) and x_(i+1) 
    first_integral = (A/(alpha + 2))*((x_iplus1)**(alpha + 2) - (x_i)**(alpha + 2))
    second_integral = (B/(alpha + 1))*((x_iplus1)**(alpha + 1) - (x_i)**(alpha + 1))

    # Sums the resulting values and returns them
    result = np.sum(first_integral + second_integral)

    return result

# Numerical solution 
I_numerical = tailored(g, alpha, N)

# Analytical solution 
I = lambda a: (2/(a + 1)) + (3/(a + 4))

# Error between numerical and analytic solution 
err = abs(I(alpha) - I_numerical)
print('Error value:')
print(err)