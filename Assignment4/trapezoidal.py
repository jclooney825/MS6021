import numpy as np 

'''

James Clooney 
MS6021
Assignment 5 
November 21, 2021 

Trapezoidal integration method implemented in 
Python. 
  
    trapezoidal(funct, a, b, N):

    funct: function handle specified by user
    a: leftmost integral bound 
    b: rightmost integral bound 
    N: number of subintervals  

For more information, see the following link: 
https://en.wikipedia.org/wiki/Trapezoidal_rule

'''

# Number of intervals 
N = 10000

# Bounds of integration
a = 0 
b = 1

# Alpha value 
alpha = 1.2

# Function handle 
funct = lambda x: x**alpha*(2 + 3*x**3)


# Trapezoidal 
def trapezoidal(funct, a, b, N):

    # Step size 
    delta_x = (b - a) / N

    # x values from bounds a to b with N subintervals 
    x = np.linspace(a, b, N)

    # Generated values 
    vals = 2*funct(x) 

    # Multiply bound values by 1/2 
    vals[0] = (1/2)*(vals[0])
    vals[-1] = (1/2)*(vals[-1])

    # Sums the resulting values and returns them
    result = (delta_x/2)*np.sum(vals)
    
    return result 

# Numerical solution 
I_numerical = trapezoidal(funct, a, b, N)

# Analytical solution 
I = lambda alpha: (2/(alpha + 1)) + (3/(alpha + 4))

# Error between numerical and analytic solution 
err = abs(I(alpha) - I_numerical)
print('Error value:')
print(err)