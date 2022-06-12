from numpy import *

'''
James Clooney 
MS6021
November 21, 2021 

The Predictor-Corrector Method for initial value 
problems. 

       d/dt(y(t)) = f(t, y),    y(0) = y0

where y can be a either a vector or a scalar. 

'''

# Predictor/Corrector Method 
def predictor_corrector(funct, tspan, y0, N):

    # Step size
    h = (tspan[1] - tspan[0]) / N

    # t vector and initalized y 
    t = linspace(tspan[0], tspan[1], N + 1)
    y = zeros((len(y0), N + 1))

    # Sets initial conditions 
    y[:,0] = y0 

    # Main foor loop that generates our solution 
    for i in range(N): 
        y_star = y[:,i] + multiply(h, funct(t[i], y[:,i]))
        y[:,i + 1] =  y[:,i] + multiply((h/2), add( funct(t[i], y[:,i]), funct(t[i + 1], y_star)) )

    return t, y 
