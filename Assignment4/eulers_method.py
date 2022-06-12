import numpy as np 

##############################################
'''

James Clooney 
MS6021
November 21, 2021 

The Explicit Euler's Method for initial value 
problems. 

       d/dt( y(t) ) = f(t, y),    y(0) = y0

where y can be a either a vector or a scalar.  

'''

##############################################

# Eulers method 
def eulers_method(funct, tspan, y0, N):

    # Step size
    h = (tspan[1] - tspan[0]) / N

    # t vector and initalized y 
    t = np.linspace(tspan[0], tspan[1], N + 1)
    y = np.zeros((len(y0), N + 1))

    # Sets initial conditions 
    y[:,0] = y0 
    
    # Main for loop that generates our solution 
    for i in range(N): 
        y[:,i + 1] =  y[:,i] + np.multiply(h, funct(t[i], y[:,i]))
   
    return t, y 
