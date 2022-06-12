import numpy as np

'''

James Clooney 
November 21, 2021 

The 4th Order Runge-Kutta Method for initial value 
problems. 

       d/dt( y(t) ) = f(t, y),    y(0) = y0

where y can be a either a vector or a scalar.  

'''

# 4th order Runge-Kutta 
def runge_kutta4(funct, tspan, y0, N):

    # Step size
    h = (tspan[1] - tspan[0]) / N

    # t vector and initalized y 
    t = np.linspace(tspan[0], tspan[1], N + 1)
    y = np.zeros((len(y0), N + 1))

    # Sets initial conditions 
    y[:,0] = y0 
    
    # Main for loop that generates our solution 
    for i in range(N): 
        
        k1 = funct(t[i], y[:,i])    
        k2 = funct(t[i] + h/2, y[:,i] + h*(k1/2))
        k3 = funct(t[i] + h/2, y[:,i] + h*(k2/2))
        k4 = funct(t[i] + h, y[:,i] + h*k3) 

        y[:,i + 1] =  y[:,i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

    return t, y
