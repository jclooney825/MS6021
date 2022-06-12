from numpy import * 
import matplotlib.pyplot as plt 
from eulers_method import eulers_method
from predictor_corrector import predictor_corrector
from runge_kutta import runge_kutta4

# Number of subintervals 
N = 20

# Initial conditions 
y0 = array([1]) 

# Function handle f(t,y) 
funct = lambda t,y: array( [-y[0] - 5*exp(-t)*sin(5*t)] )

# Interval over which the ode is solved
tspan = array([0, 3]) 
  
# Solutions with both methods
[t,y] = eulers_method(funct, tspan, y0, N)
[t_pc, y_pc] = predictor_corrector(funct, tspan, y0, N)
[t_rk, y_rk] = runge_kutta4(funct, tspan, y0, N)

# Analytical solution 
t_exact = linspace(tspan[0], tspan[1], 1000)
y_exact = exp(-t_exact)*cos(5*t_exact)

# Plot values
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

fig = plt.figure()
ax = plt.axes()

plt.plot(t , y[0], label = 'Euler')
plt.plot(t_pc, y_pc[0], label = 'P/C', color = 'red')
plt.plot(t_rk, y_rk[0], label = 'RK4')
plt.plot(t_exact,y_exact, label = 'Analytical', color = 'black')
plt.legend()
ax.set_xlabel(r'$t$', fontsize = 12)
ax.set_ylabel(r'$y(t)$', fontsize = 12,rotation=0)
plt.grid(True)
ax.set_title(r'Various Numerical Methods for Solving ODEs',fontsize=15)
plt.show()