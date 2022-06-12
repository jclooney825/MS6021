from numpy import *
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

'''
James Clooney 
MS6021
Assignment 5 
November 21, 2021 

A numerical method for solving the following PDE

                u_t = u_xx  

where u(x, 0) = sin(pi x^2) and u(0,t) = u(1,t) = 0
using an explicit Euler's method. This is the form of
a 1D heat equation where the constant of proportionality 
is 1. 

'''

# Number of position steps 
N = 64

# Position step size 
h = 1 / N

# Time step in terms of position step size 
tau = (h**2) / 4

# Number of time steps 
K = int(round(1 / tau))

# Set initial x 
x = linspace(0,1,N) 

# Set initial U 
U = zeros((N,K))

# Place in Initial condition 
U0 = sin(pi*x**2)
U[:,0] = U0

# Iterates through each time step
for i in range(0, K-1):
    # Iterates through each position
    for j in range(1, N - 1):
        # Euler's method for Parabolic PDEs
        U[j,i + 1] =  U[j,i] + (tau/h**2)*(U[j + 1,i] - 2*U[j, i] + U[j - 1,i])


def plot_results():
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
    ax.set_xlabel('$x$', fontsize = 14)
    ax.set_ylabel('$u(x,t)$', fontsize = 14)
    ax.set_title(fr'One Dimensional Heat Equation with $N = {N}$', fontsize = 15)
    ax.grid(True)

    # Plot results
    plt.plot(x, U[:, 0])
    plt.plot(x, U[:, int(K/100)])
    plt.plot(x, U[:, int(K/20)])
    plt.plot(x,U[:, int(K/5)])
    plt.plot(x,U[:, K-1])
    plt.legend([r'$t = 0$', \
                fr'$t = {round(int(K/100)*tau, 3)}$',\
                fr'$t = {round(int(K/20)*tau, 3)}$',\
                fr'$t = {round(int(K/5)*tau, 3)}$',\
                fr'$t = {round(int(K-1)*tau, 3)}$'])
    plt.show()


# Animation function 
def animation():
    # Plot settings 
    plt.rcParams['text.usetex'] = True
    fig = plt.figure()

    # Plot settings 
    ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
    ax.set_xlabel('x', fontsize = 14)
    ax.set_ylabel('u(x,t)', fontsize = 14)
    ax.set_title(f'Temperature in a 1D Rod: N = {N}', fontsize = 15)
    ax.grid(True)

    # Intialized line 
    line, = ax.plot([], [], lw=3)

    # Initialize and animate functions 
    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        line.set_data(x, U[:,5*i])
        return line,

    # Animation function 
    anim = FuncAnimation(fig, animate, init_func=init,
                               frames = 500, interval=30, blit=True)

    # Save animation 
    anim.save('HeatEquationPDE.gif')


plot_results()