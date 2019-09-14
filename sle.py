import numpy as np
import matplotlib.pyplot as plt
import cmath
import math

# Numerical simulation of chordal SLE_k over the interval [0,T]
k = 4.0
T = 1.0

# Tolerence used by the variable step size control
tol = 0.0125

# Maximum and minimum step sizes allowed in the simulation
max_stepsize = T/2048.0
min_stepsize = T/(2.0**33.0)

sqrt_k = math.sqrt(k)

# Lists containing the discretized SLE trace and Brownian path.
sle_path = [0.0+0.0j]
brownian_path = []

# Computes the SLE trace associated with a constant driver.
# This is achieved by solving z' = - 2h/z at t = 0.5.
def horizontal_trace(z0, h):
    zt = cmath.sqrt(z0**2.0 - 2.0*h)
    
    if zt.imag < 0:
        zt = -zt
        
    return zt


# Computes the SLE trace associated with a vertical driver.
# This is achieved by solving z' = sqrt_k*increment at t = 1.
def vertical_trace(z0, increment):

    return z0 + sqrt_k*increment


# Computes the SLE trace using a step of the Ninomiya-Victoir scheme.
# This is equivalent to solving the backward Loewner equation driven
# by a piecewise linear path that has horizontal and vertical pieces.
def ninomiya_victoir(z0, h, brownian_increment):
    
    return horizontal_trace(vertical_trace(horizontal_trace(z0, h),
                                           brownian_increment), h)


# Propagates the numerical SLE with a variable step-size so that
#
# |SLE_{i+1} - SLE_{i}| < tol    for all neighbouring discretization points.
#
# This type of step size control was proposed for chordal SLE in the paper:
# Tom Kennedy, Numerical Computations for the Schramm-Loewner Evolution,
# Journal of Statistical Physics, 137:839, 2009.
def sle_step(time_increment, brownian_increment):
    l = len(brownian_path)
    
    # Approximate the next SLE point using the Ninomiya-Victoir scheme.
    # SLE traces are obtained by solving Loewner's differential equation
    # backwards in time starting from zero (hence the use of reversed()).
    zt = ninomiya_victoir(0.0+0.0j, time_increment, brownian_increment)
    
    for m in reversed(range(l)):
        zt = ninomiya_victoir(zt, brownian_path[m][0], brownian_path[m][1])
        
    # Check the propsed SLE point is sufficiently close to the previous point
    if ((abs(zt - sle_path[l-1]) < tol) or (time_increment <= min_stepsize))  \
    and (time_increment <= max_stepsize):
        
        # Update the numerical SLE trace and Brownian path        
        brownian_path.append([time_increment, brownian_increment])    
        sle_path.append(zt)

    else:
        # Generate the midpoint of the Brownian path over
        # the time interval using the Brownian bridge.        
        bridge_midpoint = 0.5*np.random.normal(0.0, math.sqrt(time_increment))

        half_time_increment = 0.5*time_increment
        half_increment = 0.5*brownian_increment
        
        # Approximate the SLE trace over the two subintervals
        zt = sle_step(half_time_increment, half_increment + bridge_midpoint)
        zt = sle_step(half_time_increment, half_increment - bridge_midpoint)
        
    return zt


# Compute the numerical SLE trace over the interval [0,T] using sle_step
zt = sle_step(T, np.random.normal(0.0, math.sqrt(T)));

# Plot the numerical SLE trace using matplotlib.pyplot
no_of_points = len(sle_path)

xvalues = []
yvalues = []

for i in range(no_of_points):
    xvalues.append(sle_path[i].real)
    yvalues.append(sle_path[i].imag)

plt.figure(0, dpi=300)
plt.plot(xvalues, yvalues, linewidth=0.5)
plt.title('SLE with k = ' + str(round(k, 2)))
plt.axis('scaled')
plt.show()

print('Number of steps = ' + str(no_of_points-1))