import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

# 1 a)
T = 50.052 # in years
a = 20 # in astronomical units
G = 6.67408*10**(-11) # gravitational mass in m**3 kg**(-1) s**(-2)
SUN_MASS = 1.99*10**30 # mass of the sun in kg

def AE2M(astro_units):
    return astro_units*1.49597870700*10**11

def AGE2SEC(ages):
    return ages*365.25*24*3.6*10**3

M = 4*np.pi**2/(G*AGE2SEC(T)**2)*AE2M(a)**3
RELATIVE_MASS = M/SUN_MASS
print "The total mass of the binary star system Sirius with respect to the mass of our sun is: " + str(RELATIVE_MASS) + "(in solar masses)"


###########

# 1 b)

a_SMALL = 4.3*3 # approximated semi major axis from the figure
MASS_RATIO = a_SMALL/a


SIRIUS1_MASS = RELATIVE_MASS*MASS_RATIO
SIRIUS2_MASS = RELATIVE_MASS*(1-MASS_RATIO)

SIRIUS1_MASS_STR = "{0:.4f}".format(SIRIUS1_MASS)
SIRIUS2_MASS_STR = "{0:.4f}".format(SIRIUS2_MASS)

print "The mass ratio of Sirius A to Sirius B equals: " + str(MASS_RATIO)
print "The mass of Sirius A is: " + SIRIUS1_MASS_STR + "(in solar masses)"
print "The mass of Sirius B is: " + SIRIUS2_MASS_STR + "(in solar masses)"

###########

mass1 = 1  
mass2 = 1
G     = 1
alpha = 0.5

def vanderpol(X, t):
    r11 = X[0]
    r12 = X[1]
    z11 = X[2]
    z12 = X[3]
    r21 = X[4]
    r22 = X[5]
    z21 = X[6]
    z22 = X[7]
    r  = np.sqrt((r11-r21)**2+(r12-r22)**2)
    dr11dt = z11
    dz11dt = -G*mass2/float(r**3)*(r11-r21)
    dr12dt = z12
    dz12dt = -G*mass2/float(r**3)*(r12-r22) 
    dr21dt = z21
    dz21dt = G*mass1/float(r**3)*(r11-r21)
    dr22dt = z22
    dz22dt = G*mass1/float(r**3)*(r12-r22)
    return [dr11dt, dr12dt, dz11dt, dz12dt, dr21dt, dr22dt, dz21dt, dz22dt]

X0 = [1, 0,0,alpha,-1, 0,0,-alpha]
t = np.linspace(0, 40, 250)

sol = odeint(vanderpol, X0, t)

r11 = sol[:, 0]
r12 = sol[:, 1]
z11 = sol[:, 2]
z12 = sol[:, 3]
r21 = sol[:, 4]
r22 = sol[:, 5]
z21 = sol[:, 6]
z22 = sol[:, 7]

plt.plot(t,r11, t, r12, t, r21, t, r22)
plt.xlabel('t')
plt.legend(('r11', 'r12', 'r21', r22))
plt.savefig('images/vanderpol-1.png')

# phase portrait
plt.figure()
plt.plot(r11,r12,r21,r22)
plt.plot(r11[0], r12[0], r21[0], r22[0], 'ro')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('images/vanderpol-2.png')
