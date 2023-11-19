#!/usr/bin/env python
#Jason Chen Workshop H2O2 Rocket
import numpy as np
import matplotlib.pyplot as plt


"""Input Parameters"""

Mach_e = 4          #Mach
payload_mass = 907  #kg
k = 1.3
R = 924             #J/(kg*K)
T_t = 3839          #K
P_inf = 1           #atm
g = 9.81            #m/s^2
r_earth = 6378100   #m
r_char = 6940       #m


"""Part 1: Design of Rocket Nozzle"""

#Calculate the Chamber Pressure
#Convert from atm to kPa
P_0 = P_inf * 101.325                   #kPa
massProp = (1-0.1) * 10 * payload_mass  #kg
total_mass = payload_mass * 10          #kg
cp = (k*R) / (k-1)
Te_T_t = 1 / (1 + (k-1) / 2 * Mach_e**2)
P_chamber =  P_inf / (Te_T_t)**(cp/R)   #atm
print("Chamber Pressure = ", P_chamber, "atm")

#Calculate Throat Diameter in m
D_mk = (Te_T_t**(cp/R) * (1/Te_T_t)**0.5 * Mach_e)
D_max = (2 / (k+1))**((k+1) / (2*(k-1)))
D_hat = D_mk / D_max
A_throat = (total_mass*g / 
        (k*Mach_e**2*P_inf/P_chamber*(1/D_hat)*P_chamber*101325))
D_throat = np.sqrt(A_throat / np.pi) * 2
print("Throat Diameter = ", D_throat, "m")

#Calculate the Throat Area Ratio
E_throat = A_throat
A_Ratio = 1/D_hat
print("Exit to Throat Area Ratio = ", A_Ratio)

#Calculate Mass Flow Rate
mdot = np.sqrt(k/R) * D_mk * (P_chamber*101325)*\
        A_throat * A_Ratio / np.sqrt(T_t)
print("Mass flow rate = ", mdot, "kg/s")

#Calculate Exit Temperature
T_e = Te_T_t * T_t
#Calculate Exit Velocity
E_vel = Mach_e * np.sqrt(k * R * T_e)
print("Exit Velocity = ", E_vel, "m/s")
print("Exit Temperature = ", T_e, "K")

#Calculate Specific Impulse
Isp = E_vel / g
print("Specific Impulse = ", Isp, "s")

#Calculate Burn Time
Impulse = Isp * massProp*g
burn_time = Impulse / (mdot*E_vel)
print("Burn Time = ", burn_time, "s")
