#!/usr/bin/env python
#Jason Chen Workshop H2O2 Rocket
import numpy as np
import matplotlib.pyplot as plt
import math


"""Input Parameters"""

Mach_e = 5.                      #Mach
payload_mass = 907              #kg
k = 1.3
R = 924                         #J/(kg*K)
T_t = 3839                      #K
P_inf = 1                       #atm
g_0 = 9.81                       #m/s^2
r_earth = 6378100               #m
r_char = 6940                   #m

payload_frac = 0.05


"""End of Input Parameters"""


"""Part 1: Design of Rocket Nozzle"""
print("""____Part 1: Design of Rocket Nozzle____""", '\n')

#Calculate the Chamber Pressure
#Convert from atm to kPa
P_0 = P_inf * 101.325 
mass_prop = (1-payload_frac) / payload_frac * payload_mass
total_mass = payload_mass + mass_prop
cp = (k*R) / (k-1)
Te_T_t = 1 / (1 + (k-1) / 2 * Mach_e**2)
P_chamber =  P_inf / (Te_T_t)**(cp/R)   #atm

#Calculate Throat Diameter in m
D_mk = (Te_T_t**(cp/R) * (1/Te_T_t)**0.5 * Mach_e)
D_max = (2 / (k+1))**((k+1) / (2*(k-1)))
D_hat = D_mk / D_max
A_throat = (total_mass*g_0 / 
        (k*Mach_e**2*P_inf/P_chamber*(1/D_hat)*P_chamber*101325))
D_throat = np.sqrt(A_throat / np.pi) * 2

#Calculate the Throat Area Ratio
E_throat = A_throat
A_Ratio = 1/D_hat

#Calculate Mass Flow Rate
mdot = np.sqrt(k/R) * D_mk * (P_chamber*101325)*\
        A_throat * A_Ratio / np.sqrt(T_t)

#Calculate Exit Temperature
T_e = Te_T_t * T_t
#Calculate Exit Velocity
v_exit = Mach_e * np.sqrt(k * R * T_e)

#Calculate Specific Impulse
Isp = v_exit / g_0

#Calculate Burn Time
impulse = Isp * mass_prop*g_0
burn_time = impulse / (mdot*v_exit)

#Print the variables
print("Chamber Pressure = ","%.6g" % P_chamber, "atm")
print("Throat Diameter = ","%.6g" % D_throat, "m")
print("Exit to Throat Area Ratio = ","%.6g" % A_Ratio)
print("Mass flow rate = ","%.6g" % mdot, "kg/s")
print("Exit Velocity = ","%.6g" % v_exit, "m/s")
print("Exit Temperature = ","%.6g" % T_e, "K")
print("Specific Impulse = ","%.6g" % Isp, "s")
print("Burn Time = ","%.6g" % burn_time, "s")

"""Part 2: Vertical Trajectory"""

#Function for the gravitational force on rocket using distance and mass
def gravityAccel(alt):
        gravity = g_0 / ((1 + alt/r_earth)**2)
        return gravity

#Function for the pressure of air at altitude
def airPressure(alt):
        P0_kPa = P_inf * 101.325
        pressure_air = P0_kPa * np.exp(-alt / r_char)
        return pressure_air

def getThrust(pressure):
        return mdot * v_exit +\
                (P_0 - pressure) * 1000 * A_throat * A_Ratio

#Calculate time step
num_steps = 10**6
time_step = burn_time / num_steps

#Marching Variables
time = [0]
v_r = [0]       #m/s
z_r = [0]       #m

#Non-Marching variables
m_total = [total_mass]
thrust = [mdot * v_exit]#N
a_r = [0]               #m*s^-2
g = [g_0]               #m*s^-2
p = [P_0]               #kPa

for i in range(num_steps):
        #time value
        time_now = time[i] + time_step
        time.append(time_now) 

        #Calculate velocity
        v_r_now = v_r[i] + a_r[i] * time_step
        v_r.append(v_r_now)

        #Calculate altitude
        z_r_now = z_r[i] + v_r[i] * time_step
        z_r.append(z_r_now)

        #Calculate change in mass
        m_total_now = m_total[i] - mdot * time_step
        m_total.append(m_total_now)

        #Air Pressure
        p_now = airPressure(z_r[i])
        p.append(p_now)

        #Calculate thrust and gravity
        g_now = gravityAccel(z_r[i])
        g.append(g_now)

        thrust_now = getThrust(p[i]) 
        thrust.append(thrust_now)

        #Calculate acceleration
        a_r_now = thrust[i] / m_total[i] - g[i]
        a_r.append(a_r_now)

#Generate plots for altitude, velocity, and height
burn_time_array = np.linspace(0, burn_time, num_steps + 1)

plt.figure(figsize=(12, 4))
plt.subplot(2,2,1)
plt.plot(burn_time_array, a_r)
plt.title("Acceleration vs Time")
plt.ylabel("Acceleration [m/s^2]")
plt.xlabel("Time [s]")

plt.subplot(2,2,2)
plt.plot(burn_time_array, v_r)
plt.title("Velocity vs Time")
plt.ylabel("Velocity [m/s]")
plt.xlabel("Time [s]")

plt.subplot(2,2,3)
plt.plot(burn_time_array, z_r)
plt.title("Height vs Time")
plt.ylabel("Height [m]")
plt.xlabel("Time [s]")        

plt.tight_layout()

#Calculate and print the desired output variables
a_g0_ratio = a_r[-1] / g_0
end_burn_vel = v_r[-1]
end_burn_alt = z_r[-1] / 1000   #km

#Escape velocity 
esc_vel = np.sqrt(2 * g_0 * r_earth / (1 + z_r[-1] / r_earth))
end_vel_esc_vel_ratio = end_burn_vel / esc_vel

#Maximum altitude
z_max = (z_r[-1]+(end_burn_vel**2)*(1+z_r[-1]/r_earth)/2/g_0) / 1000 #km

#Gravity ratio
g_g0_ratio = gravityAccel(z_max*1000) / g_0
print("""____Part 2: Vertical Trajectory____""", '\n')
print("Number of steps", num_steps)
print('a/g0 at end of burn =',"%.6g" % a_g0_ratio)
print('Speed at end of burn =',"%.6g" % end_burn_vel,'m/s')
print('Altitude at end of burn =',"%.6g" % end_burn_alt,'km')
print("V_final/V_escape = ","%.6g" % end_vel_esc_vel_ratio)

if(end_vel_esc_vel_ratio >= 1):
        hyper_vel = np.sqrt(end_burn_vel**2 - 2 * g_0 * r_earth /
                            (1 + z_r[-1] / r_earth))
        print("\nThe rocket escapes Earth.")
        print("The hyperbolic excess speed is:", "%.6g" % hyper_vel, 'm/s' )
else:
        print("\nThe rocket fails to escape Earth.")
        print("The maximum altitude reached is:", "%.6g" % z_max, 'km')
        print("g/g0 at final altitude =","%.6g" % g_g0_ratio)

plt.show()
