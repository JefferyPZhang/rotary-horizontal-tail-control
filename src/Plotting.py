# ChatGPT-Generated Python Plotting Program

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data
rudderless = pd.read_csv("C:/Users/xuz-t/OneDrive/Documents/Horizontal_Rotary_Tail/build/output/rudderless_sim_history.csv", header=None)
ruddered   = pd.read_csv("C:/Users/xuz-t/OneDrive/Documents/Horizontal_Rotary_Tail/build/output/ruddered_sim_history.csv", header=None)

# Extract time vector
dt = 0.01
N = len(rudderless)
t = np.linspace(0, dt*(N-1), N)

# Extract states
# Navigation coordinates
PN_rl, PE_rl, PD_rl = rudderless[9],  rudderless[10],  rudderless[11]
PN_r, PE_r, PD_r  = ruddered[9],    ruddered[10],    ruddered[11]

# Attitude angles
phi_rl, theta_rl, psi_rl = rudderless[6], rudderless[7], rudderless[8]
phi_r, theta_r, psi_r  = ruddered[6],   ruddered[7],   ruddered[8]

# Extract inputs
# For rudderless: [da, de, tail, throttle]
da_rl, de_rl, tail_rl, dt_rl = rudderless[12], rudderless[13], rudderless[14], rudderless[15]
# For ruddered: [da, de, rudder, throttle]
da_r, de_r, dr_r, dt_r = ruddered[12], ruddered[13], ruddered[14], ruddered[15]

plt.style.use("seaborn-v0_8-darkgrid")

# 1. NAVIGATION COORDINATES P_N, P_E, P_D
plt.figure(figsize=(12,8))
plt.suptitle("Navigation Coordinate Comparison (Instability Dynamics)")

plt.subplot(3,1,1)
plt.plot(t, PN_rl, label="Rudderless")
plt.plot(t, PN_r, label="Ruddered", linestyle="--")
plt.ylabel("P_N [m]")
plt.legend()

plt.subplot(3,1,2)
plt.plot(t, PE_rl)
plt.plot(t, PE_r, linestyle="--")
plt.ylabel("P_E [m]")

plt.subplot(3,1,3)
plt.plot(t, -PD_rl, label="Rudderless")  
plt.plot(t, -PD_r, label="Ruddered", linestyle="--")
plt.ylabel("Altitude [m]")
plt.xlabel("Time [s]")

plt.tight_layout()
plt.subplots_adjust(top=0.94)
plt.show()

# 2. ATTITUDE ANGLES phi, theta, psi
plt.figure(figsize=(12,8))
plt.suptitle("Attitude Angles Comparison (Instability Dynamics)")

plt.subplot(3,1,1)
plt.plot(t, np.rad2deg(phi_rl), label="Rudderless")
plt.plot(t, np.rad2deg(phi_r), linestyle="--", label="Ruddered")
plt.ylabel("phi [deg]")
plt.legend()

plt.subplot(3,1,2)
plt.plot(t, np.rad2deg(theta_rl))
plt.plot(t, np.rad2deg(theta_r), linestyle="--")
plt.ylabel("theta [deg]")

plt.subplot(3,1,3)
plt.plot(t, np.rad2deg(psi_rl))
plt.plot(t, np.rad2deg(psi_r), linestyle="--")
plt.ylabel("psi [deg]")
plt.xlabel("Time [s]")

plt.tight_layout()
plt.subplots_adjust(top=0.94)
plt.show()

# 3. CONTROL INPUTS
plt.figure(figsize=(12,10))
plt.suptitle("Feedforward + Feedback Input Response Comparison")

# Aileron
plt.subplot(4,1,1)
plt.plot(t, np.rad2deg(da_rl), label="Rudderless")
plt.plot(t, np.rad2deg(da_r), linestyle="--", label="Ruddered")
plt.ylabel("delta_a [deg]")
plt.legend()

# Elevator
plt.subplot(4,1,2)
plt.plot(t, np.rad2deg(de_rl))
plt.plot(t, np.rad2deg(de_r), linestyle="--")
plt.ylabel("delta_e [deg]")

# Rudder or Tail
plt.subplot(4,1,3)
plt.plot(t, np.rad2deg(tail_rl), label="Tail (Rudderless)")
plt.plot(t, np.rad2deg(dr_r), linestyle="--", label="Rudder (Ruddered)")
plt.ylabel("delta_tail or delta_r [deg]")
plt.legend()

# Throttle
plt.subplot(4,1,4)
plt.plot(t, dt_rl)
plt.plot(t, dt_r, linestyle="--")
plt.ylabel("Throttle")
plt.xlabel("Time [s]")

plt.tight_layout()
plt.subplots_adjust(top=0.93)
plt.show()