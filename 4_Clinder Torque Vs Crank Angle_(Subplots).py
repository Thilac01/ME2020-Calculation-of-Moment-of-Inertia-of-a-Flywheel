import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parameters ----------------
E_number = 395
R1 = (E_number % 6) + 1
R2 = (E_number % 3) + 1
R3 = (E_number % 4) + 1

omega_rpm = 2000       # crank speed in RPM
omega = omega_rpm * 2 * np.pi / 60  # rad/s

Bore = 0.08            # m
Stroke = 0.11          # m
Crank_radius = Stroke / 2
Conn_rod_len = 0.235   # m
Area_piston = np.pi * (Bore / 2)**2
Compression_ratio = 8
Mass_recip = 1.5 + R1/10  # kg

# ---------------- Read Data ----------------
data = pd.read_excel("data.xlsx")
crank_angle = data["Crank_Angle"]  # 0-360 deg

pressure_1 = data["R3  = 1"] * 1e5 
pressure_2 = data["R3 = 2"] * 1e5 
pressure_3 = data["R3= 3"] * 1e5 
pressure_4 = data["R3 = 4"] * 1e5 

theta_rad = np.radians(crank_angle)

# ---------------- Functions ----------------
def inertia_force(theta):
    return Mass_recip * omega**2 * Crank_radius * (
        np.cos(theta) + np.cos(2 * theta) / Compression_ratio
    )

def torque(p, theta):
    sin_theta = np.sin(theta)
    term = 2 * Crank_radius * np.sqrt(
        Conn_rod_len**2 - (Crank_radius * sin_theta)**2
    )
    T = (p * Area_piston - inertia_force(theta)) * Crank_radius * (
        sin_theta + (Crank_radius * sin_theta) / term
    )
    return T

# ---------------- Calculate Torque per Cylinder ----------------
TOR_1 = np.array([torque(p, t) for p, t in zip(pressure_1, theta_rad)])
TOR_2 = np.array([torque(p, t) for p, t in zip(pressure_2, theta_rad)])
TOR_3 = np.array([torque(p, t) for p, t in zip(pressure_3, theta_rad)])
TOR_4 = np.array([torque(p, t) for p, t in zip(pressure_4, theta_rad)])

# ---------------- Firing Order Shifts ----------------
firing_order_torques = [TOR_1, TOR_3, TOR_4, TOR_2]
firing_order_labels = ["Cylinder 1", "Cylinder 3", "Cylinder 4", "Cylinder 2"]
angle_shifts_deg = [0, 180, 360, 540]

n_points = len(crank_angle)
points_per_deg = n_points / 360

def shift_torque(T, shift_deg):
    shift_points = int(shift_deg * points_per_deg)
    return np.roll(T, shift_points)

# Apply shifts
shifted_torques = [shift_torque(T, s) for T, s in zip(firing_order_torques, angle_shifts_deg)]

# ---------------- Plot in Subplots ----------------
fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
for i, ax in enumerate(axs):
    ax.plot(crank_angle, shifted_torques[i], label=firing_order_labels[i], color='C'+str(i))
    ax.set_ylabel("Torque (Nm)")
    ax.set_title(f"{firing_order_labels[i]} Torque (Shifted)")
    ax.grid(True)
    ax.legend()

axs[-1].set_xlabel("Crank Angle (deg)")
plt.tight_layout()
plt.show()
