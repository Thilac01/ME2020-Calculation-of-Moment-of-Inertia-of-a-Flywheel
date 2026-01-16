import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parameters ----------------
E_number = 395
R1 = (E_number % 6) + 1

omega_rpm = 1660
omega = omega_rpm * 2 * np.pi / 60  # rad/s

Bore = 0.08                 # m
Stroke = 0.11               # m
Crank_radius = Stroke / 2
Conn_rod_len = 0.235        # m
Area_piston = np.pi * (Bore / 2)**2
Compression_ratio = 8
Mass_recip = 1.5 + R1 / 10  # kg

# ---------------- Read Data ----------------
df = pd.read_excel("data.xlsx")

# Clean column names
df.columns = df.columns.str.strip()

print("Columns in Excel:", df.columns)

# --------- CHANGE COLUMN NAME IF NEEDED ---------
crank_angle_deg = df["Crank_Angle"].to_numpy()
pressure_bar = df["R3 = 4"].to_numpy()
# -----------------------------------------------

theta_rad = np.radians(crank_angle_deg)
pressure_Pa = pressure_bar * 1e5  # bar → Pa

# ---------------- Inertia Force ----------------
def inertia_force(theta):
    return Mass_recip * omega**2 * Crank_radius * (
        np.cos(theta) + np.cos(2 * theta) / Compression_ratio
    )

# ---------------- Torque Function ----------------
def torque(p, theta):
    sin_theta = np.sin(theta)

    term = 2 * Crank_radius * np.sqrt(
        Conn_rod_len**2 - (Crank_radius * sin_theta)**2
    )

    T = (p * Area_piston - inertia_force(theta)) * Crank_radius * (
        sin_theta + (Crank_radius * sin_theta) / term
    )
    return T

# ---------------- Calculate Torque ----------------
TOR = np.array([torque(p, t) for p, t in zip(pressure_Pa, theta_rad)])

# ---------------- Work & Mean Torque ----------------
work_per_cycle = np.trapz(TOR, theta_rad)        # J
mean_torque = work_per_cycle / (4 * np.pi)       # N·m (4-stroke)

print(f"\nWork per cycle     : {work_per_cycle:.2f} J")
print(f"Mean torque        : {mean_torque:.2f} N·m")

# ---------------- Plot ----------------
plt.figure(figsize=(12, 6))
plt.plot(crank_angle_deg, TOR, lw=2, label="Instantaneous Torque")

# Zero & mean torque lines
plt.axhline(0, lw=1, label="Zero Torque")
plt.axhline(mean_torque, lw=2, linestyle='--',
            label=f"Mean Torque = {mean_torque:.2f} N·m")

# Stroke lines
strokes = {
    "Intake": 0,
    "Compression": 180,
    "Power": 360,
    "Exhaust": 540,
    "Next Intake": 720
}

for name, angle in strokes.items():
    plt.axvline(angle, linestyle='--', lw=1)
    plt.text(angle + 5, max(TOR) * 0.8, name,
             rotation=90, va='center')

plt.xlabel("Crank Angle (deg)")
plt.ylabel("Torque (N·m)")
plt.title("Torque vs Crank Angle")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
