import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parameters ----------------
E_number = 395
R1 = (E_number % 6) + 1
R2 = (E_number % 3) + 1
R3 = (E_number % 4) + 1  # not used for dataframe

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

# Clean column names (important)
df.columns = df.columns.str.strip()

print("Columns found in Excel:", df.columns)

# ---- CHANGE ONLY THIS LINE IF YOUR COLUMN NAME IS DIFFERENT ----
pressure_bar = df["R3 = 4"].to_numpy()     # <-- column name here
# ---------------------------------------------------------------

crank_angle_deg = df["Crank_Angle"].to_numpy()
theta_rad = np.radians(crank_angle_deg)

# Convert bar → Pa
pressure_Pa = pressure_bar * 1e5

# ---------------- Inertia Force Function ----------------
def inertia_force(theta):
    return Mass_recip * omega**2 * Crank_radius * (
        np.cos(theta) + np.cos(2 * theta) / Compression_ratio
    )

# ---------------- Force on Piston ----------------
Force_on_piston = Area_piston * pressure_Pa

# Print sample values
print("\nSample piston forces:")
for ang, F in zip(crank_angle_deg[:5], Force_on_piston[:5]):
    print(f"Crank angle = {ang:.1f}°, Force = {F:.2f} N")

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

# ---------------- Plot ----------------
plt.figure(figsize=(12, 6))
plt.plot(crank_angle_deg, TOR, lw=2, label="Torque")

# Zero torque line
plt.axhline(0, lw=1, label="Zero Torque")

# Stroke divisions
strokes = {
    "Intake": 0,
    "Compression": 180,
    "Power": 360,
    "Exhaust": 540,
    "Next Intake": 720
}

for name, angle in strokes.items():
    plt.axvline(angle, linestyle="--", lw=1)
    plt.text(angle + 5, max(TOR) * 0.8, name,
             rotation=90, va="center")

# ---------------- Vertical line at 570° ----------------
angle_570 = 570
torque_570 = np.interp(angle_570, crank_angle_deg, TOR)

plt.axvline(angle_570, linestyle="-.", lw=2,
            label=f"θ = 570°, T = {torque_570:.2f} N·m")
plt.plot(angle_570, torque_570, "o")

# Labels
plt.xlabel("Crank Angle (deg)")
plt.ylabel("Torque (N·m)")
plt.title("Torque vs Crank Angle")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
