import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parameters ----------------
E_number = 395
R1 = (E_number % 6) + 1
R2 = (E_number % 3) + 1
R3_val = (E_number % 4) + 1

omega_rpm = 2000       
omega = omega_rpm * 2 * np.pi / 60  

Bore = 0.08            
Stroke = 0.11          
Crank_radius = Stroke / 2
Conn_rod_len = 0.235   
Area_piston = np.pi * (Bore / 2)**2
Mass_recip = 1.5 + R1/10  

# ---------------- Read & Clean Data ----------------
data = pd.read_excel("data.xlsx")
data.columns = data.columns.str.strip()

try:
    p1 = data["R3  = 1"].values * 1e5
    p2 = data["R3 = 2"].values * 1e5
    p3 = data["R3= 3"].values * 1e5
    p4 = data["R3 = 4"].values * 1e5
    crank_angle = data["Crank_Angle"].values
except KeyError as e:
    print(f"Error: Could not find column {e} in Excel.")
    print(f"Available columns are: {list(data.columns)}")
    exit()

theta_rad = np.radians(crank_angle)

# ---------------- Functions ----------------
def inertia_force(theta):
    return Mass_recip * omega**2 * Crank_radius * (
        np.cos(theta) + np.cos(2 * theta) / (Conn_rod_len / Crank_radius)
    )

def torque(p, theta):
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    l_r_ratio = Conn_rod_len / Crank_radius
    term = np.sqrt(l_r_ratio**2 - sin_theta**2)
    force_net = (p * Area_piston) - inertia_force(theta)
    return force_net * Crank_radius * sin_theta * (1 + cos_theta / term)

# ---------------- Extend & Shift ----------------
def extend_and_shift(T, shift_deg):
    T_720 = np.tile(T[:-1], 2)
    shift_indices = int(shift_deg * (len(T)-1) / 360)
    return np.roll(T_720, shift_indices)

t1_base = np.array([torque(p, t) for p, t in zip(p1, theta_rad)])
t2_base = np.array([torque(p, t) for p, t in zip(p2, theta_rad)])
t3_base = np.array([torque(p, t) for p, t in zip(p3, theta_rad)])
t4_base = np.array([torque(p, t) for p, t in zip(p4, theta_rad)])

# Firing order: 1-3-4-2
T1 = extend_and_shift(t1_base, 0)
T3 = extend_and_shift(t3_base, 180)
T4 = extend_and_shift(t4_base, 360)
T2 = extend_and_shift(t2_base, 540)

Total_T = T1 + T2 + T3 + T4
angles_720 = np.linspace(0, 720, len(Total_T))

# ---------------- Mean Torque ----------------
mean_torque = np.mean(Total_T)

# ---------------- Work per cycle ----------------
theta_rad_720 = np.radians(angles_720)
work_total = np.trapz(Total_T, theta_rad_720)

# ---------------- Segment Area Calculation ----------------
diff = Total_T - mean_torque
crossings = np.where(np.diff(np.sign(diff)))[0]
segment_indices = np.concatenate(([0], crossings+1, [len(diff)-1]))

areas = []
cum_energy = 0
cum_list = []

for i in range(len(segment_indices)-1):
    start = segment_indices[i]
    end = segment_indices[i+1]

    seg_theta = theta_rad_720[start:end+1]
    seg_diff = diff[start:end+1]

    area = np.trapz(seg_diff, seg_theta)

    areas.append(area)
    cum_energy += area
    cum_list.append(cum_energy)

area_table = pd.DataFrame({
    "Segment": np.arange(1, len(areas)+1),
    "Start Angle (deg)": angles_720[segment_indices[:-1]],
    "End Angle (deg)": angles_720[segment_indices[1:]],
    "Area (J)": areas,
    "Cumulative Energy (J)": cum_list
})

# ---------------- Maximum Energy Fluctuation ----------------
E_max = max(cum_list)
E_min = min(cum_list)
delta_E = E_max - E_min

print("\n--- Energy Fluctuation Table ---")
print(area_table)

print(f"\nMean Torque = {mean_torque:.2f} Nm")
print(f"Work done per cycle = {work_total:.2f} J")
print(f"Maximum Energy Fluctuation ΔE = {delta_E:.2f} J")

# ---------------- Save to Excel ----------------
area_table.to_excel("Energy_Fluctuation_Table.xlsx", index=False)

summary = pd.DataFrame({
    "Quantity": ["Mean Torque", "Work per cycle", "Max Energy Fluctuation ΔE"],
    "Value": [mean_torque, work_total, delta_E]
})
summary.to_excel("Summary_Results.xlsx", index=False)

# ---------------- Plot ----------------
plt.figure(figsize=(12,6))

plt.plot(angles_720, Total_T, linewidth=2, label='Resultant Torque')
plt.axhline(mean_torque, linestyle='--', linewidth=2, label=f"Mean Torque = {mean_torque:.1f} Nm")

plt.fill_between(angles_720, Total_T, mean_torque,
                 where=(Total_T > mean_torque),
                 alpha=0.3, label='Above Mean')

plt.fill_between(angles_720, Total_T, mean_torque,
                 where=(Total_T < mean_torque),
                 alpha=0.3, label='Below Mean')

plt.title("Turning Moment Diagram & Energy Fluctuation")
plt.xlabel("Crank Angle (Degrees)")
plt.ylabel("Torque (Nm)")
plt.xlim(0,720)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

print("\nFiles created:")
print(" - Energy_Fluctuation_Table.xlsx")
print(" - Summary_Results.xlsx")
