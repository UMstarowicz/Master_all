import numpy as np
import matplotlib.pyplot as plt

# Load data (skip comments and metadata)
data = np.loadtxt('energy.xvg', comments=['@', '#'])

time = data[:, 0]
temp = data[:, 1]

# Compute stats
mean_T = np.mean(temp)
std_T = np.std(temp)

# Skip first part (e.g., first 5 ps)
mask = time > 1

time_eq = time[mask]
temp_eq = temp[mask]

print("Equilibrated Mean:", np.mean(temp_eq))
print("Equilibrated Std:", np.std(temp_eq))

# Plot
plt.figure()
plt.plot(time, temp)
plt.axhline(mean_T, linestyle='--', label=f"Mean = {mean_T:.2f} K")
plt.xlabel("Time (ps)")
plt.ylabel("Temperature (K)")
plt.title("Stability of the TIP4P/2005f")
plt.legend()

# ---- SAVE FIGURE ----
plt.savefig("scp_fw_temperature.png", dpi=300, bbox_inches="tight")

print(f"Mean: {mean_T:.2f} K")
print(f"Std: {std_T:.2f} K")

plt.show()