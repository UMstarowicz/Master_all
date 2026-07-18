import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("msd.xvg", comments=['@', '#'])

t = data[:,0]
msd = data[:,1]

# Linear fit (ignore first part → ballistic regime)
start = int(len(t)*0.2)   # skip first 20%
coef = np.polyfit(t[start:], msd[start:], 1)

D = coef[0] / 6  # diffusion coefficient

print(f"Diffusion coefficient: {D:.5e} nm^2/ps")

# Plot
plt.plot(t, msd, label="MSD")
plt.plot(t, coef[0]*t + coef[1], '--', label="fit")
plt.xlabel("Time (ps)")
plt.ylabel("MSD (nm^2)")
plt.legend()
plt.grid()
plt.show()

