import numpy as np

data = np.loadtxt("energy.xvg", comments=["@", "#"])
time = data[:,0]
temp = data[:,1]

print("Mean T:", np.mean(temp))
print("Std T:", np.std(temp))
