import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

# =========================
# LOAD RDF
# =========================
rdf = np.loadtxt("rdf_OO.xvg", comments=["@", "#"])
r = rdf[:, 0]
g = rdf[:, 1]

# smooth RDF
g_s = gaussian_filter1d(g, sigma=1)

# =========================
# FIND FIRST PEAK
# =========================
peaks, _ = find_peaks(g_s)

if len(peaks) == 0:
    raise ValueError("No RDF peak found")

p = peaks[0]

# =========================
# FIND FIRST MINIMUM AFTER PEAK
# =========================
r_after = r[p:]
g_after = g_s[p:]

minima, _ = find_peaks(-g_after)

if len(minima) == 0:
    raise ValueError("No minimum found after first peak")

m = minima[0] + p

r_cut = r[m]
g_min = g_s[m]

# =========================
# LOAD CN
# =========================
cn = np.loadtxt("cn_OO.xvg", comments=["@", "#"])
r_cn = cn[:, 0]
N = cn[:, 1]

# interpolate coordination number
coord_number = np.interp(r_cut, r_cn, N)

# =========================
# OUTPUT
# =========================
print(f"First-shell cutoff (Å): {r_cut:.3f}")
print(f"g(min): {g_min:.3f}")
print(f"Coordination number: {coord_number:.3f}")
