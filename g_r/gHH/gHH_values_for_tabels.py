import numpy as np
from scipy.signal import find_peaks

# =========================
# SETTINGS
# =========================
models_order = ["TIP3P", "TIP4P/2005", "TIP5P", "OPC"]

files = {
    "TIP3P": "gHH_tip3p_260K_NVT.dat", 
    "TIP4P/2005": "gHH_tip4p2005_260K_NVT.dat", 
    "TIP5P": "gHH_tip5p_260K_NVT.dat", 
    "OPC": "gHH_opc_260K_NVT.dat" 
}

# =========================
# FUNCTION: FIRST HH PEAK
# =========================
def extract_hh_peak(fname):
    data = np.loadtxt(fname, skiprows=1)
    r = data[:, 0]
    g = data[:, 1]

    # ---- remove trivial zero region ----
    mask = r > 0.8   # HH peaks usually appear later than OH
    r = r[mask]
    g = g[mask]

    # ---- find peaks ----
    peaks, props = find_peaks(g, height=0.2)

    if len(peaks) == 0:
        print(f"Warning: No peak found in {fname}")
        return None

    # first peak = lowest r position
    p1 = peaks[np.argmin(r[peaks])]

    r1 = r[p1]

    return r1


# =========================
# COLLECT RESULTS
# =========================
results = {}

for model in models_order:
    val = extract_hh_peak(files[model])
    if val is not None:
        results[model] = val


# =========================
# PRINT LATEX ROW
# =========================
row = "1st peak position"
for model in models_order:
    row += f" & {results[model]:.3f}"
row += r" \\"

print("\n% ===== LaTeX row for HH RDF =====\n")
print(row)