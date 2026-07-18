import numpy as np
from scipy.signal import find_peaks

# =========================
# SETTINGS
# =========================
models_order = ["TIP3P", "TIP4P/2005", "TIP5P", "OPC"]

files = {
    "TIP3P": "gOH_tip3p_295K_NPT.dat", 
    "TIP4P/2005": "gOH_tip4p2005_295K_NPT.dat", 
    "TIP5P": "gOH_tip5p_295K_NPT.dat", 
    "OPC": "gOH_opc_295K_NPT.dat" 
}

# =========================
# FUNCTION: EXTRACT OH PEAKS
# =========================
def extract_oh_peaks(fname):
    data = np.loadtxt(fname, skiprows=1)
    r = data[:, 0]
    g = data[:, 1]

    # ---- remove unphysical zero region ----
    mask = r > 1.2   # critical for OH RDF
    r = r[mask]
    g = g[mask]

    # ---- find peaks ----
    peaks, props = find_peaks(g, height=0.5, distance=10)

    if len(peaks) < 2:
        print(f"Warning: Less than 2 peaks found in {fname}")
        return None

    # sort peaks by position (not height!)
    peaks_sorted = peaks[np.argsort(r[peaks])]

    # first and second peak
    p1 = peaks_sorted[0]
    p2 = peaks_sorted[1]

    r1 = r[p1]
    r2 = r[p2]

    return r1, r2


# =========================
# COLLECT RESULTS
# =========================
results = {}

for model in models_order:
    vals = extract_oh_peaks(files[model])
    if vals:
        results[model] = vals

# =========================
# PRINT LATEX ROWS
# =========================
def print_row(name, index):
    row = name
    for model in models_order:
        val = results[model][index]
        row += f" & {val:.3f}"
    row += r" \\"
    print(row)

print("\n% ===== LaTeX rows for OH RDF =====\n")

print_row("1st peak position", 0)
print_row("2nd peak position", 1)