import numpy as np
from scipy.signal import find_peaks

# =========================
# SETTINGS
# =========================
models_order = ["TIP3P", "TIP4P/2005", "TIP5P", "OPC"]

files = {
    "TIP3P": "gOO_tip3p_260K_NVT.dat",
    "TIP4P/2005": "gOO_tip4p2005_260K_NVT.dat",
    "TIP5P": "gOO_tip5p_260K_NVT.dat",
    "OPC": "gOO_opc_260K_NVT.dat"
}

# =========================
# FUNCTION: EXTRACT RDF FEATURES
# =========================
def extract_rdf_features(fname):
    data = np.loadtxt(fname, skiprows=1)
    r = data[:, 0]
    g = data[:, 1]

    peaks, _ = find_peaks(g, height=1.5)
    if len(peaks) == 0:
        return None

    p = peaks[0]
    r1 = r[p]
    g1 = g[p]

    # first minimum after peak
    after_peak = g[p+1:]
    m_local = np.argmin(after_peak)
    m = p + 1 + m_local

    r_min = r[m]
    g_min = g[m]

    return r1, g1, r_min, g_min


# =========================
# COLLECT RESULTS
# =========================
results = {}

for model in models_order:
    fname = files[model]
    vals = extract_rdf_features(fname)
    if vals:
        results[model] = vals

# =========================
# PRINT LATEX TABLE ROWS
# =========================
def print_row(name, index, unit=""):
    row = name
    for model in models_order:
        val = results[model][index]
        if unit:
            row += f" & {val:.3f} {unit}"
        else:
            row += f" & {val:.3f}"
    row += r" \\"
    print(row)

print("\n% ===== LaTeX rows for table =====\n")

print_row("Peak position", 0)
print_row("1st peak height", 1)
print_row("1st well position", 2)
print_row("1st well depth", 3)