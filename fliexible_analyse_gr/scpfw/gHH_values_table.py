import numpy as np
from scipy.signal import find_peaks

# =========================
# SETTINGS
# =========================
models_order = ["TIP3P"]

files = {
    "TIP3P": "rdf_HH.xvg",
}

# =========================
# FUNCTION: READ XVG
# =========================
def read_xvg(fname):
    r = []
    g = []

    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                r.append(float(parts[0]))
                g.append(float(parts[1]))

    r = np.array(r)
    g = np.array(g)

    # convert nm → Å
    r = r * 10.0

    return r, g


# =========================
# FUNCTION: EXTRACT HH PEAK
# =========================
def extract_hh_peak(fname):
    r, g = read_xvg(fname)

    # remove zero / near-zero region
    mask = r > 0.8
    r = r[mask]
    g = g[mask]

    # optional smoothing (helps a lot for HH)
    from scipy.ndimage import gaussian_filter1d
    g = gaussian_filter1d(g, sigma=2)

    # find peaks
    peaks, _ = find_peaks(g, height=0.2)

    if len(peaks) == 0:
        print(f"Warning: No peak found in {fname}")
        return None

    # first peak = smallest r
    p1 = peaks[np.argmin(r[peaks])]

    return r[p1]


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

print("\n% ===== LaTeX row for HH RDF (.xvg) =====\n")
print(row)