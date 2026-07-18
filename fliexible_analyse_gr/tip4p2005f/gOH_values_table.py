import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

# =========================
# SETTINGS
# =========================
models_order = ["TIP3P"]

files = {
    "TIP3P": "rdf_OH.xvg",
}

# =========================
# FUNCTION: READ XVG
# =========================
def read_xvg(fname):
    r = []
    g = []

    with open(fname, 'r') as f:
        for line in f:
            # skip metadata/comments
            if line.startswith('#') or line.startswith('@'):
                continue

            parts = line.split()

            if len(parts) >= 2:
                r.append(float(parts[0]))
                g.append(float(parts[1]))

    r = np.array(r)
    g = np.array(g)

    # convert nm -> Å
    r = r * 10.0

    return r, g


# =========================
# FUNCTION: EXTRACT OH PEAKS
# =========================
def extract_oh_peaks(fname):

    r, g = read_xvg(fname)

    # remove initial zero region
    mask = r > 1.0
    r = r[mask]
    g = g[mask]

    # smooth data slightly
    g_smooth = gaussian_filter1d(g, sigma=2)

    # detect peaks
    peaks, properties = find_peaks(
        g_smooth,
        height=0.3,
        distance=20
    )

    if len(peaks) < 2:
        print(f"Warning: Less than 2 peaks found in {fname}")
        return None

    # sort by radial position
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

    if vals is not None:
        results[model] = vals


# =========================
# PRINT LATEX ROWS
# =========================
row1 = "1st peak position"
for model in models_order:
    row1 += f" & {results[model][0]:.3f}"
row1 += r" \\"

row2 = "2nd peak position"
for model in models_order:
    row2 += f" & {results[model][1]:.3f}"
row2 += r" \\"

print("\n% ===== LaTeX rows for OH RDF (.xvg) =====\n")

print(row1)
print(row2)