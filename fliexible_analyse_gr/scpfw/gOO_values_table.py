import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

# =========================
# SETTINGS
# =========================
models_order = ["TIP3P"]

files = {
    "TIP3P": "rdf_OO.xvg"
}

# =========================
# FUNCTION: READ XVG
# =========================
def read_xvg(fname):

    r = []
    g = []

    with open(fname, 'r') as f:

        for line in f:

            # skip metadata
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
# FUNCTION: EXTRACT OO FEATURES
# =========================
def extract_oo_features(fname):

    r, g = read_xvg(fname)

    # remove unphysical short-distance region
    mask = r > 1.5

    r = r[mask]
    g = g[mask]

    # smooth RDF slightly
    g_smooth = gaussian_filter1d(g, sigma=2)

    # =========================
    # FIND FIRST PEAK
    # =========================

    peaks, _ = find_peaks(g_smooth)

    if len(peaks) == 0:
        print("No peaks found.")
        return None

    first_peak_idx = peaks[0]

    r_peak = r[first_peak_idx]
    g_peak = g_smooth[first_peak_idx]

    # =========================
    # FIND FIRST MINIMUM
    # =========================

    # search only AFTER first peak
    r_after_peak = r[first_peak_idx:]
    g_after_peak = g_smooth[first_peak_idx:]

    # optional: restrict search window
    mask_min = (r_after_peak > 2.8) & (r_after_peak < 4.5)

    r_min_search = r_after_peak[mask_min]
    g_min_search = g_after_peak[mask_min]

    minima, _ = find_peaks(-g_min_search)

    if len(minima) == 0:
        print("No minima found.")
        return None

    first_min_idx = minima[0]

    r_min = r_min_search[first_min_idx]
    g_min = g_min_search[first_min_idx]

    # =========================
    # PRINT RESULTS
    # =========================

    print(f"\nResults for {fname}")
    print(f"First peak position : {r_peak:.3f} Å")
    print(f"First peak height   : {g_peak:.3f}")
    print(f"First minimum pos   : {r_min:.3f} Å")
    print(f"First minimum depth : {g_min:.3f}")

    return r_peak, g_peak, r_min, g_min


# =========================
# COLLECT RESULTS
# =========================
results = {}

for model in models_order:

    vals = extract_oo_features(files[model])

    if vals is not None:
        results[model] = vals


# =========================
# PRINT LATEX TABLE ROWS
# =========================

row1 = "Peak position"
for model in models_order:
    row1 += f" & {results[model][0]:.3f}"
row1 += r" \\"

row2 = "1st peak height"
for model in models_order:
    row2 += f" & {results[model][1]:.3f}"
row2 += r" \\"

row3 = "1st well position"
for model in models_order:
    row3 += f" & {results[model][2]:.3f}"
row3 += r" \\"

row4 = "1st well depth"
for model in models_order:
    row4 += f" & {results[model][3]:.3f}"
row4 += r" \\"

print("\n% ===== LaTeX rows =====\n")

print(row1)
print(row2)
print(row3)
print(row4)