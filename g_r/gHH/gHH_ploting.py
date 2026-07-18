import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

import numpy as np
from scipy.signal import find_peaks

# =========================
# SETTINGS
# =========================
models_order = ["260K", "277K", "295K"] #["TIP3P", "TIP4P/2005", "TIP5P", "OPC"]

files = {
    "260K": "gHH_opc_260K_NPT.dat",
    "277K": "gHH_opc_277K_NPT.dat",
    "295K": "gHH_opc_295K_NPT.dat"
    #"TIP3P": "gHH_tip3p_260K_NPT.dat",
    #"TIP4P/2005": "gHH_tip4p2005_260K_NPT.dat",
    #"TIP5P": "gHH_tip5p_260K_NPT.dat",
    #"OPC": "gHH_opc_260K_NPT.dat"
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
# PLOT RDFs
# =========================

plt.figure(figsize=(8,6))

for model in models_order:

    # load data
    data = np.loadtxt(files[model], skiprows=1)

    r = data[:,0]
    g = data[:,1]

    # optional smoothing
    g_smooth = gaussian_filter1d(g, sigma=1)

    # plot RDF
    plt.plot(
        r,
        g_smooth,
        linewidth=2,
        label=model
    )


# =========================
# FIGURE SETTINGS
# =========================

plt.xlabel(r"r [$\AA$]", fontsize=14)
plt.ylabel(r"$g_{HH}(r)$", fontsize=14)

plt.title(
    r"Hydrogen-Hydrogen RDF OPC (NPT)",
    fontsize=15
)

plt.xlim(1.5, 8)
plt.ylim(bottom=0)

plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()

# =========================
# SAVE FIGURE
# =========================

plt.savefig(
    "gHH_OPC_NPT_rigid.png",
    dpi=300,
    bbox_inches='tight'
)

plt.show()