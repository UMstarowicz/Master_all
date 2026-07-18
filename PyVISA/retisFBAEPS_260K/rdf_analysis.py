import numpy as np
import matplotlib.pyplot as plt

def load_xvg(filename):
    data = np.loadtxt(filename, comments=['@', '#'])
    return data[:,0], data[:,1]

def analyze_rdf(r, g):
    # First peak
    peak_idx = np.argmax(g)
    r_peak = r[peak_idx]
    g_peak = g[peak_idx]

    # First minimum after peak
    g_after = g[peak_idx:]
    min_idx_rel = np.argmin(g_after)
    min_idx = peak_idx + min_idx_rel
    r_min = r[min_idx]
    g_min = g[min_idx]

    # Area under first peak
    area = np.trapz(g[:min_idx], r[:min_idx])

    return {
        "r_peak": r_peak,
        "g_peak": g_peak,
        "r_min": r_min,
        "g_min": g_min,
        "area": area
    }

def plot_multiple(rdf_files, labels):
    plt.figure()

    for f, label in zip(rdf_files, labels):
        r, g = load_xvg(f)
        plt.plot(r, g, label=label)

    plt.xlabel("r (nm)")
    plt.ylabel("g(r)")
    plt.legend()
    plt.title("Radial Distribution Function")
    plt.grid()
    plt.show()

# === RUN ANALYSIS ===

files = ["rdf_OO.xvg", "rdf_OH.xvg", "rdf_HH.xvg"]
labels = ["O-O", "O-H", "H-H"]

for f in files:
    r, g = load_xvg(f)
    results = analyze_rdf(r, g)

    print(f"\nResults for {f}:")
    for k, v in results.items():
        print(f"{k}: {v:.4f}")

plot_multiple(files, labels)
