import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# =========================
# SETTINGS
# =========================
models_order = ["TIP3P"]

files = {
    "TIP3P": "rdf_OO.xvg"
}
rdf_type = "OO"   # change to OO / OH / HH
temperature = "295K"
ensemble = "NPT"

# =========================
# FUNCTION: READ XVG
# =========================
def read_xvg(fname):

    r = []
    g = []

    with open(fname, 'r') as f:

        for line in f:

            # skip comments and metadata
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
# PLOT RDF
# =========================
plt.figure(figsize=(8,6))

for model in models_order:

    r, g = read_xvg(files[model])

    # optional smoothing
    g_smooth = gaussian_filter1d(g, sigma=1)

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
plt.ylabel(r"$g(r)$", fontsize=14)

plt.title(
    f"{rdf_type} RDF at {temperature} ({ensemble})",
    fontsize=15
)

plt.legend()
plt.grid(alpha=0.3)

plt.xlim(0, 8)

plt.tight_layout()

# =========================
# SAVE FIGURE
# =========================
filename = f"rdf_{rdf_type}_{temperature}_{ensemble}.png"

plt.savefig(
    filename,
    dpi=300,
    bbox_inches='tight'
)

plt.show()

print(f"\nFigure saved as: {filename}")