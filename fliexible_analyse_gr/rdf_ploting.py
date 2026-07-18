import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# =====================================================
# SETTINGS
# =====================================================

# Main directory containing model folders
main_folder = r"C:\Users\ulast\Desktop\water_models_sims\fliexible_analyse_gr"

# RDF types to plot
rdf_types = ["OO", "OH", "HH"]

# Expected file names inside each model folder
rdf_files = {
    "OO": "rdf_OO.xvg",
    "OH": "rdf_OH.xvg",
    "HH": "rdf_HH.xvg"
}

# =====================================================
# FUNCTION: READ XVG
# =====================================================

def read_xvg(filepath):

    r = []
    g = []

    with open(filepath, "r") as f:

        for line in f:

            # Skip comments and metadata
            if line.startswith("#") or line.startswith("@"):
                continue

            parts = line.split()

            if len(parts) >= 2:
                r.append(float(parts[0]))
                g.append(float(parts[1]))

    r = np.array(r)
    g = np.array(g)

    # Convert nm -> Å
    r = r * 10.0

    return r, g

# =====================================================
# FIND MODEL FOLDERS
# =====================================================

model_folders = sorted([
    folder for folder in os.listdir(main_folder)
    if os.path.isdir(os.path.join(main_folder, folder))
])

print("\nDetected models:")
for model in model_folders:
    print(model)

# =====================================================
# PLOT EACH RDF TYPE
# =====================================================

for rdf in rdf_types:

    plt.figure(figsize=(8,6))

    for model in model_folders:

        filepath = os.path.join(
            main_folder,
            model,
            rdf_files[rdf]
        )

        # skip if file missing
        if not os.path.exists(filepath):
            print(f"Missing: {filepath}")
            continue

        # load RDF
        r, g = read_xvg(filepath)

        # smooth RDF slightly
        g_smooth = gaussian_filter1d(g, sigma=1)

        # plot
        plt.plot(
            r,
            g_smooth,
            linewidth=2,
            label=model
        )

    # =================================================
    # LABELS
    # =================================================

    plt.xlabel(r"r [$\AA$]", fontsize=14)

    if rdf == "OO":
        plt.ylabel(r"$g_{OO}(r)$", fontsize=14)
        title = "Oxygen-Oxygen RDF"

    elif rdf == "OH":
        plt.ylabel(r"$g_{OH}(r)$", fontsize=14)
        title = "Oxygen-Hydrogen RDF"

    elif rdf == "HH":
        plt.ylabel(r"$g_{HH}(r)$", fontsize=14)
        title = "Hydrogen-Hydrogen RDF"

    plt.title(title, fontsize=15)

    plt.xlim(1.5, 8)
    plt.ylim(bottom=0)

    plt.grid(alpha=0.3)
    plt.legend()

    plt.tight_layout()

    # =================================================
    # SAVE FIGURE
    # =================================================

    output_name = f"g{rdf}_comparison.png"

    plt.savefig(
        output_name,
        dpi=300,
        bbox_inches='tight'
    )

    print(f"\nSaved: {output_name}")

    plt.show()