import MDAnalysis as mda
u = mda.Universe('top.psf','traj.dcd')
donors = u.select_atoms('resname TIP3 and name H1 H2')  # wodory
acceptors = u.select_atoms('resname TIP3 and name O')  # tlenu
# Utwórz słownik na bieżące stany H-bond: (d_i,a_j) -> obecność
hbond_history = {}  # mapowanie (H_idx, O_idx) -> list of 0/1
for ts in u.trajectory:
    # znajdź wszystkie pary spełniające kryteria (PBC)
    for H in donors:
        # znajdź akceptory w odległości <=0.35 nm
        close_O = acceptors.select_atoms(f'around 0.35 index {H.index}')
        for O in close_O:
            angle = calc_angle(H.position, H.bonded_atom.position, O.position)  # H-D-O
            if angle >= 150:  # kąt prosty 180-150=30
                pair = (H.index, O.index)
                hbond_history.setdefault(pair, []).append(1)
    # dla wszystkich par (H,O) które nie znalazły się powyżej, dopisz 0
    for pair, series in hbond_history.items():
        if len(series) < ts.frame+1:
            series.append(0)
# Oblicz C(tau) jako średnią produktów (sumy s_i*s_i przesunięte) itd.
