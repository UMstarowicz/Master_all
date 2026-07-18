import numpy as np

for ens_id, ens_name in [(0, '0-'), (1, '0+'), (2, '1+')]:
    fname = f'00{ens_id}/pathensemble.txt'
    total = 0
    crosses_next = 0
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            if len(cols) < 8:
                continue
            status = cols[7]
            if status == 'ACC':
                total += 1
                end = cols[5]
                if end == 'R':
                    crosses_next += 1
    prob = crosses_next / total if total > 0 else 0
    print(f'Ensemble [{ens_name}]: {total} ACC paths, {crosses_next} cross next interface, P = {prob:.4f}')
