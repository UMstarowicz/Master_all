import h5py
import matplotlib.pyplot as plt
import numpy as np

file_path = "pyvisa_compressed_data.hdf5"

with h5py.File(file_path, "r") as f:
    print("--- DEKODOWANIE STRUKTURY VALUES ---")
    
    # Wyciągamy główny kontener z danymi
    values_dataset = f['data']['values']
    
    # Ponieważ shape=(1,) i dtype=object, dobieramy się do pierwszego elementu
    inner_data = values_dataset[0]
    
    print(f"Typ danych wewnątrz 'values[0]': {type(inner_data)}")
    
    # Jeśli to tablica numpy, sprawdzamy jej wymiary
    if isinstance(inner_data, np.ndarray):
        print(f"Kształt (shape): {inner_data.shape}")
        print(f"Typ danych (dtype): {inner_data.dtype}")
        
        # Jeśli dane są zapisane jako surowe bajty (stringi/bajty z kompresji), 
        # spróbujemy podejrzeć pierwsze elementy
        print("\nPierwsze elementy (podgląd surowy):")
        print(inner_data[:5] if len(inner_data.shape) > 0 else inner_data)
        
        # PRÓBA WIZUALIZACJI:
        plt.figure(figsize=(10, 5))
        
        try:
            # Scenariusz A: Macierz 2D (np. wiersze = kroki, kolumny = interfejsy/zmienne)
            if len(inner_data.shape) == 2:
                for i in range(min(5, inner_data.shape[1])): # narysuj max 5 kolumn
                    plt.plot(inner_data[:, i], alpha=0.7, label=f"Zmienna {i}")
            
            # Scenariusz B: Zwykły wektor 1D
            elif len(inner_data.shape) == 1:
                # Jeśli to są liczby, rysujemy bezpośrednio
                if np.issubdtype(inner_data.dtype, np.number):
                    plt.plot(inner_data, alpha=0.7, label="Dane CV")
                else:
                    # Czasami w dtype=object kryją się kolejne tablice (tablica tablic)
                    print("\nWykryto tablicę obiektów (prawdopodobnie trajektorie dla różnych ścieżek):")
                    for idx, sub_arr in enumerate(inner_data[:5]):
                        if isinstance(sub_arr, (np.ndarray, list)):
                            print(f"  Ścieżka {idx}: kształt = {np.shape(sub_arr)}")
                            plt.plot(sub_arr, alpha=0.6, label=f"Ścieżka {idx}")
            
            plt.title("Wykres CV z odkodowanego 'values[0]' (16% symulacji)")
            plt.xlabel("Krok")
            plt.ylabel("Wartość")
            plt.legend()
            plt.grid(True)
            plt.savefig("wykres_odkodowany.png", dpi=300)
            print("\n[SUKCES] Wykres został zapisany jako 'wykres_odkodowany.png'!")
            
        except Exception as e:
            print(f"\n[BŁĄD RYSOWANIA]: {e}")
            
    else:
        print("Dane nie są tablicą numpy. Nie można ich narysować w ten sposób.")
