import h5py
import pickle
import matplotlib.pyplot as plt
import numpy as np

file_path = "pyvisa_compressed_data.hdf5"

print("--- EKSTRAKCJA PARAMETRU PORZĄDKU Z POLA FRAMES ---")
with h5py.File(file_path, "r") as f:
    raw_bytes = f['data']['values'][0]
    byte_stream = raw_bytes.tobytes()

try:
    # Ładujemy główny słownik kroków
    main_dict = pickle.loads(byte_stream)[0]
    sorted_steps = sorted(main_dict.keys())
    
    cv_values = []
    steps_out = []
    
    print("--> Przetwarzam kroki i wyciągam Order Parameter...")
    
    for step in sorted_steps:
        obj = main_dict[step]
        
        # Dobieramy się do obiektu 'frames'
        frames_data = obj.frames
        
        # Próbujemy wyciągnąć wartość. Ponieważ 'frames' wyświetliło się jak tabela,
        # sprawdzamy czy to Pandas DataFrame, słownik, czy obiekt z kolumnami.
        try:
            if hasattr(frames_data, 'get'):
                # Jeśli zachowuje się jak słownik
                val = frames_data.get('Order Parameter')
            elif hasattr(frames_data, 'columns'):
                # Jeśli to Pandas DataFrame
                val = frames_data['Order Parameter'].values
            else:
                # Jeśli to obiekt PyRETIS, najprawdopodobniej można go indeksować tekstowo
                # lub wyciągnąć atrybut przez słownik wewnętrzny
                val = frames_data['Order Parameter']
            
            # Jeśli 'val' to cała seria/lista klatek dla danej ścieżki, bierzemy średnią, 
            # maksimum lub pierwszą wartość. Zazwyczaj dla punktu bierzemy wartość maksymalną lub pierwszą.
            if isinstance(val, (np.ndarray, list)) or hasattr(val, 'values'):
                # Spłaszczamy na wypadek obiektów Pandas i bierzemy np. średnią wartość z tej ścieżki
                # (Możesz zmienić np. na np.max(val) jeśli szukasz punktu zwrotnego)
                actual_value = np.mean(val)
            else:
                actual_value = float(val)
                
            cv_values.append(actual_value)
            steps_out.append(step)
            
        except Exception as inner_e:
            # Awaryjna ścieżka jeśli powyższe indeksowanie nie zadziała:
            # Próbujemy wyciągnąć surowy tekst/wartość z reprezentacji tekstowej lub atrybutów
            try:
                # PyRETIS trzyma czasami dane w słowniku pod dziwnymi kluczami
                if hasattr(frames_data, '__dict__'):
                    for k, v in frames_data.__dict__.items():
                        if 'order' in k.lower():
                            cv_values.append(v[0] if hasattr(v, '__len__') else v)
                            steps_out.append(step)
                            break
            except:
                pass

    # --- GENEROWANIE WYKRESU ---
    if len(cv_values) > 0:
        plt.figure(figsize=(12, 6))
        
        # Jeśli punktów jest bardzo dużo, rysujemy co np. 10, żeby wykres był czytelny
        stride = max(1, len(cv_values) // 2000)
        
        plt.plot(steps_out[::stride], cv_values[::stride], color='#d62728', alpha=0.8, linewidth=1.5, label='Order Parameter (Średnia ze ścieżki)')
        
        plt.title("Wykres Parametru Porządku z 16% symulacji PyVisA (Odkodowany z sukcesem!)")
        plt.xlabel("Krok symulacji (Step)")
        plt.ylabel("Wartość Order Parameter")
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.legend()
        
        output_png = "ostateczny_wykres_sukces.png"
        plt.savefig(output_png, dpi=300, bbox_inches='tight')
        print(f"\n[GENIALNIE!] Prawdziwy wykres został wygenerowany i zapisany jako: '{output_png}'")
    else:
        print("\n[INFO] Nie udało się dopasować sposobu wyciągania danych z obiektu 'frames'.")
        # Jeśli to się stanie, podejrzymy dokładny typ obiektu frames
        first_obj = main_dict[sorted_steps[0]].frames
        print(f"Dokładny typ obiektu 'frames' to: {type(first_obj)}")
        if hasattr(first_obj, '__dict__'):
            print("Jego wewnętrzne zmienne:", list(first_obj.__dict__.keys()))

except Exception as e:
    print(f"\n[BŁĄD]: {e}")
