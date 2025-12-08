import numpy as np
import pandas as pd
from numpy.linalg import eig

# Utilities
def load_matrix(path):
    return pd.read_csv(path, header=None).values

def analyze_mode(pole):
    lam = pole
    wn = np.sqrt(lam.real**2 + lam.imag**2)
    if wn == 0:
        return 1.0, 0.0, 0.0
    return wn

def print_modes(title, A):
    print("\n" + "="*70)
    print(f"  {title}")
    print("="*70)

    eigvals, eigvecs = eig(A)

    modes = []
    for lam in eigvals:
        wn = analyze_mode(lam)
        modes.append((lam, wn))

    df = pd.DataFrame({
        "Pole (lambda)": [m[0] for m in modes],
        "omega_n (rad/s)": [m[1] for m in modes],
    })

    print(df.to_string(index=False))
    print()
    return df

def main():

    print("Loading A, B matrices...\n")

    # Rudderless
    A_rl = load_matrix("C:/Users/xuz-t/OneDrive/Documents/Horizontal_Rotary_Tail/build/output/A_rudderless.csv")

    # Ruddered
    A_r = load_matrix("C:/Users/xuz-t/OneDrive/Documents/Horizontal_Rotary_Tail/build/output/A_ruddered.csv")

    # Full modal analysis
    df_rl = print_modes("MODAL ANALYSIS: RUDDERLESS AIRCRAFT", A_rl)
    df_r  = print_modes("MODAL ANALYSIS: RUDDERED AIRCRAFT",   A_r)

    # Save mode tables
    df_rl.to_csv("C:/Users/xuz-t/OneDrive/Documents/Horizontal_Rotary_Tail/build/output/modes_rudderless.csv", index=False)
    df_r.to_csv("C:/Users/xuz-t/OneDrive/Documents/Horizontal_Rotary_Tail/build/output/modes_ruddered.csv", index=False)

    print("Saved modal analysis tables to:")
    print("  output/modes_rudderless.csv")
    print("  output/modes_ruddered.csv")

if __name__ == "__main__":
    main()