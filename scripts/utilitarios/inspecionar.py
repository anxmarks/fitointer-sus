from pathlib import Path
import pandas as pd

# Caminhos a partir da raiz do projeto, para rodar de qualquer lugar
BASE_DIR = Path(__file__).resolve().parents[2]

plantas = pd.read_csv(BASE_DIR / "dados/manuais/fitodata - Página1.csv")
banco_1 = pd.read_csv(BASE_DIR / "dados/baixados/coconut_csv-11-2026.csv")
banco_2 = pd.read_csv(BASE_DIR / "dados/baixados/coconut_csv_lite-11-2026.csv")

print("PLANTAS:", plantas.columns.tolist())
print(plantas.head())
print()

print("BANCO 1:", banco_1.columns.tolist())
print(banco_1.head())
print()

print("BANCO 2:", banco_2.columns.tolist())
print(banco_2.head())
print("BANCO 2:", banco_2.columns.tolist())