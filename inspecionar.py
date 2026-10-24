import pandas as pd

plantas = pd.read_csv("dados/fitodata - Página1.csv")
banco_1 = pd.read_csv("dados/coconut_csv-10-2026.csv")
banco_2 = pd.read_csv("dados/coconut_csv_lite-10-2026.csv")

print("PLANTAS:", plantas.columns.tolist())
print(plantas.head())
print()

print("BANCO 1:", banco_1.columns.tolist())
print(banco_1.head())
print()

print("BANCO 2:", banco_2.columns.tolist())
print(banco_2.head())
print("BANCO 2:", banco_2.columns.tolist())