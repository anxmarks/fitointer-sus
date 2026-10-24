import pandas as pd

df = pd.read_csv("resultados/ativos_principais.csv")
df["nome_exibido"] = df["name"].fillna(df["iupac_name"]).astype(str).str[:50]

for planta, g in df.groupby("Nome popular"):
    print(f"\n{planta} ({len(g)})")
    print(g[["nome_exibido", "cas", "encontrado_por"]].to_string(index=False))