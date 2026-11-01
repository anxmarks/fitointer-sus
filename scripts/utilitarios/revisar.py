from pathlib import Path
import pandas as pd

# Caminhos a partir da raiz do projeto, para rodar de qualquer lugar
BASE_DIR = Path(__file__).resolve().parents[2]

df = pd.read_csv(BASE_DIR / "resultados/1_plantas/ativos_principais.csv")
df["nome_exibido"] = df["name"].fillna(df["iupac_name"]).astype(str).str[:50]

for planta, g in df.groupby("Nome popular"):
    print(f"\n{planta} ({len(g)})")
    print(g[["nome_exibido", "cas", "encontrado_por"]].to_string(index=False))