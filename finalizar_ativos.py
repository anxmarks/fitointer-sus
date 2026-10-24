import re
import pandas as pd

# ------------------------------------------------------------
# Limpa a tabela de ativos: tira moléculas que não são o ativo,
# preenche nomes vazios e junta duplicatas.
# Entrada : resultados/ativos_principais.csv
# Saída   : resultados/ativos_principais_final.csv
# ------------------------------------------------------------

df = pd.read_csv("resultados/ativos_principais.csv", low_memory=False)
antes = len(df)

# 1) Nome do composto: usa "name"; se vazio, tenta o CAS; senão o nome IUPAC
cas_nome = {                      # confirme no PubChem
    "458-37-7": "curcumin",
    "22608-11-3": "demethoxycurcumin",
    "21593-77-1": "S-allyl-L-cysteine",
}
def nome_composto(r):
    if pd.notna(r["name"]):
        return r["name"]
    for c in str(r.get("cas", "")).split("|"):
        if c in cas_nome:
            return cas_nome[c]
    return r["iupac_name"] if pd.notna(r["iupac_name"]) else r["identifier"]

df["composto"] = df.apply(nome_composto, axis=1)

# 2) Excluir o que bateu na busca mas não é o ativo principal
excluir = ["terpecurcumin", "bisabolocurcumin"]     # híbridos sesquiterpeno-curcumina
# excluir += ["acemannan"]                           # polissacarídeo grande (decisão de vocês)
mascara = df["composto"].str.lower().apply(lambda n: any(e in n for e in excluir))
print("Excluídos:", int(mascara.sum()))
df = df[~mascara]

# 3) Duplicatas: mesma estrutura na mesma planta
def chave_nome(n):
    n = str(n).lower().strip()
    return re.sub(r"^(\([+\-\u00b1rsRSEZ0-9,]*\)-|\[\d+\]-|[ld]-)", "", n)

if "standard_inchi_key" in df.columns:
    df["_dup"] = df["standard_inchi_key"].str[:14].fillna(df["composto"].apply(chave_nome))
else:
    df["_dup"] = df["composto"].apply(chave_nome)   # sem InChIKey: compara pelo nome

df = df.drop_duplicates(subset=["Nome popular", "_dup"]).drop(columns="_dup")

cols = [c for c in ["Nome popular", "Nome cientifico", "composto", "cas", "identifier",
                    "canonical_smiles", "molecular_formula", "encontrado_por"] if c in df.columns]
df[cols].to_csv("resultados/ativos_principais_final.csv", index=False)

print(f"{antes} -> {len(df)} moléculas")
print(df.groupby("Nome popular").size().to_string())