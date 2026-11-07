import numpy as np
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Cruza plantas e fármacos da RENAME pela enzima (CYP3A4, CYP2C9, CYP2D6, P-gp):
#   composto da planta inibe a enzima + fármaco é substrato dela
#       -> aumento da exposição ao fármaco
#   composto da planta induz a enzima + fármaco é substrato dela
#       -> redução da exposição ao fármaco
# CYP1A2 e CYP2C19 ficam de fora: a tabela da FDA usada na RENAME só tem
# as quatro enzimas acima.
#
# Nível de evidência de cada par composto x fármaco:
#   alta  : ativo principal da planta, efeito forte ou moderado (até 10 uM),
#           medido no próprio composto, e fármaco sensível (FDA). Substratos
#           da P-gp não têm grau de sensibilidade na FDA e ficam no máximo em media
#   media : ativo principal com efeito forte ou moderado, ou qualquer
#           composto com efeito forte medido nele mesmo
#   baixa : o resto (efeito fraco, composto secundário ou dado de outro
#           estereoisômero)
# alta e media exigem também que os registros ativos sejam pelo menos tantos
# quanto os inativos (dados contraditórios ficam em baixa).
# Entrada : resultados/2_chembl/compostos_enzimas.csv       (padronizar_atividades.py)
#           resultados/2_chembl/atividades_padronizadas.csv (padronizar_atividades.py)
#           resultados/3_rename/resultados_RENAME_FDA/rename_substratos_FDA.csv
# Saída   : resultados/4_interacoes/interacoes_composto_farmaco.csv
#           resultados/4_interacoes/interacoes_planta_farmaco.csv
#           resultados/4_interacoes/matriz_planta_farmaco.csv
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]

ARQUIVO_COMPOSTOS_ENZIMAS = BASE_DIR / "resultados" / "2_chembl" / "compostos_enzimas.csv"
ARQUIVO_ATIVIDADES = BASE_DIR / "resultados" / "2_chembl" / "atividades_padronizadas.csv"
ARQUIVO_SUBSTRATOS = (
    BASE_DIR / "resultados" / "3_rename" / "resultados_RENAME_FDA" / "rename_substratos_FDA.csv"
)

PASTA_SAIDA = BASE_DIR / "resultados" / "4_interacoes"
PASTA_SAIDA.mkdir(parents=True, exist_ok=True)

MECANISMO = {
    "inibicao": "inibição da {enzima} -> aumento da exposição ao fármaco",
    "inducao": "indução da {enzima} -> redução da exposição ao fármaco",
}
ORDEM_NIVEL = {"alta": 3, "media": 2, "baixa": 1}


# ==========================================
# 2. LEITURA DOS ARQUIVOS
# ==========================================

comp_enz = pd.read_csv(ARQUIVO_COMPOSTOS_ENZIMAS, encoding="utf-8-sig")
ativ = pd.read_csv(ARQUIVO_ATIVIDADES, encoding="utf-8-sig", low_memory=False)
substratos = pd.read_csv(ARQUIVO_SUBSTRATOS, encoding="utf-8-sig")

# Efeitos da planta que mudam o metabolismo do fármaco, com evidência ativa
efeitos = comp_enz[
    comp_enz["efeito"].isin(MECANISMO) & (comp_enz["evidencia"] == "ativo")
].copy()

print("Pares composto x enzima com inibição ou indução ativa:", len(efeitos))
print("Fármacos da RENAME substratos (FDA):", substratos["farmaco"].nunique())


# ==========================================
# 3. REFERÊNCIAS DE CADA COMPOSTO x ENZIMA
# ==========================================

# Documentos do ChEMBL (ou texto da literatura) dos registros ativos; quando há
# dado do próprio composto, os de outros estereoisômeros ficam de fora, como
# no resumo de padronizar_atividades.py
chave = ["planta", "nome", "enzima", "efeito"]
ativos = ativ[ativ["resultado"] == "ativo"].copy()
tem_exato = ativos.groupby(chave)["correspondencia"].transform(
    lambda s: (s != "esqueleto").any()).astype(bool)
ativos = ativos[~tem_exato | (ativos["correspondencia"] != "esqueleto")]
referencias = ativos.groupby(chave, as_index=False).agg(
    referencias=("referencia", lambda s: "|".join(sorted(set(s.dropna().astype(str))))),
    ensaios=("assay_description", lambda s: " | ".join(sorted(set(s.dropna()))[:3])),
)
efeitos = efeitos.merge(referencias, on=chave, how="left")


# ==========================================
# 4. CRUZAMENTO PELA ENZIMA
# ==========================================

farmacos = substratos.rename(columns={
    "forca": "sensibilidade_farmaco", "observacao": "observacao_farmaco",
    "nomes_rename": "nomes_rename",
})[["farmaco", "nomes_rename", "enzima", "sensibilidade_farmaco", "observacao_farmaco"]]

inter = efeitos.merge(farmacos, on="enzima")
inter["mecanismo"] = [MECANISMO[e].format(enzima=z) for e, z in zip(inter["efeito"], inter["enzima"])]

proprio = inter["correspondencia"].str.split("|").apply(lambda c: c != ["esqueleto"])
efeito_relevante = inter["forca"].isin(["forte", "moderada"])
# Linhas da literatura podem vir sem ativo_principal preenchido (NaN viraria True)
principal = inter["ativo_principal"].fillna(False).astype(bool)
consistente = inter["ativos"] >= inter["inativos"]

inter["nivel_evidencia"] = np.select(
    [consistente & principal & efeito_relevante & proprio
     & (inter["sensibilidade_farmaco"] == "sensivel"),
     consistente & ((principal & efeito_relevante) | ((inter["forca"] == "forte") & proprio))],
    ["alta", "media"], "baixa",
)
inter["fonte_farmaco"] = "FDA (Drug Development and Drug Interactions, Table 1)"

colunas = [
    "planta", "nome", "inchikey", "ativo_principal", "enzima", "efeito", "forca", "melhor_pX",
    "ativos", "inativos", "correspondencia", "farmaco", "nomes_rename",
    "sensibilidade_farmaco", "observacao_farmaco", "mecanismo", "nivel_evidencia",
    "fontes", "referencias", "ensaios", "fonte_farmaco",
]
inter = inter[colunas].rename(columns={
    "nome": "composto", "forca": "forca_efeito", "fontes": "fonte_efeito",
    "ativos": "registros_ativos", "inativos": "registros_inativos",
})
inter["_nivel"] = inter["nivel_evidencia"].map(ORDEM_NIVEL)
inter = inter.sort_values(
    ["_nivel", "planta", "farmaco", "melhor_pX"], ascending=[False, True, True, False]
)


# ==========================================
# 5. RESUMO PLANTA x FÁRMACO
# ==========================================

def juntar(s):
    return "; ".join(dict.fromkeys(s.dropna().astype(str)))

resumo = inter.groupby(["planta", "farmaco"], as_index=False).agg(
    nomes_rename=("nomes_rename", "first"),
    nivel_evidencia=("_nivel", "max"),
    mecanismos=("mecanismo", juntar),
    enzimas=("enzima", lambda s: "|".join(sorted(set(s)))),
    compostos=("composto", "nunique"),
    ativos_principais=("composto", lambda s: juntar(
        s[inter.loc[s.index, "ativo_principal"].fillna(False).astype(bool)])),
    principais_compostos=("composto", lambda s: juntar(s.head(5))),
)
resumo["nivel_evidencia"] = resumo["nivel_evidencia"].map({v: k for k, v in ORDEM_NIVEL.items()})
resumo["_nivel"] = resumo["nivel_evidencia"].map(ORDEM_NIVEL)
resumo = resumo.sort_values(["_nivel", "planta", "farmaco"], ascending=[False, True, True]).drop(columns="_nivel")

matriz = resumo.pivot(index="planta", columns="farmaco", values="nivel_evidencia").fillna("")


# ==========================================
# 6. SALVAR RESULTADOS
# ==========================================

inter.drop(columns="_nivel").to_csv(
    PASTA_SAIDA / "interacoes_composto_farmaco.csv", index=False, encoding="utf-8-sig")
resumo.to_csv(PASTA_SAIDA / "interacoes_planta_farmaco.csv", index=False, encoding="utf-8-sig")
matriz.to_csv(PASTA_SAIDA / "matriz_planta_farmaco.csv", encoding="utf-8-sig")

print("\nPares composto x fármaco:", len(inter))
print("Pares planta x fármaco:", len(resumo))
print(resumo["nivel_evidencia"].value_counts().to_string())
print("\nPares planta x fármaco por enzima:")
print(inter.drop_duplicates(["planta", "farmaco", "enzima"]).groupby("enzima").size().to_string())
print("\nNível de evidência mais alto por planta:")
print(resumo.groupby("planta")["nivel_evidencia"].value_counts().unstack(fill_value=0).to_string())
print("\nArquivos salvos em:", PASTA_SAIDA)
