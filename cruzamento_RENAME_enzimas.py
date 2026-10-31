from pathlib import Path
import pandas as pd

# ------------------------------------------------------------
# Liga os medicamentos da RENAME (já cruzados com o ClinPGx) às
# enzimas estudadas, usando as relações gene–fármaco do ClinPGx.
# Entrada : resultados/medicamentos_RENAME_consolidados_DDI.csv
#           dados/drugs.tsv
#           dados/relationships.tsv  (ClinPGx > Downloads > relationships.zip)
# Saída   : resultados/resultados_RENAME_enzimas/
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS E PARÂMETROS
# ==========================================

BASE_DIR = Path(__file__).resolve().parent

ARQUIVO_RENAME_DDI = BASE_DIR / "resultados" / "medicamentos_RENAME_consolidados_DDI.csv"
ARQUIVO_DRUGS = BASE_DIR / "dados" / "drugs.tsv"
ARQUIVO_RELACOES = BASE_DIR / "dados" / "relationships.tsv"

PASTA_SAIDA = BASE_DIR / "resultados" / "resultados_RENAME_enzimas"
PASTA_SAIDA.mkdir(parents=True, exist_ok=True)

# Gene no ClinPGx -> nome usado no projeto
ENZIMAS = {
    "CYP3A4": "CYP3A4",
    "CYP2C9": "CYP2C9",
    "CYP2D6": "CYP2D6",
    "ABCB1": "P-gp",
}

# "not associated" fica de fora; "ambiguous" entra, marcado na tabela
ASSOCIACOES_ACEITAS = ["associated", "ambiguous"]


# ==========================================
# 2. LEITURA DOS ARQUIVOS
# ==========================================

if not ARQUIVO_RELACOES.exists():
    raise SystemExit(
        f"Arquivo não encontrado: {ARQUIVO_RELACOES}\n"
        "Baixe relationships.zip na página de downloads do ClinPGx "
        "e salve o relationships.tsv dentro da pasta dados/."
    )

rename = pd.read_csv(ARQUIVO_RENAME_DDI, encoding="utf-8-sig", dtype=str)
drugs = pd.read_csv(ARQUIVO_DRUGS, sep="\t", encoding="utf-8-sig", dtype=str)
relacoes = pd.read_csv(ARQUIVO_RELACOES, sep="\t", encoding="utf-8-sig", dtype=str)

print("RENAME cruzada com ClinPGx:", len(rename))
print("Relações gene–fármaco:", len(relacoes))


# ==========================================
# 3. CÓDIGO DO CLINPGX PARA CADA MEDICAMENTO
# ==========================================

# ddi_nome é exatamente o "Name" do drugs.tsv
codigo_por_nome = dict(zip(drugs["Name"].str.lower(), drugs["PharmGKB Accession Id"]))
rename["pharmgkb_id"] = rename["ddi_nome"].str.lower().map(codigo_por_nome)

sem_codigo = rename["pharmgkb_id"].isna().sum()
print("Medicamentos sem código do ClinPGx:", sem_codigo)


# ==========================================
# 4. RELAÇÕES COM AS ENZIMAS ESTUDADAS
# ==========================================

# O arquivo traz as duas direções (Gene–Chemical e Chemical–Gene):
# deixa sempre o gene na coluna 1 e o fármaco na coluna 2
colunas_1 = ["Entity1_id", "Entity1_name", "Entity1_type"]
colunas_2 = ["Entity2_id", "Entity2_name", "Entity2_type"]

invertidas = relacoes["Entity1_type"].eq("Chemical") & relacoes["Entity2_type"].eq("Gene")
relacoes.loc[invertidas, colunas_1 + colunas_2] = (
    relacoes.loc[invertidas, colunas_2 + colunas_1].to_numpy()
)

relacoes = relacoes[
    relacoes["Entity1_type"].eq("Gene")
    & relacoes["Entity2_type"].eq("Chemical")
    & relacoes["Entity1_name"].isin(ENZIMAS)
    & relacoes["Association"].isin(ASSOCIACOES_ACEITAS)
].drop_duplicates(["Entity1_name", "Entity2_id"])

relacoes = relacoes.rename(columns={
    "Entity1_name": "gene",
    "Entity2_id": "pharmgkb_id",
    "Association": "associacao_clinpgx",
    "Evidence": "evidencia_clinpgx",
    "PK": "farmacocinetica",
    "PD": "farmacodinamica",
    "PMIDs": "pmids",
})
relacoes["enzima"] = relacoes["gene"].map(ENZIMAS)

print("Relações com as enzimas estudadas:", len(relacoes))


# ==========================================
# 5. CRUZAMENTO RENAME x ENZIMAS
# ==========================================

colunas_relacao = [
    "pharmgkb_id", "enzima", "gene", "associacao_clinpgx",
    "farmacocinetica", "farmacodinamica", "evidencia_clinpgx", "pmids",
]

cruzado = rename.merge(relacoes[colunas_relacao], on="pharmgkb_id", how="inner")

# Uma linha por fármaco x enzima; os nomes da RENAME (sais, formas) ficam juntos
tabela = (
    cruzado
    .groupby(["ddi_nome", "pharmgkb_id", "enzima"], as_index=False)
    .agg(
        nomes_rename=("nome_original", lambda s: "; ".join(sorted(set(s)))),
        origem_correspondencia=("origem_correspondencia", "first"),
        validacao_manual=("validacao_manual", "first"),
        gene=("gene", "first"),
        associacao_clinpgx=("associacao_clinpgx", "first"),
        farmacocinetica=("farmacocinetica", "first"),
        farmacodinamica=("farmacodinamica", "first"),
        evidencia_clinpgx=("evidencia_clinpgx", "first"),
        pmids=("pmids", "first"),
    )
    .sort_values(["enzima", "ddi_nome"])
)
tabela["fonte"] = "ClinPGx"

# Visão por fármaco: uma coluna por enzima (True/False)
por_farmaco = (
    pd.crosstab(tabela["ddi_nome"], tabela["enzima"])
    .reindex(columns=list(ENZIMAS.values()), fill_value=0)
    .gt(0)
    .reset_index()
)

# Pendentes que importam: só os que têm relação com alguma enzima
revisao = (
    tabela[tabela["origem_correspondencia"].eq("similaridade_alta_pendente")]
    .groupby(["ddi_nome", "nomes_rename"], as_index=False)
    .agg(enzimas=("enzima", lambda s: ", ".join(sorted(s))))
)
revisao["validacao_manual"] = "pendente"


# ==========================================
# 6. SALVAR RESULTADOS
# ==========================================

tabela.to_csv(PASTA_SAIDA / "rename_farmaco_enzima.csv", index=False, encoding="utf-8-sig")
por_farmaco.to_csv(PASTA_SAIDA / "rename_farmacos_por_enzima.csv", index=False, encoding="utf-8-sig")
revisao.to_csv(PASTA_SAIDA / "revisao_manual_prioritaria.csv", index=False, encoding="utf-8-sig")

confirmados = tabela["origem_correspondencia"].eq("match_exato")
auditoria = pd.DataFrame([
    ("Medicamentos RENAME cruzados com ClinPGx", len(rename)),
    ("Sem código do ClinPGx", sem_codigo),
    ("Fármacos com relação a alguma enzima", tabela["ddi_nome"].nunique()),
    ("  com correspondência confirmada", tabela.loc[confirmados, "ddi_nome"].nunique()),
    ("  com correspondência pendente", tabela.loc[~confirmados, "ddi_nome"].nunique()),
    *[(f"Fármacos com relação a {e}", (tabela["enzima"] == e).sum()) for e in ENZIMAS.values()],
], columns=["indicador", "quantidade"])
auditoria.to_csv(PASTA_SAIDA / "auditoria_RENAME_enzimas.csv", index=False, encoding="utf-8-sig")

print()
print(auditoria.to_string(index=False))
print("\nArquivos salvos em:", PASTA_SAIDA)
