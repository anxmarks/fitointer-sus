import re
import numpy as np
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Padroniza as atividades do ChEMBL contra CYPs e P-gp:
#   efeito    : inibicao, ativacao, inducao, substrato, ligacao ou outro
#   resultado : ativo, inativo, inconclusivo ou indeterminado
#   forca     : forte (pX >= 6, até 1 uM), moderada (5 a 6, até 10 uM) ou fraca
# pX é o pchembl_value (-log10 da concentração molar); quando ele falta,
# é calculado do valor em nM.
# Dados da literatura entram por dados/manuais/atividades_literatura.csv,
# com fonte = "literatura".
# Entrada : resultados/1_plantas/compostos.csv
#           resultados/2_chembl/compostos_chembl.csv  (coleta_chembl.py)
#           resultados/2_chembl/atividades_cyp.csv    (coleta_chembl.py)
#           dados/manuais/atividades_literatura.csv
# Saída   : resultados/2_chembl/atividades_padronizadas.csv  (um registro por linha)
#           resultados/2_chembl/compostos_enzimas.csv        (resumo composto x enzima)
#           resultados/2_chembl/resumo_plantas_enzimas.csv
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]

ARQUIVO_COMPOSTOS = BASE_DIR / "resultados" / "1_plantas" / "compostos.csv"
ARQUIVO_MAPA = BASE_DIR / "resultados" / "2_chembl" / "compostos_chembl.csv"
ARQUIVO_ATIVIDADES = BASE_DIR / "resultados" / "2_chembl" / "atividades_cyp.csv"
ARQUIVO_LITERATURA = BASE_DIR / "dados" / "manuais" / "atividades_literatura.csv"

PASTA_SAIDA = BASE_DIR / "resultados" / "2_chembl"

# Até 10 uM conta como ativo (mesmo corte dos painéis do PubChem)
LIMITE_ATIVO_NM = 10_000


# ==========================================
# 2. REGRAS
# ==========================================

TIPOS_CONCENTRACAO = {"IC50", "Ki", "AC50", "Potency", "EC50", "Kd", "Km", "INH"}
TIPOS_SUBSTRATO = {"Km", "Vmax", "Log 1/Km", "Drug metabolism", "Papp", "Ratio_Papp"}
TIPOS_INIBICAO = {"IC50", "Ki", "INH", "Inhibition", "Ratio IC50", "Kinact",
                  "mechanism based inhibition", "AC50"}
TIPOS_LIGACAO = {"Kd", "Delta Tm"}

COMENTARIOS_ATIVO = ("active", "inhibitor [+]", "substrate [+]", "dose-dependent effect",
                     "note: corresponding ic50 reported as active",
                     "cytochrome p450 mechanistic inhibitor")
COMENTARIOS_INATIVO = ("not active", "inhibition < 50%", "is not a substrate [-]")
COMENTARIOS_INCONCLUSIVO = ("inconclusive", "not determined")


def efeito(r):
    tipo = r["standard_type"]
    desc = str(r["assay_description"]).lower()
    coment = str(r["activity_comment"]).lower()

    if tipo in TIPOS_SUBSTRATO or "substrate [" in coment or "drug activation" in desc \
            or "relevant substrates" in desc:
        return "substrato"
    if "activator" in desc:
        return "ativacao"
    # "induction of CYP3A4", "inducers" (e não "docetaxel-induced cytotoxicity")
    if re.search(r"\binducers?\b|induction of (human )?cyp|mrna", desc):
        return "inducao"
    if tipo in TIPOS_INIBICAO or "inhibit" in desc or "kinact" in desc or "inhibitor [" in coment \
            or any(p in desc for p in ("reversal", "accumulation", "modulation", "resistance")):
        return "inibicao"
    if any(p in desc for p in ("transport", "efflux", "uptake", "metaboli")):
        return "substrato"
    if tipo in TIPOS_LIGACAO:
        return "ligacao"
    return "outro"


def valor_nm(r):
    """Valor em nM, só para medidas de concentração."""
    if r["standard_type"] not in TIPOS_CONCENTRACAO or pd.isna(r["standard_value"]):
        return np.nan
    fator = {"nM": 1, "uM": 1_000}.get(r["standard_units"])
    return r["standard_value"] * fator if fator else np.nan


def resultado(r):
    coment = str(r["activity_comment"]).lower()
    if coment.startswith(COMENTARIOS_INATIVO):
        return "inativo"
    if coment.startswith(COMENTARIOS_ATIVO):
        return "ativo"
    if coment.startswith(COMENTARIOS_INCONCLUSIVO):
        return "inconclusivo"
    # Listas curadas ("Clinically relevant substrates/inducers of ... CYP3A4")
    if str(r["assay_description"]).lower().startswith("clinically relevant"):
        return "ativo"

    nm = r["valor_nM"]
    if pd.notna(nm):
        if r["standard_relation"] in (">", ">="):
            return "inativo" if nm >= LIMITE_ATIVO_NM else "indeterminado"
        return "ativo" if nm <= LIMITE_ATIVO_NM else "inativo"

    # % de inibição numa concentração única
    if r["standard_type"] == "Inhibition" and r["standard_units"] == "%" \
            and pd.notna(r["standard_value"]):
        return "ativo" if r["standard_value"] >= 50 else "inativo"

    return "indeterminado"


def forca(px):
    if pd.isna(px):
        return None
    if px >= 6:
        return "forte"
    if px >= 5:
        return "moderada"
    return "fraca"


# ==========================================
# 3. LEITURA DOS ARQUIVOS
# ==========================================

compostos = pd.read_csv(ARQUIVO_COMPOSTOS, encoding="utf-8-sig")
mapa = pd.read_csv(ARQUIVO_MAPA)
ativ = pd.read_csv(ARQUIVO_ATIVIDADES)

if "activity_comment" not in ativ.columns:
    raise SystemExit(
        "atividades_cyp.csv não tem a coluna activity_comment. "
        "Rode de novo scripts/2_chembl/coleta_chembl.py."
    )

print("Registros do ChEMBL:", len(ativ))


# ==========================================
# 4. PADRONIZAÇÃO
# ==========================================

# Dados marcados pelo ChEMBL como suspeitos (ex.: "Outside typical range")
suspeitos = ativ["data_validity_comment"].notna()
ativ = ativ[~suspeitos].copy()

ativ["efeito"] = ativ.apply(efeito, axis=1)
ativ["valor_nM"] = ativ.apply(valor_nm, axis=1)
ativ["resultado"] = ativ.apply(resultado, axis=1)

px_calculado = 9 - np.log10(ativ["valor_nM"].where(ativ["standard_relation"].isin(["=", "<", "<="])))
ativ["pX"] = ativ["pchembl_value"].fillna(px_calculado).round(2)
ativ["forca"] = ativ["pX"].where(ativ["resultado"] == "ativo").apply(forca)

# Liga cada registro aos compostos das plantas
ligacao = mapa.merge(compostos[["planta", "nome", "inchikey", "ativo_principal"]], on="inchikey")
chembl = ligacao.merge(ativ, left_on="chembl_id", right_on="molecule_chembl_id")
# Compostos sem nome no COCONUT ficam com o InChIKey
chembl["nome"] = chembl["nome"].fillna(chembl["inchikey"])
chembl["fonte"] = "ChEMBL"
chembl["referencia"] = chembl["document_chembl_id"]

colunas = [
    "planta", "nome", "inchikey", "ativo_principal", "enzima", "efeito", "resultado",
    "forca", "pX", "standard_type", "standard_relation", "standard_value", "standard_units",
    "activity_comment", "assay_type", "assay_description", "chembl_id", "correspondencia",
    "assay_chembl_id", "fonte", "referencia",
]
tabela = chembl[colunas]


# ==========================================
# 5. DADOS DA LITERATURA
# ==========================================

literatura = pd.DataFrame()
if ARQUIVO_LITERATURA.exists():
    literatura = pd.read_csv(ARQUIVO_LITERATURA, encoding="utf-8-sig")
if len(literatura):
    literatura["fonte"] = "literatura"
    literatura["correspondencia"] = "literatura"
    if "forca" not in literatura or literatura["forca"].isna().all():
        literatura["forca"] = literatura["pX"].apply(forca)
    tabela = pd.concat([tabela, literatura.reindex(columns=colunas)], ignore_index=True)
    print("Registros da literatura:", len(literatura))


# ==========================================
# 6. RESUMO COMPOSTO x ENZIMA
# ==========================================

# Quando há dado do próprio composto (InChIKey exato ou literatura), os de outros
# estereoisômeros (correspondência por esqueleto) não entram no resumo
chave = ["planta", "nome", "enzima", "efeito"]
tem_exato = tabela.groupby(chave)["correspondencia"].transform(lambda s: (s != "esqueleto").any()).astype(bool)
base = tabela[~tem_exato | (tabela["correspondencia"] != "esqueleto")]

resumo = base.groupby(chave, as_index=False, dropna=False).agg(
    inchikey=("inchikey", "first"),
    ativo_principal=("ativo_principal", "first"),
    registros=("resultado", "size"),
    ativos=("resultado", lambda s: (s == "ativo").sum()),
    inativos=("resultado", lambda s: (s == "inativo").sum()),
    melhor_pX=("pX", lambda s: s[base.loc[s.index, "resultado"] == "ativo"].max()),
    correspondencia=("correspondencia", lambda s: "|".join(sorted(set(s)))),
    fontes=("fonte", lambda s: "|".join(sorted(set(s)))),
)

# Evidência: basta um registro ativo (triagem de risco); o número de ativos e
# inativos fica ao lado para a revisão
resumo["evidencia"] = np.select(
    [resumo["ativos"] > 0, resumo["inativos"] > 0], ["ativo", "inativo"], "indeterminado"
)
resumo["forca"] = resumo["melhor_pX"].where(resumo["evidencia"] == "ativo").apply(forca)
resumo = resumo[chave[:2] + ["inchikey", "ativo_principal"] + chave[2:] + [
    "evidencia", "forca", "melhor_pX", "ativos", "inativos", "registros",
    "correspondencia", "fontes",
]].sort_values(["planta", "ativo_principal", "nome", "enzima"], ascending=[True, False, True, True])

# Plantas x enzimas: compostos com inibição ativa
inib = resumo[(resumo["efeito"] == "inibicao") & (resumo["evidencia"] == "ativo")]
por_planta = (
    inib.groupby(["planta", "enzima"])
    .agg(compostos_inibidores=("nome", "nunique"),
         ativos_principais_inibidores=("ativo_principal", "sum"),
         fortes=("forca", lambda s: (s == "forte").sum()))
    .reset_index()
)


# ==========================================
# 7. SALVAR RESULTADOS
# ==========================================

tabela.to_csv(PASTA_SAIDA / "atividades_padronizadas.csv", index=False, encoding="utf-8-sig")
resumo.to_csv(PASTA_SAIDA / "compostos_enzimas.csv", index=False, encoding="utf-8-sig")
por_planta.to_csv(PASTA_SAIDA / "resumo_plantas_enzimas.csv", index=False, encoding="utf-8-sig")

print("Descartados por suspeita do ChEMBL:", int(suspeitos.sum()))
print("\nRegistros por efeito e resultado:")
print(pd.crosstab(tabela["efeito"], tabela["resultado"]).to_string())
print("\nCompostos com inibição ativa, por planta e enzima:")
print(por_planta.pivot(index="planta", columns="enzima", values="compostos_inibidores")
      .fillna(0).astype(int).to_string())
print("\nArquivos salvos em:", PASTA_SAIDA)
