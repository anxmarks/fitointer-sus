import re
import unicodedata
from pathlib import Path
import pandas as pd

# ------------------------------------------------------------
# Monta a lista final de compostos das plantas a partir do COCONUT:
# padroniza o nome das plantas, junta duplicatas pelo InChIKey e marca
# os ativos principais (finalizar_ativos.py).
# Compara-se o primeiro bloco do InChIKey (esqueleto da molécula), para
# a mesma substância com estereoquímica anotada de forma diferente não
# virar duas.
# Entrada : resultados/1_plantas/plantas_coconut.csv        (cruzar_plantas.py)
#           resultados/1_plantas/ativos_principais_final.csv (finalizar_ativos.py)
#           dados/manuais/fitodata - Página1.csv
# Saída   : resultados/1_plantas/plantas.csv
#           resultados/1_plantas/compostos.csv              (usado por coleta_chembl.py)
#           resultados/1_plantas/resumo_compostos_por_planta.csv
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]

ARQUIVO_COCONUT = BASE_DIR / "resultados" / "1_plantas" / "plantas_coconut.csv"
ARQUIVO_ATIVOS = BASE_DIR / "resultados" / "1_plantas" / "ativos_principais_final.csv"
ARQUIVO_FITODATA = BASE_DIR / "dados" / "manuais" / "fitodata - Página1.csv"

PASTA_SAIDA = BASE_DIR / "resultados" / "1_plantas"


# ==========================================
# 2. NOMES PADRONIZADOS DAS PLANTAS
# ==========================================

# Chave sem acento -> nome usado no projeto
PLANTAS = {
    "alho": "Alho",
    "babosa": "Babosa",
    "boldo": "Boldo-brasileiro",
    "camomila": "Camomila",
    "curcuma": "Cúrcuma",
    "garra do diabo": "Garra-do-diabo",
    "gengibre": "Gengibre",
    "guarana": "Guaraná",
    "hortela": "Hortelã-pimenta",
    "maracuja": "Maracujá",
    "pata de vaca": "Pata-de-vaca",
    "unha de gato": "Unha-de-gato",
}


def chave_planta(nome):
    nome = unicodedata.normalize("NFKD", str(nome).lower().strip())
    nome = "".join(c for c in nome if not unicodedata.combining(c))
    return re.sub(r"[\s-]+", " ", nome).strip()


def planta_padrao(nome):
    chave = chave_planta(nome)
    if chave not in PLANTAS:
        raise SystemExit(f"Planta sem nome padronizado: {nome!r}. Acrescente em PLANTAS.")
    return PLANTAS[chave]


def juntar_valores(s):
    valores = sorted(set(s.dropna()))
    return "|".join(valores) if valores else None


# ==========================================
# 3. LEITURA DOS ARQUIVOS
# ==========================================

coconut = pd.read_csv(ARQUIVO_COCONUT, dtype=str, low_memory=False)
ativos = pd.read_csv(ARQUIVO_ATIVOS, dtype=str)
fitodata = pd.read_csv(ARQUIVO_FITODATA, dtype=str)
fitodata.columns = fitodata.columns.str.strip()

if "standard_inchi_key" not in coconut.columns:
    raise SystemExit(
        "resultados/1_plantas/plantas_coconut.csv não tem a coluna standard_inchi_key.\n"
        "Rode de novo cruzar_plantas.py (que agora guarda o InChIKey), "
        "filtrar_ativos.py e finalizar_ativos.py."
    )

print("COCONUT:", len(coconut), "linhas")
print("Ativos principais:", len(ativos), "linhas")


# ==========================================
# 4. COMPOSTOS
# ==========================================

# Nome: o do COCONUT; se vazio, o já escolhido em finalizar_ativos.py (ex.: curcumina
# pelo CAS); senão o nome IUPAC
nome_ativo = coconut["identifier"].map(ativos.drop_duplicates("identifier").set_index("identifier")["composto"])

coconut = pd.DataFrame({
    "planta": coconut["Nome popular"].apply(planta_padrao),
    "inchikey": coconut["standard_inchi_key"],
    "nome": coconut["name"].fillna(nome_ativo).fillna(coconut["iupac_name"]),
    "coconut_id": coconut["identifier"],
    "smiles": coconut.get("canonical_smiles"),
    "formula": coconut["molecular_formula"],
    "classe_quimica": coconut["chemical_class"],
})

sem_chave = coconut["inchikey"].isna()
coconut = coconut[~sem_chave].copy()
coconut["inchikey_base"] = coconut["inchikey"].str[:14]

# Uma linha por planta x esqueleto
compostos = coconut.groupby(["planta", "inchikey_base"], as_index=False).agg(
    inchikey=("inchikey", "first"),
    nome=("nome", "first"),
    coconut_id=("coconut_id", juntar_valores),
    smiles=("smiles", "first"),
    formula=("formula", "first"),
    classe_quimica=("classe_quimica", "first"),
)

# Ativos principais: identificados pelo código do COCONUT
ids_ativos = set(ativos["identifier"])
compostos["ativo_principal"] = compostos["coconut_id"].fillna("").apply(
    lambda ids: any(i in ids_ativos for i in ids.split("|"))
)

colunas = [
    "planta", "inchikey", "inchikey_base", "nome", "ativo_principal",
    "coconut_id", "smiles", "formula", "classe_quimica",
]
compostos = compostos[colunas].sort_values(
    ["planta", "ativo_principal", "nome"], ascending=[True, False, True]
)


# ==========================================
# 5. TABELA DE PLANTAS E RESUMO
# ==========================================

plantas = pd.DataFrame({
    "planta": fitodata["Nome popular"].apply(planta_padrao),
    "nome_cientifico": fitodata["Nome cientifico"].str.strip(),
    "id_wiki": fitodata["id_wiki"],
})

resumo = (
    compostos.groupby("planta")
    .agg(compostos=("inchikey_base", "size"),
         ativos_principais=("ativo_principal", "sum"))
    .reindex(plantas["planta"], fill_value=0)
    .reset_index()
)

# Ativos principais que ficaram sem InChIKey no COCONUT
ativos_achados = set(compostos.loc[compostos["ativo_principal"], "coconut_id"].str.split("|").explode())
ativos_sem_chave = ativos[~ativos["identifier"].isin(ativos_achados)]


# ==========================================
# 6. SALVAR RESULTADOS
# ==========================================

plantas.to_csv(PASTA_SAIDA / "plantas.csv", index=False, encoding="utf-8-sig")
compostos.to_csv(PASTA_SAIDA / "compostos.csv", index=False, encoding="utf-8-sig")
resumo.to_csv(PASTA_SAIDA / "resumo_compostos_por_planta.csv", index=False, encoding="utf-8-sig")

print()
print(resumo.to_string(index=False))
print("\nTotal de compostos (planta x molécula):", len(compostos))
print("Linhas do COCONUT sem InChIKey (descartadas):", int(sem_chave.sum()))
if len(ativos_sem_chave):
    print("\nAtivos principais sem InChIKey (não entraram):")
    print(ativos_sem_chave[["Nome popular", "composto", "identifier"]].to_string(index=False))
