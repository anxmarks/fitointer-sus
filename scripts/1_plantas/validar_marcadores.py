import re
import unicodedata
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Monta a tabela de validação dos marcadores: uma linha por planta x marcador
# listado em dados/manuais/plantas_ativos.csv, dizendo o que o projeto achou
# para cada um (compostos no COCONUT, InChIKey, se tem atividade medida).
# As colunas de conferência (fonte_monografia, status, observacao) são
# preenchidas à mão e ficam preservadas quando o script roda de novo.
# Entrada : dados/manuais/plantas_ativos.csv
#           resultados/1_plantas/ativos_principais_final.csv (finalizar_ativos.py)
#           resultados/1_plantas/compostos.csv               (padronizar_compostos.py)
#           resultados/2_chembl/atividades_padronizadas.csv  (padronizar_atividades.py)
# Saída   : dados/manuais/validacao_marcadores.csv
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]

ARQUIVO_ATIVOS_CSV = BASE_DIR / "dados" / "manuais" / "plantas_ativos.csv"
ARQUIVO_ATIVOS_FINAL = BASE_DIR / "resultados" / "1_plantas" / "ativos_principais_final.csv"
ARQUIVO_COMPOSTOS = BASE_DIR / "resultados" / "1_plantas" / "compostos.csv"
ARQUIVO_ATIVIDADES = BASE_DIR / "resultados" / "2_chembl" / "atividades_padronizadas.csv"

ARQUIVO_SAIDA = BASE_DIR / "dados" / "manuais" / "validacao_marcadores.csv"

# Preenchidas à mão; o script nunca sobrescreve o que já estiver aqui
COLUNAS_MANUAIS = ["fonte_monografia", "status", "observacao"]


def sem_acento(txt):
    txt = unicodedata.normalize("NFKD", str(txt))
    return txt.encode("ascii", "ignore").decode().lower().strip()


def chave_planta(nome):
    nome = re.sub(r"\(.*?\)", "", str(nome))
    return re.sub(r"\s+", " ", sem_acento(nome).replace("-", " ")).strip()


def casar(chave, chaves):
    """Mesma tolerância do filtrar_ativos.py: cada arquivo chama a planta de um
    jeito ("boldo", "boldo da terra", "boldo brasileiro")."""
    if chave in chaves:
        return chave
    for outra in chaves:
        if outra.startswith(chave) or chave.startswith(outra):
            return outra
    # "boldo da terra" x "boldo brasileiro": cai para a primeira palavra
    primeira = chave.split()[0]
    iguais = [o for o in chaves if o.split()[0] == primeira]
    return iguais[0] if len(iguais) == 1 else None


# ==========================================
# 2. LEITURA DOS ARQUIVOS
# ==========================================

pesq = pd.read_csv(ARQUIVO_ATIVOS_CSV, sep=",", encoding="utf-8-sig").iloc[:, :2]
pesq.columns = ["planta_csv", "ativos"]

ativos = pd.read_csv(ARQUIVO_ATIVOS_FINAL, encoding="utf-8-sig")
if "marcador" not in ativos.columns:
    raise SystemExit(
        "ativos_principais_final.csv não tem a coluna marcador.\n"
        "Rode de novo filtrar_ativos.py e finalizar_ativos.py."
    )

compostos = pd.read_csv(ARQUIVO_COMPOSTOS, encoding="utf-8-sig")
atividades = pd.read_csv(ARQUIVO_ATIVIDADES, encoding="utf-8-sig", low_memory=False)

# InChIKeys com atividade medida (ativo ou inativo, não só indeterminado)
com_atividade = set(
    atividades.loc[atividades["resultado"].isin(["ativo", "inativo"]), "inchikey"]
)

ativos["_k"] = ativos["Nome popular"].apply(chave_planta)
compostos["_k"] = compostos["planta"].apply(chave_planta)
# Nome padronizado da planta (padronizar_compostos.py), para a tabela final
nome_padrao = dict(zip(compostos["_k"], compostos["planta"]))


# ==========================================
# 3. UMA LINHA POR PLANTA x MARCADOR
# ==========================================

linhas = []
for _, r in pesq.iterrows():
    k = chave_planta(r["planta_csv"])
    # Mesma quebra do filtrar_ativos.py: vírgula, ponto e vírgula, " e " e parênteses
    texto = re.sub(r"[()]", ",", str(r["ativos"]))
    termos = [sem_acento(p) for p in re.split(r"[,;]| e ", texto) if p.strip()]

    da_planta = ativos[ativos["_k"] == casar(k, set(ativos["_k"]))]
    for termo in dict.fromkeys(termos):
        # O marcador guardado em filtrar_ativos.py é a chave do dicionário de
        # tradução, que é um pedaço do termo do CSV
        casa_marcador = da_planta["marcador"].fillna("").map(
            lambda m: any(c and c in termo for c in m.split("|")))
        achados = da_planta[casa_marcador.astype(bool)]
        chaves = sorted(set(achados["standard_inchi_key"].dropna()))
        linhas.append({
            "planta": nome_padrao.get(casar(k, set(nome_padrao)), r["planta_csv"]),
            "planta_csv": r["planta_csv"],
            "marcador": termo,
            "tipo": "substancia" if len(achados) else "classe ou sem tradução",
            "compostos_encontrados": len(achados),
            "nomes_compostos": "; ".join(sorted(achados["composto"].dropna().astype(str))),
            "inchikeys": "; ".join(chaves),
            "com_atividade_medida": sum(1 for c in chaves if c in com_atividade),
        })

tabela = pd.DataFrame(linhas)
for coluna in COLUNAS_MANUAIS:
    tabela[coluna] = None
tabela.loc[:, "status"] = "a_conferir"


# ==========================================
# 4. PRESERVAR A CONFERÊNCIA JÁ FEITA
# ==========================================

if ARQUIVO_SAIDA.exists():
    antes = pd.read_csv(ARQUIVO_SAIDA, encoding="utf-8-sig", dtype=str)
    guardar = ["planta", "marcador"] + [c for c in COLUNAS_MANUAIS if c in antes.columns]
    tabela = (
        tabela.drop(columns=[c for c in COLUNAS_MANUAIS if c in antes.columns])
        .merge(antes[guardar], on=["planta", "marcador"], how="left")
    )
    tabela["status"] = tabela["status"].fillna("a_conferir")
    print("Linhas já conferidas mantidas:", int((tabela["status"] == "conferido").sum()))

tabela.to_csv(ARQUIVO_SAIDA, index=False, encoding="utf-8-sig")


# ==========================================
# 5. RESUMO
# ==========================================

print("Marcadores listados:", len(tabela))
print(tabela["tipo"].value_counts().to_string())
print("\nMarcadores sem composto no COCONUT (precisam de decisão):")
print(tabela.loc[tabela["compostos_encontrados"].eq(0), ["planta", "marcador"]]
      .to_string(index=False))
print("\nPor planta:")
print(tabela.groupby("planta").agg(
    marcadores=("marcador", "size"),
    com_composto=("compostos_encontrados", lambda s: (s > 0).sum()),
    compostos=("compostos_encontrados", "sum"),
    com_atividade=("com_atividade_medida", "sum"),
).to_string())
print("\nArquivo salvo em:", ARQUIVO_SAIDA)
