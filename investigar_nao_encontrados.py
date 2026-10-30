
import pandas as pd
import unicodedata
import re
from pathlib import Path
from rapidfuzz import process, fuzz

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parent

ARQUIVO_DDI = BASE_DIR / "dados" / "drugs.tsv"
ARQUIVO_RENAME = (
    BASE_DIR / "resultados" / "resultados_medicamentos_RENAME_DDI" /
    "medicamentos_RENAME_nao_encontrados_DDI.csv"
)

PASTA_SAIDA = BASE_DIR / "resultados" / "resultados_investigacao_nao_encontrados_DDI"
PASTA_SAIDA.mkdir(parents=True, exist_ok=True)


# ==========================================
# 2. LEITURA
# ==========================================

ddi = pd.read_csv(
    ARQUIVO_DDI,
    sep="\t",
    encoding="utf-8-sig",
    dtype=str
)

nao_encontrados = pd.read_csv(
    ARQUIVO_RENAME,
    encoding="utf-8-sig",
    dtype=str
)

ddi.columns = ddi.columns.str.strip()


# ==========================================
# 3. NORMALIZAÇÃO
# ==========================================

def normalizar(texto):

    if pd.isna(texto):
        return ""

    texto = str(texto).lower().strip()

    texto = unicodedata.normalize("NFKD", texto)

    texto = "".join(
        c for c in texto
        if not unicodedata.combining(c)
    )

    texto = re.sub(r"[^a-z0-9]+", " ", texto)

    return re.sub(r"\s+", " ", texto).strip()


# ==========================================
# 4. SINÔNIMOS DO PHARMGKB
# ==========================================

def extrair_sinonimos(valor):

    if pd.isna(valor) or not str(valor).strip():
        return []

    try:
        return [
            item.strip().strip('"').strip("'")
            for item in next(
                __import__("csv").reader(
                    [str(valor)],
                    skipinitialspace=True
                )
            )
            if item.strip()
        ]
    except Exception:
        return [str(valor).strip()]


aliases = []

for indice, medicamento in ddi.iterrows():

    nome = medicamento["Name"]

    if pd.notna(nome) and normalizar(nome):

        aliases.append({
            "ddi_indice": indice,
            "ddi_nome": nome,
            "nome_alias": nome,
            "tipo_alias": "Name",
            "nome_normalizado": normalizar(nome)
        })

    for sinonimo in extrair_sinonimos(
        medicamento["Generic Names"]
    ):

        if normalizar(sinonimo):

            aliases.append({
                "ddi_indice": indice,
                "ddi_nome": nome,
                "nome_alias": sinonimo,
                "tipo_alias": "Generic Names",
                "nome_normalizado": normalizar(sinonimo)
            })


base_aliases = pd.DataFrame(aliases)

base_aliases = base_aliases.drop_duplicates(
    subset=["ddi_indice", "nome_normalizado"]
)

# Dicionário de busca:
# chave = nome normalizado
# valor = registros correspondentes

nomes_busca = base_aliases[
    "nome_normalizado"
].drop_duplicates().tolist()

mapa_aliases = {}

for _, alias in base_aliases.iterrows():

    chave = alias["nome_normalizado"]

    mapa_aliases.setdefault(chave, []).append(
        alias.to_dict()
    )


# ==========================================
# 5. BUSCA DE CANDIDATOS
# ==========================================

resultados = []

for indice, medicamento in nao_encontrados.iterrows():

    nome_original = medicamento.get("nome_original", "")

    principio = medicamento.get(
        "principio_ativo_base", ""
    )

    componente = medicamento.get(
        "componente", ""
    )

    # Prioriza o princípio ativo-base
    nome_alvo = normalizar(principio)

    if not nome_alvo:
        nome_alvo = normalizar(componente)

    if not nome_alvo:
        continue

    # Busca os 10 nomes mais semelhantes
    candidatos = process.extract(
        nome_alvo,
        nomes_busca,
        scorer=fuzz.ratio,
        limit=10,
        score_cutoff=50
    )

    for nome_candidato, pontuacao, _ in candidatos:

        if pontuacao >= 90:
            classificacao = "similaridade_alta"

        elif pontuacao >= 75:
            classificacao = "similaridade_moderada"

        elif pontuacao >= 50:
            classificacao = "similaridade_baixa"

        else:
            continue

        for alias in mapa_aliases[nome_candidato]:

            dados_ddi = ddi.loc[
                alias["ddi_indice"]
            ]

            linha = medicamento.to_dict()

            linha.update({
                "nome_busca": nome_alvo,
                "candidato_pharmgkb": alias["ddi_nome"],
                "sinonimo_encontrado": alias["nome_alias"],
                "campo_ddi": alias["tipo_alias"],
                "similaridade": pontuacao,
                "classificacao": classificacao,
                "ddi_indice": alias["ddi_indice"],
                "SMILES": dados_ddi.get("SMILES", ""),
                "InChI": dados_ddi.get("InChI", ""),
                "RxNorm Identifiers": dados_ddi.get(
                    "RxNorm Identifiers", ""
                ),
                "PubChem Compound Identifiers": dados_ddi.get(
                    "PubChem Compound Identifiers", ""
                ),
                "validacao_manual": "pendente"
            })

            resultados.append(linha)


# ==========================================
# 6. RESULTADO
# ==========================================

df = pd.DataFrame(resultados)

if not df.empty:

    df = df.drop_duplicates(
        subset=[
            "nome_original",
            "principio_ativo_base",
            "ddi_indice",
            "sinonimo_encontrado"
        ]
    )

    df = df.sort_values(
        by="similaridade",
        ascending=False
    )

arquivo_saida = (
    PASTA_SAIDA /
    "candidatos_revisao_manual_DDI.csv"
)

df.to_csv(
    arquivo_saida,
    index=False,
    encoding="utf-8-sig"
)


# ==========================================
# 7. RELATÓRIO
# ==========================================

print("\n========== INVESTIGAÇÃO DOS NÃO ENCONTRADOS ==========")

print(
    "Registros analisados:",
    len(nao_encontrados)
)

print(
    "Candidatos identificados:",
    len(df)
)

if not df.empty:

    print("\nClassificação dos candidatos:")

    print(df["classificacao"].value_counts())

    print("\nExemplos de candidatos:")

    print(
        df[
            [
                "principio_ativo_base",
                "candidato_pharmgkb",
                "sinonimo_encontrado",
                "similaridade",
                "classificacao"
            ]
        ].head(30).to_string(index=False)
    )

else:
    print("Nenhum candidato encontrado.")

print("\nArquivo salvo em:")
print(arquivo_saida)
