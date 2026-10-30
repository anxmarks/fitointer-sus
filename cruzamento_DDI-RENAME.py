
import csv
import pandas as pd
import re
import unicodedata
from pathlib import Path

# ==========================================
# 1. CAMINHOS DOS ARQUIVOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parent

ARQUIVO_DDI = BASE_DIR / "dados" / "drugs.tsv"
ARQUIVO_RENAME = BASE_DIR / "dados" / "rename_principios_ativos.csv"

PASTA_SAIDA = BASE_DIR / "resultados" / "resultados_medicamentos_RENAME_DDI"
PASTA_SAIDA.mkdir(parents=True, exist_ok=True)


# ==========================================
# 2. LEITURA DOS BANCOS
# ==========================================

ddi = pd.read_csv(
    ARQUIVO_DDI,
    sep="\t",
    encoding="utf-8-sig",
    dtype=str
)

rename = pd.read_csv(
    ARQUIVO_RENAME,
    sep=",",
    engine="python",
    encoding="utf-8-sig",
    dtype=str
)

ddi.columns = ddi.columns.str.strip()
rename.columns = rename.columns.str.strip()

print("DDI:", ddi.shape)
print("RENAME:", rename.shape)


# ==========================================
# 3. CONFERÊNCIA DAS COLUNAS
# ==========================================

colunas_ddi = ["Name", "Generic Names"]
colunas_rename = [
    "nome_original",
    "componente",
    "principio_ativo_base",
    "associacao"
]

for coluna in colunas_ddi:
    if coluna not in ddi.columns:
        raise ValueError(f"Coluna ausente no DDI: {coluna}")

for coluna in colunas_rename:
    if coluna not in rename.columns:
        raise ValueError(f"Coluna ausente na RENAME: {coluna}")


# ==========================================
# 4. NORMALIZAÇÃO DOS NOMES
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
    texto = re.sub(r"\s+", " ", texto).strip()

    return texto


# ==========================================
# 5. EQUIVALÊNCIAS PORTUGUÊS-INGLÊS
# ==========================================

equivalencias = {
    "aciclovir": "acyclovir",
    "anastrozol": "anastrozole",
    "abiraterona": "abiraterone",
    "desmopressina": "desmopressin",
    "hidrocortisona": "hydrocortisone",
    "abatacepte": "abatacept",
    "abciximabe": "abciximab",
}

# Normaliza as chaves e os valores da tabela
equivalencias_normalizadas = {
    normalizar(chave): normalizar(valor)
    for chave, valor in equivalencias.items()
}


# ==========================================
# 6. EXTRAÇÃO DOS SINÔNIMOS DO PHARMGKB
# ==========================================

def extrair_sinonimos(valor):

    if pd.isna(valor):
        return []

    valor = str(valor).strip()

    if not valor:
        return []

    try:
        nomes = next(
            csv.reader(
                [valor],
                skipinitialspace=True
            )
        )

        return [
            nome.strip().strip('"').strip("'")
            for nome in nomes
            if nome.strip()
        ]

    except Exception:
        return [valor]


aliases = []

for indice, medicamento in ddi.iterrows():

    # Nome principal
    nome = medicamento["Name"]

    if pd.notna(nome):

        aliases.append({
            "ddi_indice": indice,
            "ddi_nome": nome,
            "nome_ddi": nome,
            "nome_normalizado": normalizar(nome),
            "tipo_correspondencia": "nome"
        })

    # Nomes genéricos e sinônimos
    sinonimos = extrair_sinonimos(
        medicamento["Generic Names"]
    )

    for sinonimo in sinonimos:

        nome_normalizado = normalizar(sinonimo)

        if not nome_normalizado:
            continue

        aliases.append({
            "ddi_indice": indice,
            "ddi_nome": medicamento["Name"],
            "nome_ddi": sinonimo,
            "nome_normalizado": nome_normalizado,
            "tipo_correspondencia": "nome_generico"
        })


base_aliases = pd.DataFrame(aliases)

base_aliases = base_aliases[
    base_aliases["nome_normalizado"] != ""
]

base_aliases = base_aliases.drop_duplicates(
    subset=["ddi_indice", "nome_normalizado"]
)

# Dicionário para acelerar a busca
mapa_aliases = {}

for _, alias in base_aliases.iterrows():

    chave = alias["nome_normalizado"]

    mapa_aliases.setdefault(chave, []).append(alias.to_dict())


print("\nTotal de nomes e sinônimos:", len(base_aliases))
print(
    "Sinônimos únicos:",
    base_aliases["nome_normalizado"].nunique()
)


# ==========================================
# 7. CRUZAMENTO RENAME X PHARMGKB
# ==========================================

campos_rename = [
    "principio_ativo_base",
    "componente"
]

campos_quimicos = [
    "SMILES",
    "InChI",
    "RxNorm Identifiers",
    "PubChem Compound Identifiers",
    "External Vocabulary",
    "ATC Identifiers"
]

correspondencias = []
nao_encontrados = []

indices_encontrados = set()

for indice, medicamento in rename.iterrows():

    encontrou = False

    for coluna in campos_rename:

        valor = medicamento[coluna]

        if pd.isna(valor):
            continue

        nome_normalizado = normalizar(valor)

        if not nome_normalizado:
            continue

        # Primeiro: tenta correspondência exata normalizada
        resultados = mapa_aliases.get(nome_normalizado, [])

        metodo = "exata_normalizada"

        # Segundo: se não encontrou e o campo é
        # principio_ativo_base, tenta equivalência
        if (
            not resultados
            and coluna == "principio_ativo_base"
            and nome_normalizado in equivalencias_normalizadas
        ):

            nome_equivalente = equivalencias_normalizadas[
                nome_normalizado
            ]

            resultados = mapa_aliases.get(nome_equivalente, [])

            metodo = "equivalencia_portugues_ingles"

        for resultado in resultados:

            encontrou = True
            indices_encontrados.add(indice)

            linha = medicamento.to_dict()

            linha.update({
                "ddi_nome": resultado["ddi_nome"],
                "ddi_nome_correspondente": resultado["nome_ddi"],
                "campo_rename_utilizado": coluna,
                "campo_ddi_correspondente":
                    resultado["tipo_correspondencia"],
                "metodo_correspondencia": metodo,
                "nome_rename_comparado": valor
            })

            dados_ddi = ddi.loc[resultado["ddi_indice"]]

            for campo in campos_quimicos:

                if campo in ddi.columns:
                    linha[campo] = dados_ddi[campo]

            correspondencias.append(linha)

    if not encontrou:
        linha_nao_encontrada = medicamento.to_dict()
        linha_nao_encontrada["status"] = "nao_encontrado"
        nao_encontrados.append(linha_nao_encontrada)


# ==========================================
# 8. ORGANIZAÇÃO DOS RESULTADOS
# ==========================================

df_correspondencias = pd.DataFrame(correspondencias)
df_nao_encontrados = pd.DataFrame(nao_encontrados)

if not df_correspondencias.empty:

    df_correspondencias = df_correspondencias.drop_duplicates()

    df_correspondencias = df_correspondencias.sort_values(
        by=[
            "nome_original",
            "metodo_correspondencia",
            "ddi_nome"
        ]
    )

# Salva correspondências
df_correspondencias.to_csv(
    PASTA_SAIDA / "medicamentos_RENAME_encontrados_DDI.csv",
    index=False,
    encoding="utf-8-sig"
)

# Salva não encontrados
df_nao_encontrados.to_csv(
    PASTA_SAIDA / "medicamentos_RENAME_nao_encontrados_DDI.csv",
    index=False,
    encoding="utf-8-sig"
)


# ==========================================
# 9. ANÁLISE DAS ASSOCIAÇÕES PARCIAIS
# ==========================================

# Considera os registros da RENAME que possuem associação
# e verifica quantos componentes foram encontrados.

associacoes_parciais = []

if "associacao" in rename.columns:

    grupos = rename.groupby("nome_original", dropna=False)

    for nome, grupo in grupos:

        associacao = grupo["associacao"].astype(str).str.lower().isin(
            ["true", "1", "sim", "yes"]
        ).any()

        if not associacao:
            continue

        indices_grupo = set(grupo.index)

        encontrados_grupo = indices_grupo.intersection(
            indices_encontrados
        )

        if encontrados_grupo and len(encontrados_grupo) < len(indices_grupo):

            linha = {
                "nome_original": nome,
                "total_componentes": len(indices_grupo),
                "componentes_encontrados": len(encontrados_grupo),
                "componentes_nao_encontrados":
                    len(indices_grupo) - len(encontrados_grupo),
                "status": "associacao_parcial"
            }

            associacoes_parciais.append(linha)

df_associacoes_parciais = pd.DataFrame(associacoes_parciais)

df_associacoes_parciais.to_csv(
    PASTA_SAIDA / "associacoes_parcialmente_encontradas.csv",
    index=False,
    encoding="utf-8-sig"
)


# ==========================================
# 10. AUDITORIA DO CRUZAMENTO
# ==========================================

total_rename = len(rename)
total_encontrados = len(indices_encontrados)
total_nao_encontrados = total_rename - total_encontrados

percentual_encontrado = (
    total_encontrados / total_rename * 100
    if total_rename else 0
)

percentual_nao_encontrado = (
    total_nao_encontrados / total_rename * 100
    if total_rename else 0
)

if not df_correspondencias.empty:

    exatas = (
        df_correspondencias["metodo_correspondencia"]
        == "exata_normalizada"
    ).sum()

    equivalencias = (
        df_correspondencias["metodo_correspondencia"]
        == "equivalencia_portugues_ingles"
    ).sum()

    ativos_unicos = (
        df_correspondencias["principio_ativo_base"]
        .dropna()
        .map(normalizar)
        .nunique()
    )

else:
    exatas = 0
    equivalencias = 0
    ativos_unicos = 0


auditoria = pd.DataFrame([
    {"indicador": "Registros originais RENAME", "quantidade": total_rename},
    {"indicador": "Registros únicos encontrados", "quantidade": total_encontrados},
    {"indicador": "Registros não encontrados", "quantidade": total_nao_encontrados},
    {"indicador": "Correspondências exatas", "quantidade": exatas},
    {"indicador": "Correspondências por equivalência", "quantidade": equivalencias},
    {"indicador": "Princípios ativos únicos encontrados", "quantidade": ativos_unicos},
    {"indicador": "Associações parcialmente encontradas", "quantidade": len(df_associacoes_parciais)},
    {"indicador": "Cobertura percentual", "quantidade": round(percentual_encontrado, 2)},
    {"indicador": "Não encontrados percentual", "quantidade": round(percentual_nao_encontrado, 2)}
])

auditoria.to_csv(
    PASTA_SAIDA / "auditoria_cruzamento_DDI_RENAME.csv",
    index=False,
    encoding="utf-8-sig"
)


# ==========================================
# 11. RELATÓRIO NO TERMINAL
# ==========================================

print("\n========== RESULTADO DO CRUZAMENTO ==========")

print("Registros RENAME:", total_rename)
print("Registros únicos encontrados:", total_encontrados)
print("Registros não encontrados:", total_nao_encontrados)

print(f"Cobertura: {percentual_encontrado:.2f}%")
print(f"Não encontrados: {percentual_nao_encontrado:.2f}%")

print("\nCorrespondências exatas:", exatas)
print("Correspondências por equivalência:", equivalencias)
print("Princípios ativos únicos encontrados:", ativos_unicos)

print(
    "Associações parcialmente encontradas:",
    len(df_associacoes_parciais)
)

print("\n========== CORRESPONDÊNCIAS POR MÉTODO ==========")

if not df_correspondencias.empty:
    print(
        df_correspondencias[
            "metodo_correspondencia"
        ].value_counts()
    )

print("\n========== EXEMPLOS ENCONTRADOS ==========")

if not df_correspondencias.empty:

    print(
        df_correspondencias[
            [
                "nome_original",
                "principio_ativo_base",
                "ddi_nome",
                "metodo_correspondencia"
            ]
        ].head(20).to_string(index=False)
    )

print("\n========== EXEMPLOS NÃO ENCONTRADOS ==========")

if not df_nao_encontrados.empty:

    print(
        df_nao_encontrados[
            [
                "nome_original",
                "componente",
                "principio_ativo_base"
            ]
        ].head(20).to_string(index=False)
    )

print("\n========== ASSOCIAÇÕES PARCIAIS ==========")

if not df_associacoes_parciais.empty:
    print(df_associacoes_parciais.head(20).to_string(index=False))
else:
    print("Nenhuma associação parcialmente encontrada.")

print("\nArquivos salvos em:", PASTA_SAIDA)
