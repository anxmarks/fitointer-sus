

from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent

PASTA_RESULTADOS = BASE_DIR / "resultados"

ARQUIVO_ENCONTRADOS = (
    PASTA_RESULTADOS /
    "resultados_medicamentos_RENAME_DDI" /
    "medicamentos_RENAME_encontrados_DDI.csv"
)

ARQUIVO_CANDIDATOS = (
    PASTA_RESULTADOS /
    "resultados_investigacao_nao_encontrados_DDI" /
    "candidatos_revisao_manual_DDI.csv"
)

ARQUIVO_FINAL = (
    PASTA_RESULTADOS /
    "medicamentos_RENAME_consolidados_DDI.csv"
)

ARQUIVO_AUDITORIA = (
    PASTA_RESULTADOS /
    "auditoria_consolidacao_DDI_RENAME.csv"
)

LIMITE_SIMILARIDADE = 90



# ==========================================
# 2. LEITURA DOS ARQUIVOS
# ==========================================

encontrados = pd.read_csv(
    ARQUIVO_ENCONTRADOS,
    encoding="utf-8-sig",
    dtype=str
)

candidatos = pd.read_csv(
    ARQUIVO_CANDIDATOS,
    encoding="utf-8-sig",
    dtype=str
)

print("Arquivos carregados.")
print("Encontrados inicialmente:", len(encontrados))
print("Candidatos para análise:", len(candidatos))


# ==========================================
# 3. IDENTIFICAR CANDIDATOS DE ALTA
# ==========================================

candidatos["similaridade"] = pd.to_numeric(
    candidatos["similaridade"],
    errors="coerce"
)

alta_similaridade = candidatos[
    candidatos["similaridade"] >= LIMITE_SIMILARIDADE
].copy()

print(
    "\nCandidatos com similaridade >= 90%:",
    len(alta_similaridade)
)


# ==========================================
# 4. IDENTIFICAR A ORIGEM DOS MATCHES
# ==========================================

encontrados = encontrados.copy()
alta_similaridade = alta_similaridade.copy()

encontrados["origem_correspondencia"] = (
    encontrados.get(
        "tipo_correspondencia",
        pd.Series(
            ["match_exato"] * len(encontrados),
            index=encontrados.index
        )
    )
)

alta_similaridade["origem_correspondencia"] = (
    "similaridade_alta_pendente"
)


# ==========================================
# 5. PADRONIZAR COLUNAS
# ==========================================

# União das colunas dos dois arquivos
todas_colunas = list(
    dict.fromkeys(
        list(encontrados.columns) +
        list(alta_similaridade.columns)
    )
)

encontrados = encontrados.reindex(
    columns=todas_colunas
)

alta_similaridade = alta_similaridade.reindex(
    columns=todas_colunas
)


# ==========================================
# 6. UNIR OS RESULTADOS
# ==========================================

consolidado = pd.concat(
    [encontrados, alta_similaridade],
    ignore_index=True
)


# ==========================================
# 7. REMOVER DUPLICATAS
# ==========================================

# Usa a combinação do registro RENAME e do registro
# PharmGKB, quando esses identificadores estão disponíveis.

chaves = []

if "rename_indice" in consolidado.columns:
    chaves.append("rename_indice")

if "ddi_indice" in consolidado.columns:
    chaves.append("ddi_indice")

if len(chaves) == 2:

    consolidado = consolidado.drop_duplicates(
        subset=chaves,
        keep="first"
    )

else:

    # Alternativa caso os arquivos não tenham índices.
    colunas_identificacao = [
        coluna for coluna in [
            "nome_original",
            "principio_ativo_base",
            "candidato_pharmgkb",
            "nome_DDI",
            "nome_pharmgkb"
        ]
        if coluna in consolidado.columns
    ]

    if colunas_identificacao:
        consolidado = consolidado.drop_duplicates(
            subset=colunas_identificacao,
            keep="first"
        )


# ==========================================
# 8. ORGANIZAR RESULTADOS
# ==========================================

if "similaridade" in consolidado.columns:

    consolidado["similaridade"] = pd.to_numeric(
        consolidado["similaridade"],
        errors="coerce"
    )

    consolidado = consolidado.sort_values(
        by="similaridade",
        ascending=False,
        na_position="last"
    )


# ==========================================
# 9. GERAR AUDITORIA
# ==========================================

auditoria = pd.DataFrame({
    "indicador": [
        "Registros encontrados no cruzamento inicial",
        "Candidatos com similaridade >= 90%",
        "Total de registros após consolidação",
        "Registros com correspondência exata/equivalência",
        "Registros com similaridade alta pendente"
    ],
    "quantidade": [
        len(encontrados),
        len(alta_similaridade),
        len(consolidado),
        (
            consolidado["origem_correspondencia"]
            .isin(["match_exato", "equivalencia_pt_en"])
            .sum()
        ),
        (
            consolidado["origem_correspondencia"]
            == "similaridade_alta_pendente"
        ).sum()
    ]
})


# ==========================================
# ENRIQUECER E PADRONIZAR DADOS PHARMGKB
# ==========================================

# Ler novamente a base original do PharmGKB
ARQUIVO_DDI_ORIGINAL = BASE_DIR / "dados" / "drugs.tsv"

ddi_original = pd.read_csv(
    ARQUIVO_DDI_ORIGINAL,
    sep="\t",
    encoding="utf-8-sig",
    dtype=str
)

ddi_original.columns = ddi_original.columns.str.strip()

# Relacionar os candidatos ao registro original
# do PharmGKB usando o índice preservado.

colunas_pharmgkb = [
    "Name",
    "Generic Names",
    "Trade Names",
    "SMILES",
    "InChI",
    "RxNorm Identifiers",
    "PubChem Compound Identifiers",
    "External Vocabulary",
    "ATC Identifiers"
]

for coluna in colunas_pharmgkb:

    if coluna not in ddi_original.columns:
        continue

    if coluna not in consolidado.columns:
        consolidado[coluna] = pd.NA

# Recuperar os dados originais para cada candidato
if "ddi_indice" in consolidado.columns:

    consolidado["ddi_indice"] = pd.to_numeric(
        consolidado["ddi_indice"],
        errors="coerce"
    )

    for indice, linha in consolidado.iterrows():

        indice_ddi = linha["ddi_indice"]

        if pd.isna(indice_ddi):
            continue

        indice_ddi = int(indice_ddi)

        if indice_ddi not in ddi_original.index:
            continue

        registro_ddi = ddi_original.loc[indice_ddi]

        for coluna in colunas_pharmgkb:

            valor_atual = consolidado.at[indice, coluna]

            if pd.isna(valor_atual) or str(valor_atual).strip() == "":
                consolidado.at[indice, coluna] = (
                    registro_ddi.get(coluna, pd.NA)
                )

# Padronizar os nomes principais
if "Name" in consolidado.columns:
    consolidado["ddi_nome"] = consolidado["ddi_nome"].fillna(
        consolidado["Name"]
    )

if "sinonimo_encontrado" in consolidado.columns:
    consolidado["ddi_nome_correspondente"] = (
        consolidado["ddi_nome_correspondente"].fillna(
            consolidado["sinonimo_encontrado"]
        )
    )

if "candidato_pharmgkb" in consolidado.columns:
    consolidado["ddi_nome"] = (
        consolidado["ddi_nome"].fillna(
            consolidado["candidato_pharmgkb"]
        )
    )

if "sinonimo_encontrado" in consolidado.columns:
    consolidado["ddi_nome_correspondente"] = (
        consolidado["ddi_nome_correspondente"].fillna(
            consolidado["sinonimo_encontrado"]
        )
    )

# Padronizar o método de correspondência
consolidado["metodo_correspondencia"] = (
    consolidado["metodo_correspondencia"].fillna(
        consolidado["origem_correspondencia"]
    )
)

# Remover colunas auxiliares duplicadas
# sem perder as informações que foram padronizadas.
colunas_remover = [
    "Name",
    "Generic Names",
    "Trade Names",
    "candidato_pharmgkb",
    "sinonimo_encontrado"
]

consolidado = consolidado.drop(
    columns=[
        c for c in colunas_remover
        if c in consolidado.columns
    ]
)


# ==========================================
# 10. SALVAR ARQUIVOS
# ==========================================

consolidado.to_csv(
    ARQUIVO_FINAL,
    index=False,
    encoding="utf-8-sig"
)

auditoria.to_csv(
    ARQUIVO_AUDITORIA,
    index=False,
    encoding="utf-8-sig"
)


# ==========================================
# 11. RELATÓRIO FINAL
# ==========================================

print("\n========== CONSOLIDAÇÃO FINAL ==========")

print(auditoria.to_string(index=False))

print("\nArquivo consolidado:")
print(ARQUIVO_FINAL)

print("\nArquivo de auditoria:")
print(ARQUIVO_AUDITORIA)

print("\nProcesso concluído!")

print("\nCOLUNAS DO ARQUIVO CONSOLIDADO:")
print(consolidado.columns.tolist())

print("\nPRIMEIRAS 5 LINHAS:")
print(consolidado.head().to_string())