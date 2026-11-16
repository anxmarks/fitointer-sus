import json
import re
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Camada genética dos pares planta x fármaco que passam por CYP2C9 e CYP2D6.
# Para cada fármaco x gene, junta do ClinPGx:
#   - o melhor nível de evidência das anotações (1A > 1B > 2A > 2B > 3 > 4)
#   - os alelos de função reduzida ou nula (classificação do CPIC)
#   - as diretrizes clínicas (CPIC, DPWG, CPNDS...)
# Em quem tem alelo de função reduzida, a inibição da enzima pela planta soma
# com a genética (fenoconversão): o paciente passa a se comportar como
# metabolizador mais lento, e a exposição ao fármaco aumenta ainda mais.
# Entrada : resultados/4_interacoes/interacoes_composto_farmaco.csv
#           dados/baixados/summaryAnnotations/summary_annotations.tsv
#           dados/baixados/summaryAnnotations/summary_ann_alleles.tsv
#           dados/baixados/guidelineAnnotations.json/*.json
#           (ClinPGx > Downloads: summaryAnnotations.zip e guidelineAnnotations.json.zip)
# Saída   : resultados/4_interacoes/farmacos_genetica.csv     (fármaco x gene)
#           resultados/4_interacoes/interacoes_genetica.csv   (planta x fármaco x gene)
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]

ARQUIVO_INTERACOES = BASE_DIR / "resultados" / "4_interacoes" / "interacoes_composto_farmaco.csv"
PASTA_SUMMARY = BASE_DIR / "dados" / "baixados" / "summaryAnnotations"
PASTA_DIRETRIZES = BASE_DIR / "dados" / "baixados" / "guidelineAnnotations.json"

PASTA_SAIDA = BASE_DIR / "resultados" / "4_interacoes"

GENES = ["CYP2C9", "CYP2D6"]
ORDEM_NIVEL = ["1A", "1B", "2A", "2B", "3", "4"]
FUNCAO_REDUZIDA = ["Decreased function", "No function"]

for caminho in [PASTA_SUMMARY / "summary_annotations.tsv", PASTA_DIRETRIZES]:
    if not caminho.exists():
        raise SystemExit(
            f"Arquivo não encontrado: {caminho}\n"
            "Baixe summaryAnnotations.zip e guidelineAnnotations.json.zip em "
            "clinpgx.org > Downloads e extraia em dados/baixados/ (veja o README)."
        )


# ==========================================
# 2. ANOTAÇÕES DO CLINPGX
# ==========================================

interacoes = pd.read_csv(ARQUIVO_INTERACOES, encoding="utf-8-sig")
anotacoes = pd.read_csv(PASTA_SUMMARY / "summary_annotations.tsv", sep="\t", dtype=str)
alelos = pd.read_csv(PASTA_SUMMARY / "summary_ann_alleles.tsv", sep="\t", dtype=str)

pares = interacoes[interacoes["enzima"].isin(GENES)]
farmacos = sorted(pares["farmaco"].unique())
print("Fármacos da RENAME com interação via", "/".join(GENES) + ":", ", ".join(farmacos))

# Uma linha por anotação x fármaco
anotacoes["farmaco"] = anotacoes["Drug(s)"].str.lower().str.split(";")
anotacoes = anotacoes.explode("farmaco")
anotacoes = anotacoes[anotacoes["Gene"].isin(GENES) & anotacoes["farmaco"].isin(farmacos)]
anotacoes["ordem"] = anotacoes["Level of Evidence"].map({n: i for i, n in enumerate(ORDEM_NIVEL)})


def ordenar_alelos(lista):
    return sorted(set(lista), key=lambda a: [int(n) for n in re.findall(r"\d+", a)] or [0])


# Alelos de função reduzida ou nula, só das anotações de nível 1 e 2
fortes = anotacoes[anotacoes["ordem"] <= ORDEM_NIVEL.index("2B")]
alelos_red = alelos.merge(fortes[["Summary Annotation ID", "Gene", "farmaco"]], on="Summary Annotation ID")
alelos_red = alelos_red[alelos_red["Allele Function"].isin(FUNCAO_REDUZIDA)]
alelos_red["alelo"] = alelos_red["Gene"] + alelos_red["Genotype/Allele"]
alelos_red = (
    alelos_red.groupby(["farmaco", "Gene"])["alelo"]
    .agg(lambda s: ", ".join(ordenar_alelos(s)))
    .rename("alelos_funcao_reduzida").reset_index()
)

resumo_anot = anotacoes.groupby(["farmaco", "Gene"], as_index=False).agg(
    melhor_nivel=("ordem", lambda s: ORDEM_NIVEL[int(s.min())]),
    anotacoes=("Summary Annotation ID", "nunique"),
    categorias=("Phenotype Category", lambda s: "|".join(sorted(set(s.dropna())))),
    urls_nivel1=("URL", lambda s: " ".join(sorted(set(s[anotacoes.loc[s.index, "ordem"] <= 1])))),
)


# ==========================================
# 3. DIRETRIZES CLÍNICAS
# ==========================================

diretrizes = []
for arquivo in sorted(PASTA_DIRETRIZES.glob("*.json")):
    g = json.loads(arquivo.read_text(encoding="utf-8"))["guideline"]
    genes = {x.get("symbol") for x in g.get("relatedGenes", [])}
    drogas = {x.get("name", "").lower() for x in g.get("relatedChemicals", [])}
    for gene in genes & set(GENES):
        for farmaco in drogas & set(farmacos):
            diretrizes.append({
                "farmaco": farmaco, "Gene": gene, "fonte": g.get("source"),
                "recomendacao": bool(g.get("recommendation")),
                "url": f"https://www.clinpgx.org/guidelineAnnotation/{g.get('id')}",
            })
diretrizes = pd.DataFrame(diretrizes, columns=["farmaco", "Gene", "fonte", "recomendacao", "url"])

resumo_dir = diretrizes.groupby(["farmaco", "Gene"], as_index=False).agg(
    diretrizes=("fonte", lambda s: "; ".join(sorted(set(s)))),
    diretriz_com_recomendacao=("recomendacao", "any"),
    urls_diretrizes=("url", lambda s: " ".join(sorted(set(s)))),
)


# ==========================================
# 4. FÁRMACO x GENE
# ==========================================

genetica = (
    pd.DataFrame([(f, g) for f in farmacos for g in GENES], columns=["farmaco", "Gene"])
    .merge(pares[["farmaco", "enzima"]].drop_duplicates(), left_on=["farmaco", "Gene"],
           right_on=["farmaco", "enzima"])
    .drop(columns="enzima")
    .merge(resumo_anot, on=["farmaco", "Gene"], how="left")
    .merge(alelos_red, on=["farmaco", "Gene"], how="left")
    .merge(resumo_dir, on=["farmaco", "Gene"], how="left")
    .rename(columns={"Gene": "gene"})
)
genetica["anotacoes"] = genetica["anotacoes"].fillna(0).astype(int)
genetica["diretriz_com_recomendacao"] = genetica["diretriz_com_recomendacao"].fillna(False).astype(bool)

# Relevância da genética para o par: diretriz com recomendação por genótipo
# ou anotação de nível 1 => alta; nível 2 => moderada; só 3 ou 4 => baixa
nivel = genetica["melhor_nivel"].fillna("")
genetica["relevancia_genetica"] = "baixa"
genetica.loc[nivel.str.startswith("2"), "relevancia_genetica"] = "moderada"
genetica.loc[nivel.str.startswith("1") | genetica["diretriz_com_recomendacao"], "relevancia_genetica"] = "alta"
genetica.loc[genetica["anotacoes"].eq(0) & genetica["diretrizes"].isna(), "relevancia_genetica"] = "sem dados"


# ==========================================
# 5. PLANTA x FÁRMACO x GENE
# ==========================================

def alerta(r):
    if r["relevancia_genetica"] in ("sem dados", "baixa") or pd.isna(r["alelos_funcao_reduzida"]):
        return None
    return (f"Portadores de alelos de função reduzida ou nula ({r['alelos_funcao_reduzida']}) "
            f"já metabolizam menos {r['farmaco']}; a inibição da {r['gene']} pela planta "
            "pode somar ao efeito genético (fenoconversão) e aumentar a exposição.")


planta_farmaco = pares.groupby(["planta", "farmaco", "enzima"], as_index=False).agg(
    nomes_rename=("nomes_rename", "first"),
    nivel_evidencia_planta=("nivel_evidencia", lambda s: max(s, key={"alta": 3, "media": 2, "baixa": 1}.get)),
    compostos=("composto", lambda s: "; ".join(dict.fromkeys(s))),
)
inter_gen = planta_farmaco.merge(genetica, left_on=["farmaco", "enzima"], right_on=["farmaco", "gene"])
inter_gen["alerta"] = inter_gen.apply(alerta, axis=1)
inter_gen = inter_gen.drop(columns="gene").sort_values(["relevancia_genetica", "planta", "farmaco"])


# ==========================================
# 6. SALVAR RESULTADOS
# ==========================================

genetica.to_csv(PASTA_SAIDA / "farmacos_genetica.csv", index=False, encoding="utf-8-sig")
inter_gen.to_csv(PASTA_SAIDA / "interacoes_genetica.csv", index=False, encoding="utf-8-sig")

print()
print(genetica[["farmaco", "gene", "melhor_nivel", "anotacoes", "diretrizes",
                "relevancia_genetica", "alelos_funcao_reduzida"]].to_string(index=False))
print("\nPares planta x fármaco com camada genética:", len(inter_gen))
print(pd.crosstab(inter_gen["nivel_evidencia_planta"], inter_gen["relevancia_genetica"]).to_string())
print("\nArquivos salvos em:", PASTA_SAIDA)
