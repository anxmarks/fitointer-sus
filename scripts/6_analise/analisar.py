import sqlite3
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Tabelas-resumo da análise, a partir do banco (montar_banco.py).
# A análise principal usa só os ativos principais de cada planta: os compostos
# secundários do COCONUT incluem moléculas que não são da planta (ex.: paclitaxel
# listado na cúrcuma). As tabelas trazem as duas contagens para comparação.
# Entrada : resultados/fitointer.db
#           dados/manuais/rename_principios_ativos.csv
# Saída   : resultados/6_analise/ranking_plantas.csv
#           resultados/6_analise/ranking_farmacos.csv
#           resultados/6_analise/resumo_enzimas.csv
#           resultados/6_analise/origem_evidencia.csv
#           resultados/6_analise/origem_evidencia_total.csv
#           resultados/6_analise/lacunas.csv
#           resultados/6_analise/risco_somado_genetica.csv
#           resultados/6_analise/numeros_do_estudo.csv
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]

ARQUIVO_BANCO = BASE_DIR / "resultados" / "fitointer.db"
ARQUIVO_RENAME = BASE_DIR / "dados" / "manuais" / "rename_principios_ativos.csv"

PASTA_SAIDA = BASE_DIR / "resultados" / "6_analise"
PASTA_SAIDA.mkdir(parents=True, exist_ok=True)

# Fármacos de índice terapêutico estreito entre os substratos da RENAME:
# uma variação pequena na exposição já muda o efeito ou a toxicidade
INDICE_ESTREITO = ["warfarin", "phenytoin", "digoxin", "tacrolimus", "sirolimus", "everolimus"]

if not ARQUIVO_BANCO.exists():
    raise SystemExit(
        f"Arquivo não encontrado: {ARQUIVO_BANCO}\n"
        "Rode antes scripts/5_banco/montar_banco.py (veja o README)."
    )

con = sqlite3.connect(ARQUIVO_BANCO)


# ==========================================
# 2. LEITURA DO BANCO
# ==========================================

inter = pd.read_sql("SELECT * FROM vw_interacoes", con)
compostos = pd.read_sql("SELECT * FROM compostos", con)
atividades = pd.read_sql("SELECT a.*, c.planta, c.nome, c.ativo_principal"
                         " FROM atividades a JOIN compostos c USING (id_composto)", con)
genetica = pd.read_sql("SELECT * FROM farmacos_genetica", con)
substratos = pd.read_sql("SELECT * FROM farmacos_enzimas WHERE papel = 'substrato'", con)
rename = pd.read_csv(ARQUIVO_RENAME, encoding="utf-8-sig", dtype=str)

inter["ativo_principal"] = inter["ativo_principal"].astype(bool)
principais = inter[inter["ativo_principal"]]


ORDEM_NIVEL = {"alta": 3, "media": 2, "baixa": 1}


def melhor_nivel(df):
    """Um nível por par planta x fármaco: o mais alto entre os compostos."""
    pares = df.sort_values("nivel_evidencia", key=lambda s: s.map(ORDEM_NIVEL), ascending=False)
    return pares.drop_duplicates(["planta", "farmaco"])


def contar(df, por, nome):
    """Pares planta x fármaco distintos, pelo nível mais alto de cada par."""
    return (
        melhor_nivel(df)
        .pivot_table(index=por, columns="nivel_evidencia", values="composto",
                     aggfunc="count", fill_value=0)
        .reindex(columns=["alta", "media", "baixa"], fill_value=0)
        .add_prefix(f"{nome}_")
    )


# ==========================================
# 3. RANKING DE PLANTAS
# ==========================================

ranking_plantas = (
    principais.groupby("planta").agg(
        farmacos=("farmaco", "nunique"),
        compostos_envolvidos=("composto", "nunique"),
        enzimas=("enzima", lambda s: "|".join(sorted(set(s)))),
        inibicao=("efeito", lambda s: (s == "inibicao").sum()),
        inducao=("efeito", lambda s: (s == "inducao").sum()),
    )
    .join(contar(principais, "planta", "pares"))
    .join(compostos.groupby("planta").agg(compostos_no_total=("id_composto", "size"),
                                          ativos_principais=("ativo_principal", "sum")))
    .join(principais[principais["farmaco"].isin(INDICE_ESTREITO)]
          .groupby("planta")["farmaco"].nunique().rename("farmacos_indice_estreito"))
    .fillna({"farmacos_indice_estreito": 0})
    .sort_values(["pares_alta", "farmacos"], ascending=False)
    .reset_index()
)
ranking_plantas["farmacos_indice_estreito"] = ranking_plantas["farmacos_indice_estreito"].astype(int)


# ==========================================
# 4. RANKING DE FÁRMACOS
# ==========================================

ranking_farmacos = (
    principais.groupby("farmaco").agg(
        plantas=("planta", "nunique"),
        nomes_rename=("nomes_rename", "first"),
        enzimas=("enzima", lambda s: "|".join(sorted(set(s)))),
        sensibilidade=("sensibilidade_farmaco", "first"),
        plantas_nomes=("planta", lambda s: "; ".join(sorted(set(s)))),
    )
    .join(contar(principais, "farmaco", "pares"))
    .reset_index()
)
ranking_farmacos["indice_estreito"] = ranking_farmacos["farmaco"].isin(INDICE_ESTREITO)
ranking_farmacos = ranking_farmacos.merge(
    genetica[["farmaco", "gene", "melhor_nivel", "relevancia_genetica"]], on="farmaco", how="left"
).sort_values(["pares_alta", "plantas"], ascending=False)


# ==========================================
# 5. ENZIMAS
# ==========================================

resumo_enzimas = (
    principais.groupby("enzima").agg(
        plantas=("planta", "nunique"),
        farmacos=("farmaco", "nunique"),
        compostos=("composto", "nunique"),
    )
    .join(contar(principais, "enzima", "pares"))
    .join(substratos.groupby("enzima")["farmaco"].nunique().rename("substratos_rename"))
    .sort_values("pares_alta", ascending=False)
    .reset_index()
)


# ==========================================
# 6. ORIGEM DA EVIDÊNCIA
# ==========================================

# A triagem do PubChem (qHTS) cobre muitos compostos mas é menos confiável que
# um ensaio dedicado; a correspondência por esqueleto vem de outro estereoisômero
origem = atividades["assay_description"].fillna("").str.startswith("PUBCHEM_BIOASSAY")
atividades["origem"] = origem.map({True: "triagem_pubchem", False: "ensaio_dedicado"})
atividades.loc[atividades["fonte"] == "literatura", "origem"] = "literatura"

ativos_medidos = atividades[atividades["resultado"] == "ativo"]
resumo_origem = dict(
    registros=("id_atividade", "size"),
    compostos=("nome", "nunique"),
    por_esqueleto=("correspondencia", lambda s: (s == "esqueleto").sum()),
)

origem_evidencia = ativos_medidos.groupby(["efeito", "origem"]).agg(**resumo_origem).reset_index()
# Total por origem à parte: somar a coluna compostos do detalhe contaria o mesmo
# composto uma vez por efeito
origem_total = ativos_medidos.groupby("origem").agg(**resumo_origem).reset_index()


# ==========================================
# 7. LACUNAS
# ==========================================

# Ativos principais sem nenhuma atividade medida nas enzimas do estudo
com_dados = set(atividades.loc[atividades["resultado"].isin(["ativo", "inativo"]), "id_composto"])
sem_dados = compostos[compostos["ativo_principal"].astype(bool) & ~compostos["id_composto"].isin(com_dados)]

lacunas = (
    compostos[compostos["ativo_principal"].astype(bool)].groupby("planta")
    .agg(ativos_principais=("id_composto", "size"))
    .join(sem_dados.groupby("planta").agg(
        sem_dados=("id_composto", "size"),
        nomes_sem_dados=("nome", lambda s: "; ".join(sorted(s.dropna()))),
    ))
    .join(principais.groupby("planta")["farmaco"].nunique().rename("farmacos_afetados"))
    .fillna({"sem_dados": 0, "farmacos_afetados": 0})
    .astype({"sem_dados": int, "farmacos_afetados": int})
    .reset_index()
)
lacunas["sem_dados_pct"] = (100 * lacunas["sem_dados"] / lacunas["ativos_principais"]).round(0)


# ==========================================
# 8. RISCO SOMADO COM A GENÉTICA
# ==========================================

risco_genetica = (
    principais.merge(genetica, left_on=["farmaco", "enzima"], right_on=["farmaco", "gene"])
    .query("nivel_evidencia != 'baixa' and relevancia_genetica == 'alta'")
    .groupby(["planta", "farmaco", "enzima"], as_index=False).agg(
        nivel_evidencia=("nivel_evidencia", "first"),
        compostos=("composto", lambda s: "; ".join(sorted(set(s)))),
        nivel_clinpgx=("melhor_nivel", "first"),
        diretrizes=("diretrizes", "first"),
        alelos_funcao_reduzida=("alelos_funcao_reduzida", "first"),
    )
    .sort_values(["farmaco", "planta"])
)


# ==========================================
# 9. NÚMEROS DO ESTUDO
# ==========================================

numeros = pd.DataFrame([
    ("Plantas estudadas", compostos["planta"].nunique()),
    ("Compostos das plantas (COCONUT)", len(compostos)),
    ("Ativos principais", int(compostos["ativo_principal"].sum())),
    ("Ativos principais com atividade medida", int(compostos["ativo_principal"].astype(bool).sum() - len(sem_dados))),
    ("Registros de atividade padronizados", len(atividades)),
    ("Registros da literatura", int((atividades["fonte"] == "literatura").sum())),
    ("Princípios ativos da RENAME", rename["principio_ativo_base"].nunique()),
    ("Fármacos da RENAME substratos das enzimas (FDA)", substratos["farmaco"].nunique()),
    ("Pares planta x fármaco", len(inter.drop_duplicates(["planta", "farmaco"]))),
    ("Pares com ativo principal", len(principais.drop_duplicates(["planta", "farmaco"]))),
    ("Pares de nível alto (ativo principal)", int((melhor_nivel(principais)["nivel_evidencia"] == "alta").sum())),
    ("Pares com risco somado pela genética", len(risco_genetica)),
], columns=["indicador", "quantidade"])


# ==========================================
# 10. SALVAR RESULTADOS
# ==========================================

for nome, df in [("ranking_plantas", ranking_plantas), ("ranking_farmacos", ranking_farmacos),
                 ("resumo_enzimas", resumo_enzimas), ("origem_evidencia", origem_evidencia),
                 ("origem_evidencia_total", origem_total),
                 ("lacunas", lacunas), ("risco_somado_genetica", risco_genetica),
                 ("numeros_do_estudo", numeros)]:
    df.to_csv(PASTA_SAIDA / f"{nome}.csv", index=False, encoding="utf-8-sig")

con.close()

print(numeros.to_string(index=False))
print("\nPlantas (só ativos principais):")
print(ranking_plantas[["planta", "farmacos", "pares_alta", "pares_media", "pares_baixa",
                       "farmacos_indice_estreito", "inducao"]].to_string(index=False))
print("\nFármacos mais afetados:")
print(ranking_farmacos.head(12)[["farmaco", "plantas", "pares_alta", "enzimas",
                                 "indice_estreito", "relevancia_genetica"]].to_string(index=False))
print("\nEnzimas:")
print(resumo_enzimas.to_string(index=False))
print("\nOrigem da evidência (registros ativos):")
print(origem_evidencia.to_string(index=False))
print("\nLacunas:")
print(lacunas[["planta", "ativos_principais", "sem_dados", "sem_dados_pct", "farmacos_afetados"]].to_string(index=False))
print("\nArquivos salvos em:", PASTA_SAIDA)
