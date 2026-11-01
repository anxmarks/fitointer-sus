import re
import unicodedata
from difflib import SequenceMatcher
from pathlib import Path
import pandas as pd

# ------------------------------------------------------------
# Filtro de metabolismo: liga os medicamentos da RENAME à tabela da FDA
# de substratos, inibidores e indutores de CYP3A4, CYP2C9, CYP2D6 e P-gp.
#
# Dois caminhos de correspondência:
#   a) pelo nome em inglês já encontrado no ClinPGx (ddi_nome);
#   b) busca direta: nome em português da RENAME x nome da FDA, para
#      recuperar o que o cruzamento com o ClinPGx não encontrou
#      (varfarina, tacrolimo, fenitoína...).
#
# Entrada : dados/manuais/rename_principios_ativos.csv
#           resultados/3_rename/medicamentos_RENAME_consolidados_DDI.csv
#           dados/manuais/fda_interacoes.csv  (Table 1 da página "For Healthcare
#           Professionals | FDA's Examples of Drugs that Interact with
#           CYP Enzymes and Transporter Systems", versão de 29/05/2026)
#           dados/manuais/validacao_FDA_RENAME.csv  (opcional: revisão manual)
# Saída   : resultados/3_rename/resultados_RENAME_FDA/
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS E PARÂMETROS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]

ARQUIVO_RENAME = BASE_DIR / "dados" / "manuais" / "rename_principios_ativos.csv"
ARQUIVO_RENAME_DDI = BASE_DIR / "resultados" / "3_rename" / "medicamentos_RENAME_consolidados_DDI.csv"
ARQUIVO_FDA = BASE_DIR / "dados" / "manuais" / "fda_interacoes.csv"
ARQUIVO_VALIDACAO = BASE_DIR / "dados" / "manuais" / "validacao_FDA_RENAME.csv"

PASTA_SAIDA = BASE_DIR / "resultados" / "3_rename" / "resultados_RENAME_FDA"
PASTA_SAIDA.mkdir(parents=True, exist_ok=True)

ENZIMAS = ["CYP3A4", "CYP2C9", "CYP2D6", "P-gp"]

# Ordem usada para escolher a força mais alta quando há mais de uma
ORDEM_FORCA = {"sensivel": 0, "forte": 0, "moderado": 1, "fraco": 2, "": 3}

# Similaridade mínima para um nome virar candidato na busca direta
LIMITE_CANDIDATO = 0.80

# Sais e complementos que não fazem parte do princípio ativo
SAIS = [
    "cloridrato", "dicloridrato", "bromidrato", "sulfato", "bissulfato",
    "fosfato", "citrato", "fumarato", "hemifumarato", "acetato", "maleato",
    "mesilato", "besilato", "succinato", "tartarato", "nitrato", "estolato",
    "etilsuccinato", "valerato", "propionato", "dipropionato", "decanoato",
    "palmitato", "cloreto", "carbonato", "gliconato", "lactato", "oxalato",
]
COMPLEMENTOS = [
    "sodica", "sodico", "potassica", "potassico", "calcica", "calcico",
    "magnesica", "dissodica", "desoproxila", "monoidratada", "di-hidratada",
]


# ==========================================
# 2. LEITURA DOS ARQUIVOS
# ==========================================

rename = pd.read_csv(ARQUIVO_RENAME, encoding="utf-8-sig", dtype=str)
rename_ddi = pd.read_csv(ARQUIVO_RENAME_DDI, encoding="utf-8-sig", dtype=str)
fda = pd.read_csv(ARQUIVO_FDA, encoding="utf-8-sig", dtype=str).fillna("")

rename = rename[rename["tipo"].eq("medicamento")]

print("Componentes de medicamentos da RENAME:", len(rename))
print("RENAME cruzada com ClinPGx:", len(rename_ddi))
print("Linhas da tabela FDA:", len(fda))


# ==========================================
# 3. NORMALIZAÇÃO DOS NOMES
# ==========================================

def normalizar(texto):
    texto = unicodedata.normalize("NFKD", str(texto).lower().strip())
    texto = "".join(c for c in texto if not unicodedata.combining(c))
    texto = re.sub(r"[^a-z0-9-]+", " ", texto)
    return re.sub(r"\s+", " ", texto).strip()


def principio_ativo(texto):
    """'hemifumarato de quetiapina' -> 'quetiapina'; 'varfarina sodica' -> 'varfarina'."""
    texto = normalizar(texto)
    texto = re.sub(rf"^({'|'.join(SAIS)}) de ", "", texto)
    palavras = [p for p in texto.split() if p not in COMPLEMENTOS]
    return " ".join(palavras)


def chave_fonetica(texto):
    """Aproxima a grafia em português e em inglês: phenytoin ~ fenitoina."""
    texto = normalizar(texto).replace("-", " ")
    for de, para in [("ph", "f"), ("th", "t"), ("w", "v"), ("y", "i"),
                     ("k", "c"), ("qu", "c"), ("rr", "r"), ("ss", "s")]:
        texto = texto.replace(de, para)
    # tira a terminação: -ina/-ine/-in, -ol/-ole, -imo/-imus, -ona/-one...
    texto = re.sub(r"(a|e|o|us|um)\b", "", texto)
    return texto.replace(" ", "")


def similaridade(a, b):
    return SequenceMatcher(None, a, b).ratio()


rename["principio_ativo"] = rename["principio_ativo_base"].apply(principio_ativo)
rename["chave"] = rename["principio_ativo"].apply(chave_fonetica)

nomes_fda = pd.DataFrame({"farmaco": fda["farmaco"].unique()})
nomes_fda["chave"] = nomes_fda["farmaco"].apply(chave_fonetica)


# ==========================================
# 4. CAMINHO A: PELO NOME DO CLINPGX
# ==========================================

via_ddi = rename_ddi.assign(farmaco=rename_ddi["ddi_nome"].str.lower().str.strip())
via_ddi = via_ddi[via_ddi["farmaco"].isin(nomes_fda["farmaco"])]
via_ddi = pd.DataFrame({
    "nome_original": via_ddi["nome_original"],
    "principio_ativo": via_ddi["principio_ativo_base"],
    "farmaco": via_ddi["farmaco"],
    "metodo": "clinpgx_" + via_ddi["origem_correspondencia"],
    "similaridade": 1.0,
    "validacao_manual": via_ddi["validacao_manual"].fillna(""),
})


# ==========================================
# 5. CAMINHO B: BUSCA DIRETA PORTUGUÊS x FDA
# ==========================================

linhas = []
for principio, chave in rename[["principio_ativo", "chave"]].drop_duplicates().itertuples(index=False):
    for farmaco, chave_fda in nomes_fda.itertuples(index=False):
        nota = similaridade(chave, chave_fda)
        if nota >= LIMITE_CANDIDATO:
            linhas.append((principio, farmaco, nota))

candidatos = pd.DataFrame(linhas, columns=["principio_ativo", "farmaco", "similaridade"])

# Fica só o melhor candidato de cada princípio ativo
candidatos = (
    candidatos.sort_values("similaridade", ascending=False)
    .drop_duplicates("principio_ativo")
)
candidatos["metodo"] = "busca_direta"
candidatos["validacao_manual"] = candidatos["similaridade"].map(
    lambda n: "aceito_automatico" if n == 1 else "pendente"
)

# Revisão manual: dados/manuais/validacao_FDA_RENAME.csv com as colunas
# principio_ativo, farmaco, validacao_manual (aceito / rejeitado)
if ARQUIVO_VALIDACAO.exists():
    validacao = pd.read_csv(ARQUIVO_VALIDACAO, encoding="utf-8-sig", dtype=str)
    revisados = dict(zip(zip(validacao["principio_ativo"], validacao["farmaco"]),
                         validacao["validacao_manual"].str.strip().str.lower()))
    chaves = list(zip(candidatos["principio_ativo"], candidatos["farmaco"]))
    candidatos["validacao_manual"] = [revisados.get(c, v) for c, v in
                                      zip(chaves, candidatos["validacao_manual"])]

via_direta = rename[["nome_original", "principio_ativo"]].merge(candidatos, on="principio_ativo")


# ==========================================
# 6. JUNTAR OS DOIS CAMINHOS
# ==========================================

correspondencias = pd.concat([via_ddi, via_direta], ignore_index=True)
correspondencias = correspondencias[correspondencias["validacao_manual"].ne("rejeitado")]

# Mesmo medicamento e mesmo fármaco pelos dois caminhos: vale o mais confiável
correspondencias["_ordem"] = correspondencias["metodo"].map(
    {"clinpgx_match_exato": 0, "busca_direta": 1}
).fillna(2)
correspondencias = (
    correspondencias.sort_values(["_ordem", "similaridade"], ascending=[True, False])
    .drop_duplicates(["nome_original", "farmaco"])
    .drop(columns="_ordem")
)

confirmada = (
    correspondencias["metodo"].eq("clinpgx_match_exato")
    | correspondencias["validacao_manual"].isin(["aceito_automatico", "aceito"])
)
correspondencias["correspondencia"] = confirmada.map({True: "confirmada", False: "pendente"})


# ==========================================
# 7. CRUZAMENTO COM A TABELA DA FDA
# ==========================================

cruzado = correspondencias.merge(fda, on="farmaco")

# Uma linha por fármaco x enzima x papel; os nomes da RENAME (sais, formas) ficam juntos
tabela = (
    cruzado
    .groupby(["farmaco", "enzima", "papel"], as_index=False)
    .agg(
        nomes_rename=("nome_original", lambda s: "; ".join(sorted(set(s)))),
        forca=("forca", "first"),
        observacao=("observacao", "first"),
        correspondencia=("correspondencia", "min"),   # "confirmada" < "pendente"
        metodo=("metodo", lambda s: ", ".join(sorted(set(s)))),
    )
)
tabela["fonte"] = "FDA"
tabela["_ordem"] = tabela["forca"].map(ORDEM_FORCA)
tabela = tabela.sort_values(["papel", "enzima", "_ordem", "farmaco"]).drop(columns="_ordem")

# Substratos: são os fármacos que sofrem a interação quando a planta
# inibe ou induz a enzima. É a lista que segue para a tabela de risco.
substratos = tabela[tabela["papel"].eq("substrato")]

# Visão por fármaco: papel em cada enzima (vazio = sem dado na FDA)
por_farmaco = (
    tabela.assign(papel_forca=tabela["papel"] + (" " + tabela["forca"]).where(tabela["forca"].ne(""), ""))
    .pivot_table(index="farmaco", columns="enzima", values="papel_forca",
                 aggfunc=lambda s: "; ".join(s))
    .reindex(columns=ENZIMAS)
    .fillna("")
    .reset_index()
)

# O que falta revisar à mão (já no formato de dados/manuais/validacao_FDA_RENAME.csv)
revisao = (
    correspondencias[correspondencias["correspondencia"].eq("pendente")]
    .groupby(["principio_ativo", "farmaco"], as_index=False)
    .agg(nomes_rename=("nome_original", lambda s: "; ".join(sorted(set(s)))),
         metodo=("metodo", "first"),
         similaridade=("similaridade", "max"))
    .sort_values("similaridade", ascending=False)
)
revisao["validacao_manual"] = "pendente"


# ==========================================
# 8. SALVAR RESULTADOS
# ==========================================

tabela.to_csv(PASTA_SAIDA / "rename_farmaco_enzima_FDA.csv", index=False, encoding="utf-8-sig")
substratos.to_csv(PASTA_SAIDA / "rename_substratos_FDA.csv", index=False, encoding="utf-8-sig")
por_farmaco.to_csv(PASTA_SAIDA / "rename_farmacos_por_enzima_FDA.csv", index=False, encoding="utf-8-sig")
revisao.to_csv(PASTA_SAIDA / "revisao_manual_FDA.csv", index=False, encoding="utf-8-sig")

so_busca_direta = tabela["metodo"].eq("busca_direta")
auditoria = pd.DataFrame([
    ("Componentes de medicamentos da RENAME", len(rename)),
    ("Fármacos da RENAME presentes na tabela FDA", tabela["farmaco"].nunique()),
    ("  encontrados só pela busca direta", tabela.loc[so_busca_direta, "farmaco"].nunique()),
    ("  substratos de alguma enzima", substratos["farmaco"].nunique()),
    ("Correspondências pendentes de revisão", len(revisao)),
    *[(f"Substratos de {e}", (substratos["enzima"] == e).sum()) for e in ENZIMAS],
    *[(f"Inibidores de {e}", ((tabela["papel"] == "inibidor") & (tabela["enzima"] == e)).sum()) for e in ENZIMAS],
    *[(f"Indutores de {e}", ((tabela["papel"] == "indutor") & (tabela["enzima"] == e)).sum()) for e in ENZIMAS],
], columns=["indicador", "quantidade"])
auditoria.to_csv(PASTA_SAIDA / "auditoria_RENAME_FDA.csv", index=False, encoding="utf-8-sig")

print()
print(auditoria.to_string(index=False))
print("\nArquivos salvos em:", PASTA_SAIDA)
