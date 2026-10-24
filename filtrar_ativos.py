import re
import unicodedata
import pandas as pd

# ------------------------------------------------------------
# Filtra as moléculas do COCONUT (já separadas por planta) para ficar
# só com os ATIVOS PRINCIPAIS listados no seu CSV de pesquisa.
# Entrada : resultados/plantas_coconut.csv  +  dados/plantas_ativos.csv
# Saída   : resultados/ativos_principais.csv
# ------------------------------------------------------------

def sem_acento(txt):
    txt = unicodedata.normalize("NFKD", str(txt))
    return txt.encode("ascii", "ignore").decode().lower().strip()

def chave_planta(nome):
    # "Garra-do-diabo (Harpagophytum...)" -> "garra do diabo"
    nome = re.sub(r"\(.*?\)", "", str(nome))
    return re.sub(r"\s+", " ", sem_acento(nome).replace("-", " ")).strip()

# Dicionário português -> inglês (o COCONUT usa nomes em inglês).
# Sem acentos. Um termo pode virar vários nomes. Acrescente o que faltar.
traducao = {
    "alicina": ["allicin"],
    "aliina": ["alliin"],
    "s-alilcisteina": ["s-allylcysteine", "allylcysteine"],
    "dissulfeto de dialila": ["diallyl disulfide"],
    "aloina": ["aloin", "barbaloin"],
    "aloe-emodina": ["aloe-emodin"],
    "acemana": ["acemannan"],
    "boldina": ["boldine"],
    "ascaridol": ["ascaridole"],
    "forskolina": ["forskolin", "coleonol"],
    "coleonol": ["coleonol"],
    "barbatusina": ["barbatusin"],
    "acido rosmarinico": ["rosmarinic acid"],
    "apigenina": ["apigenin"],
    "luteolina": ["luteolin"],
    "bisabolol": ["bisabolol"],
    "matricina": ["matricin"],
    "camazuleno": ["chamazulene"],
    "curcuminoides": ["curcumin", "demethoxycurcumin", "bisdemethoxycurcumin"],
    "curcumina": ["curcumin"],
    "turmerona": ["turmerone"],
    "harpagosideo": ["harpagoside"],
    "harpagideo": ["harpagide"],
    "gingerois": ["gingerol"],
    "shogaois": ["shogaol"],
    "zingerona": ["zingerone"],
    "cafeina": ["caffeine"],
    "metilxantinas": ["caffeine", "theobromine", "theophylline"],
    "mentol": ["menthol"],
    "mentona": ["menthone"],
    "vitexina": ["vitexin"],
    "isovitexina": ["isovitexin"],
    "orientina": ["orientin"],
    "kaempferitrina": ["kaempferitrin"],
    "mitrafilina": ["mitraphylline"],
    "isomitrafilina": ["isomitraphylline"],
    "pteropodina": ["pteropodine"],
}
chaves_ordenadas = sorted(traducao, key=len, reverse=True)

# Números CAS (identificam a molécula mesmo quando o COCONUT deixa o nome vazio).
# CONFIRA cada um no PubChem antes de usar no relatório.
cas_por_chave = {
    "curcumina": ["458-37-7"],
    "curcuminoides": ["458-37-7", "22608-11-3", "24939-16-0"],  # curcumina, demetoxi-, bisdemetoxicurcumina
    "shogaois": ["555-66-8"],                                    # 6-shogaol
    "s-alilcisteina": ["21593-77-1"],
}

# 1) Ativos de cada planta a partir do seu CSV de pesquisa
pesq = pd.read_csv("dados/plantas_ativos.csv", sep=",", encoding="utf-8-sig")
print("CSV de pesquisa:", pesq.shape[0], "linhas x", pesq.shape[1], "colunas")
if pesq.shape[1] < 2 or pesq.iloc[:, 1].isna().all():
    raise SystemExit(
        "\nERRO: o arquivo dados/plantas_ativos.csv não está no formato certo "
        "(precisa de 12 linhas e 2 colunas: planta e principais ativos).\n"
        "Substitua pelo plantas_ativos.csv limpo e confira que a leitura usa sep=\",\"."
    )
pesq = pesq.iloc[:, :2]
pesq.columns = ["planta", "ativos"]

alvos = {}          # chave da planta -> lista de nomes em inglês
alvos_cas = {}      # chave da planta -> lista de CAS
sem_traducao = {}   # termos do CSV que o dicionário não conhece
for _, r in pesq.iterrows():
    k = chave_planta(r["planta"])
    texto = re.sub(r"[()]", ",", str(r["ativos"]))
    pedacos = [sem_acento(p) for p in re.split(r"[,;]| e ", texto) if p.strip()]
    nomes = []
    cas_lista = []
    for p in pedacos:
        achou = [c for c in chaves_ordenadas if c in p]
        if achou:
            for c in achou:
                nomes += traducao[c]
                cas_lista += cas_por_chave.get(c, [])
        else:
            sem_traducao.setdefault(k, []).append(p)
    alvos[k] = sorted(set(nomes))
    alvos_cas[k] = sorted(set(cas_lista))

# 2) Moléculas do COCONUT por planta -> filtrar pelos ativos
df = pd.read_csv("resultados/plantas_coconut.csv", low_memory=False)
df["_nome"] = df["name"].fillna("").str.lower()
df["_iupac"] = df["iupac_name"].fillna("").str.lower()
df["_syn"] = df["synonyms"].fillna("").str.lower() if "synonyms" in df.columns else ""
df["_k"] = df["Nome popular"].apply(chave_planta)

def achar_chave(k):
    if k in alvos:
        return k
    for outro in alvos:                       # "hortela" ~ "hortela pimenta"
        if outro.startswith(k) or k.startswith(outro):
            return outro
    return None

df["_cas"] = df["cas"].fillna("").astype(str)

partes = []
for k in df["_k"].unique():
    sub = df[df["_k"] == k]
    chave = achar_chave(k)
    termos = alvos.get(chave, [])
    cas_alvo = alvos_cas.get(chave, [])

    m_nome = pd.Series(False, index=sub.index)
    for t in termos:
        m_nome = m_nome | sub["_nome"].str.contains(t, regex=False) \
                        | sub["_iupac"].str.contains(t, regex=False) \
                        | sub["_syn"].str.contains(t, regex=False)
    m_cas = pd.Series(False, index=sub.index)
    for c in cas_alvo:
        m_cas = m_cas | sub["_cas"].str.contains(c, regex=False)

    achados = sub[m_nome | m_cas].copy()
    achados["encontrado_por"] = ["nome e cas" if a and b else ("nome" if a else "cas")
                                 for a, b in zip(m_nome[m_nome | m_cas], m_cas[m_nome | m_cas])]
    achados["ativo_buscado"] = ", ".join(termos)
    partes.append(achados)
    print(f"{sub['Nome popular'].iloc[0]:<15} {len(sub):>5} moléculas -> {len(achados):>3} ativos "
          f"(só pelo CAS: {int((m_cas & ~m_nome).sum())})")

final = pd.concat(partes, ignore_index=True).drop(columns=["_nome", "_iupac", "_syn", "_k", "_cas"])
final.to_csv("resultados/ativos_principais.csv", index=False)
print("\nTotal:", len(final))

if sem_traducao:
    print("\nTermos do seu CSV sem tradução no dicionário (adicione se forem ativos):")
    for k, v in sem_traducao.items():
        print(f"  {k}: {v}")