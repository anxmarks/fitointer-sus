import pandas as pd

# 1) Ler os arquivos
plantas = pd.read_csv("dados/fitodata - Página1.csv")
plantas.columns = plantas.columns.str.strip()

colunas = ["identifier", "name", "iupac_name", "molecular_formula",
           "chemical_class", "np_classifier_pathway", "np_classifier_class",
           "organisms", "cas"]
coco = pd.read_csv("dados/coconut_csv-10-2026.csv", usecols=colunas, low_memory=False)

# 2) Espécies de busca (chave = nome popular em minúsculo)
especies_fixas = {
    "boldo": ["plectranthus barbatus", "coleus barbatus"],
    "maracujá": ["passiflora incarnata"],
    "babosa": ["aloe vera", "aloe barbadensis"],
    "camomila": ["matricaria chamomilla", "matricaria recutita"],
    "hortelã": ["mentha piperita"],
}

def genero_especie(nome):
    nome = str(nome).lower().replace("\u00d7", " ").replace("(", " ").replace(")", " ")
    palavras = [p for p in nome.split() if p != "x"]
    return " ".join(palavras[:2])

plantas["busca"] = plantas["Nome cientifico"].apply(genero_especie)

# 3) Limpar a coluna organisms
org = coco["organisms"].fillna("").str.lower()
org = org.str.replace("\u00d7", " ", regex=False)
org = org.str.replace(r"\bx\b", " ", regex=True)
org = org.str.replace(r"\s+", " ", regex=True)
coco["organisms_lower"] = org

# 4) Procurar cada planta dentro de organisms
resultados = []
for _, linha in plantas.iterrows():
    chave_popular = str(linha["Nome popular"]).strip().lower()
    nomes = especies_fixas.get(chave_popular, [linha["busca"]])

    mascara = pd.Series(False, index=coco.index)
    for n in nomes:
        mascara = mascara | coco["organisms_lower"].str.contains(n, regex=False)

    achados = coco[mascara].drop(columns="organisms_lower").copy()
    achados.insert(0, "Nome popular", linha["Nome popular"])
    achados.insert(1, "Nome cientifico", linha["Nome cientifico"])
    achados.insert(2, "id_wiki", linha["id_wiki"])
    resultados.append(achados)

    print(f"{linha['Nome popular']:<15} {len(achados):>6} moléculas")

final = pd.concat(resultados, ignore_index=True)
final.to_csv("resultados/plantas_coconut.csv", index=False)
print("\nTotal de linhas:", len(final))