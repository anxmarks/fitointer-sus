import json
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Gera a página interativa dos pares planta x fármaco: junta os dados num JSON
# e encaixa no modelo pagina_modelo.html, no lugar de __DADOS__.
# A página é um arquivo só, que abre no navegador sem servidor.
# Entrada : scripts/7_graficos/pagina_modelo.html
#           resultados/4_interacoes/interacoes_planta_farmaco.csv
#           resultados/4_interacoes/interacoes_genetica.csv
#           resultados/3_rename/resultados_RENAME_FDA/rename_substratos_FDA.csv
#           resultados/6_analise/numeros_do_estudo.csv
# Saída   : resultados/7_graficos/interacoes.html
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]

ARQUIVO_MODELO = BASE_DIR / "scripts" / "7_graficos" / "pagina_modelo.html"
ARQUIVO_PARES = BASE_DIR / "resultados" / "4_interacoes" / "interacoes_planta_farmaco.csv"
ARQUIVO_GENETICA = BASE_DIR / "resultados" / "4_interacoes" / "interacoes_genetica.csv"
ARQUIVO_SUBSTRATOS = (
    BASE_DIR / "resultados" / "3_rename" / "resultados_RENAME_FDA" / "rename_substratos_FDA.csv"
)
ARQUIVO_NUMEROS = BASE_DIR / "resultados" / "6_analise" / "numeros_do_estudo.csv"

ARQUIVO_SAIDA = BASE_DIR / "resultados" / "7_graficos" / "interacoes.html"

MARCADOR = "__DADOS__"


def ler(arquivo):
    if not arquivo.exists():
        raise SystemExit(f"Arquivo não encontrado: {arquivo}\nRode antes as etapas anteriores (veja o README).")
    return pd.read_csv(arquivo, encoding="utf-8-sig")


# ==========================================
# 2. DADOS DA PÁGINA
# ==========================================

pares = ler(ARQUIVO_PARES)
genetica = ler(ARQUIVO_GENETICA)
substratos = ler(ARQUIVO_SUBSTRATOS)
numeros = ler(ARQUIVO_NUMEROS)

# Só os pares em que a genética muda o risco entram como alerta
alertas = {
    (r.planta, r.farmaco): {
        "gene": r.enzima, "nivel": r.melhor_nivel, "diretrizes": r.diretrizes,
        "alelos": r.alelos_funcao_reduzida,
        "alerta": r.alerta if isinstance(r.alerta, str) else "",
    }
    for r in genetica.itertuples() if r.relevancia_genetica == "alta"
}
sensibilidade = dict(zip(substratos["farmaco"], substratos["forca"].fillna("")))
observacao = dict(zip(substratos["farmaco"], substratos["observacao"].fillna("")))

# Nomes curtos: a página carrega os dados inteiros, então vale economizar
lista = []
for r in pares.itertuples():
    item = {
        "pl": r.planta, "fa": r.farmaco, "pt": r.nomes_rename, "ni": r.nivel_evidencia,
        "en": r.enzimas, "me": r.mecanismos, "nc": int(r.compostos),
        "ap": r.ativos_principais if isinstance(r.ativos_principais, str) else "",
        "co": r.principais_compostos,
        "se": sensibilidade.get(r.farmaco, ""), "ob": observacao.get(r.farmaco, ""),
    }
    if (r.planta, r.farmaco) in alertas:
        item["ge"] = alertas[(r.planta, r.farmaco)]
    lista.append(item)

dados = {"pares": lista, "numeros": numeros.to_dict("records")}


# ==========================================
# 3. ENCAIXAR NO MODELO
# ==========================================

modelo = ARQUIVO_MODELO.read_text(encoding="utf-8")
if MARCADOR not in modelo:
    raise SystemExit(f"O modelo não tem o marcador {MARCADOR}: {ARQUIVO_MODELO}")

pagina = modelo.replace(MARCADOR, json.dumps(dados, ensure_ascii=False, separators=(",", ":")))
ARQUIVO_SAIDA.parent.mkdir(parents=True, exist_ok=True)
ARQUIVO_SAIDA.write_text(pagina, encoding="utf-8")

print(f"{len(lista)} pares, {sum('ge' in i for i in lista)} com alerta genético")
print(f"Página salva em: {ARQUIVO_SAIDA} ({ARQUIVO_SAIDA.stat().st_size // 1024} KB)")
