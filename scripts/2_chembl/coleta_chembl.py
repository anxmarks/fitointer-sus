import time
from pathlib import Path
import pandas as pd
import requests

# Caminhos a partir da raiz do projeto, para rodar de qualquer lugar
BASE_DIR = Path(__file__).resolve().parents[2]

# Compostos das plantas (COCONUT), gerados por padronizar_compostos.py
ENTRADA = BASE_DIR / "resultados" / "1_plantas" / "compostos.csv"

#Enzimas e transportador humanos
ALVOS = {
    "CHEMBL340": "CYP3A4",
    "CHEMBL289": "CYP2D6",
    "CHEMBL3397": "CYP2C9",
    "CHEMBL3356": "CYP1A2",
    "CHEMBL3622": "CYP2C19",
    "CHEMBL4302": "P-gp",
}

SAIDA = BASE_DIR / "resultados" / "2_chembl"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
HEADERS = {"User-Agent": "renisus-interacoes/0.1 (projeto academico)"}


def em_lotes(itens, tamanho=50):
    for i in range(0, len(itens), tamanho):
        yield itens[i:i + tamanho]


def chembl_get(recurso, chave, params):
    """GET paginado na API REST do ChEMBL, com novas tentativas se o servidor falhar."""
    url = f"{CHEMBL_API}/{recurso}.json"
    params = {**params, "limit": 1000}
    resultados = []
    while url:
        for tentativa in range(5):
            r = requests.get(url, params=params, headers=HEADERS, timeout=120)
            if r.status_code < 500:
                break
            print(f"  ChEMBL respondeu {r.status_code}, tentando de novo...")
            time.sleep(2 ** tentativa)
        r.raise_for_status()
        dados = r.json()
        resultados.extend(dados[chave])
        proxima = dados["page_meta"]["next"]
        url = f"https://www.ebi.ac.uk{proxima}" if proxima else None
        params = None  # o link da próxima página já traz os filtros
    return resultados


def mapear_chembl(inchikeys):
    linhas = []
    for lote in em_lotes(inchikeys):
        res = chembl_get("molecule", "molecules", {
            "molecule_structures__standard_inchi_key__in": ",".join(lote),
            "only": "molecule_chembl_id,pref_name,molecule_structures",
        })
        for m in res:
            linhas.append({
                "inchikey": m["molecule_structures"]["standard_inchi_key"],
                "chembl_id": m["molecule_chembl_id"],
                "nome_chembl": m["pref_name"],
            })
    return linhas


def atividades_alvos(chembl_ids):
    campos = [
        "molecule_chembl_id", "target_chembl_id", "standard_type",
        "standard_relation", "standard_value", "standard_units",
        "pchembl_value", "assay_type", "assay_description", "document_chembl_id",
    ]
    linhas = []
    for lote in em_lotes(chembl_ids):
        linhas += chembl_get("activity", "activities", {
            "molecule_chembl_id__in": ",".join(lote),
            "target_chembl_id__in": ",".join(ALVOS),
            "only": ",".join(campos),
        })
    df = pd.DataFrame(linhas, columns=campos)
    df["enzima"] = df["target_chembl_id"].map(ALVOS)
    return df


def main():
    SAIDA.mkdir(parents=True, exist_ok=True)

    if not ENTRADA.exists():
        raise SystemExit(
            f"Arquivo não encontrado: {ENTRADA}\n"
            "Ele é gerado por padronizar_compostos.py. Antes, rode cruzar_plantas.py, "
            "filtrar_ativos.py e finalizar_ativos.py com o arquivo do COCONUT em dados/baixados/."
        )

    compostos = pd.read_csv(ENTRADA)
    print(f"{compostos['inchikey'].nunique()} compostos em {ENTRADA}")

    print("\n1/2 Mapeando compostos no ChEMBL...")
    mapa = pd.DataFrame(
        mapear_chembl(compostos["inchikey"].unique().tolist()),
        columns=["inchikey", "chembl_id", "nome_chembl"],
    )
    mapa.to_csv(SAIDA / "compostos_chembl.csv", index=False)
    print(f"{mapa['inchikey'].nunique()} de {compostos['inchikey'].nunique()} "
          "compostos encontrados no ChEMBL")

    print("\n2/2 Atividades contra CYPs e P-gp...")
    ativ = atividades_alvos(mapa["chembl_id"].unique().tolist())
    ativ.to_csv(SAIDA / "atividades_cyp.csv", index=False)
    print(f"{len(ativ)} registros de atividade")
    print(ativ.groupby("enzima").size().to_string())


if __name__ == "__main__":
    main()
