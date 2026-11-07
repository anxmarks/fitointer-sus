import sqlite3
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Junta os resultados das etapas num banco SQLite:
#   plantas, compostos, atividades          (Checkpoint 2)
#   farmacos_enzimas                        (RENAME x FDA)
#   interacoes                              (planta x fármaco)
# O banco é refeito do zero a cada execução.
# Entrada : resultados/1_plantas/plantas.csv, compostos.csv
#           resultados/2_chembl/atividades_padronizadas.csv
#           resultados/3_rename/resultados_RENAME_FDA/rename_farmaco_enzima_FDA.csv
#           resultados/4_interacoes/interacoes_composto_farmaco.csv
# Saída   : resultados/fitointer.db
# ------------------------------------------------------------

# ==========================================
# 1. CAMINHOS
# ==========================================

BASE_DIR = Path(__file__).resolve().parents[2]
RES = BASE_DIR / "resultados"

ARQUIVO_PLANTAS = RES / "1_plantas" / "plantas.csv"
ARQUIVO_COMPOSTOS = RES / "1_plantas" / "compostos.csv"
ARQUIVO_ATIVIDADES = RES / "2_chembl" / "atividades_padronizadas.csv"
ARQUIVO_FARMACOS = RES / "3_rename" / "resultados_RENAME_FDA" / "rename_farmaco_enzima_FDA.csv"
ARQUIVO_INTERACOES = RES / "4_interacoes" / "interacoes_composto_farmaco.csv"

ARQUIVO_BANCO = RES / "fitointer.db"


# ==========================================
# 2. ESTRUTURA DO BANCO
# ==========================================

ESQUEMA = """
CREATE TABLE plantas (
    planta          TEXT PRIMARY KEY,
    nome_cientifico TEXT,
    id_wiki         TEXT
);

CREATE TABLE compostos (
    id_composto     INTEGER PRIMARY KEY,
    planta          TEXT NOT NULL REFERENCES plantas(planta),
    nome            TEXT,
    inchikey        TEXT NOT NULL,
    inchikey_base   TEXT,
    ativo_principal INTEGER NOT NULL,
    coconut_id      TEXT,
    smiles          TEXT,
    formula         TEXT,
    classe_quimica  TEXT,
    UNIQUE (planta, inchikey)
);

CREATE TABLE atividades (
    id_atividade      INTEGER PRIMARY KEY,
    id_composto       INTEGER NOT NULL REFERENCES compostos(id_composto),
    enzima            TEXT NOT NULL,
    efeito            TEXT,
    resultado         TEXT,
    forca             TEXT,
    pX                REAL,
    standard_type     TEXT,
    standard_relation TEXT,
    standard_value    REAL,
    standard_units    TEXT,
    activity_comment  TEXT,
    assay_type        TEXT,
    assay_description TEXT,
    chembl_id         TEXT,
    correspondencia   TEXT,
    assay_chembl_id   TEXT,
    fonte             TEXT NOT NULL,
    referencia        TEXT
);

CREATE TABLE farmacos_enzimas (
    farmaco      TEXT NOT NULL,
    nomes_rename TEXT,
    enzima       TEXT NOT NULL,
    papel        TEXT NOT NULL,
    forca        TEXT,
    observacao   TEXT,
    fonte        TEXT,
    PRIMARY KEY (farmaco, enzima, papel)
);

CREATE TABLE interacoes (
    id_interacao         INTEGER PRIMARY KEY,
    id_composto          INTEGER NOT NULL REFERENCES compostos(id_composto),
    enzima               TEXT NOT NULL,
    efeito               TEXT NOT NULL,
    forca_efeito         TEXT,
    melhor_pX            REAL,
    registros_ativos     INTEGER,
    registros_inativos   INTEGER,
    correspondencia      TEXT,
    farmaco              TEXT NOT NULL,
    sensibilidade_farmaco TEXT,
    mecanismo            TEXT,
    nivel_evidencia      TEXT,
    fonte_efeito         TEXT,
    referencias          TEXT,
    fonte_farmaco        TEXT
);

CREATE INDEX idx_atividades_composto ON atividades(id_composto, enzima);
CREATE INDEX idx_interacoes_farmaco ON interacoes(farmaco);

-- Atalho para consultas: interações com o nome da planta e do composto
CREATE VIEW vw_interacoes AS
SELECT c.planta, c.nome AS composto, c.ativo_principal, i.enzima, i.efeito,
       i.forca_efeito, i.melhor_pX, i.farmaco, f.nomes_rename, i.sensibilidade_farmaco,
       i.mecanismo, i.nivel_evidencia, i.referencias
FROM interacoes i
JOIN compostos c ON c.id_composto = i.id_composto
LEFT JOIN farmacos_enzimas f
       ON f.farmaco = i.farmaco AND f.enzima = i.enzima AND f.papel = 'substrato';
"""


# ==========================================
# 3. LEITURA E LIGAÇÃO DAS TABELAS
# ==========================================

def ler(arquivo):
    if not arquivo.exists():
        raise SystemExit(f"Arquivo não encontrado: {arquivo}\nRode antes as etapas anteriores (veja o README).")
    return pd.read_csv(arquivo, encoding="utf-8-sig", low_memory=False)


plantas = ler(ARQUIVO_PLANTAS)
compostos = ler(ARQUIVO_COMPOSTOS)
atividades = ler(ARQUIVO_ATIVIDADES)
farmacos = ler(ARQUIVO_FARMACOS)
interacoes = ler(ARQUIVO_INTERACOES)

compostos.insert(0, "id_composto", range(1, len(compostos) + 1))
compostos["ativo_principal"] = compostos["ativo_principal"].astype(bool).astype(int)


def ligar_composto(df, col_nome):
    """Troca planta + composto pelo id_composto (pelo InChIKey; se faltar, pelo nome)."""
    por_chave = df.merge(compostos[["planta", "inchikey", "id_composto"]],
                         on=["planta", "inchikey"], how="left")
    por_nome = df.merge(
        compostos[["planta", "nome", "id_composto"]].drop_duplicates(["planta", "nome"]),
        left_on=["planta", col_nome], right_on=["planta", "nome"], how="left",
    )["id_composto"]
    por_chave["id_composto"] = por_chave["id_composto"].fillna(
        pd.Series(por_nome.values, index=por_chave.index))
    # Compostos sem nome no COCONUT ficaram com o InChIKey no lugar do nome
    sem_id = por_chave["id_composto"].isna()
    if sem_id.any():
        print(f"  {int(sem_id.sum())} linhas sem composto correspondente (descartadas):")
        print(por_chave.loc[sem_id, ["planta", col_nome]].drop_duplicates().to_string(index=False))
    return por_chave[~sem_id].astype({"id_composto": int})


atividades = ligar_composto(atividades, "nome")
interacoes = ligar_composto(interacoes, "composto")


# ==========================================
# 4. GRAVAR O BANCO
# ==========================================

ARQUIVO_BANCO.unlink(missing_ok=True)
con = sqlite3.connect(ARQUIVO_BANCO)
con.execute("PRAGMA foreign_keys = ON")
con.executescript(ESQUEMA)


def gravar(tabela, df):
    colunas = [c[1] for c in con.execute(f"PRAGMA table_info({tabela})")]
    df = df[[c for c in colunas if c in df.columns]]
    df = df.astype(object).where(df.notna(), None)
    marcadores = ",".join("?" * len(df.columns))
    con.executemany(
        f"INSERT INTO {tabela} ({','.join(df.columns)}) VALUES ({marcadores})",
        df.itertuples(index=False, name=None),
    )


gravar("plantas", plantas)
gravar("compostos", compostos)
gravar("atividades", atividades)
gravar("farmacos_enzimas", farmacos)
gravar("interacoes", interacoes)
con.commit()

for tabela in ["plantas", "compostos", "atividades", "farmacos_enzimas", "interacoes"]:
    n = con.execute(f"SELECT COUNT(*) FROM {tabela}").fetchone()[0]
    print(f"{tabela:<18} {n:>6} linhas")

print("\nExemplo: interações de nível alto ou médio com ativos principais")
print(pd.read_sql("""
    SELECT planta, composto, enzima, forca_efeito, nivel_evidencia,
           COUNT(*) AS farmacos, GROUP_CONCAT(farmaco, ', ') AS exemplos
    FROM vw_interacoes
    WHERE ativo_principal = 1 AND nivel_evidencia IN ('alta', 'media')
    GROUP BY planta, composto, enzima, forca_efeito, nivel_evidencia
    ORDER BY planta
""", con).to_string(index=False))

con.close()
print("\nBanco salvo em:", ARQUIVO_BANCO)
