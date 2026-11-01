# fitointer-sus

## Organização

```
scripts/
  1_plantas/     plantas → compostos (COCONUT)
  2_chembl/      compostos → atividade em CYPs e P-gp (ChEMBL)
  3_rename/      fármacos da RENAME → enzimas (ClinPGx e FDA)
  utilitarios/   scripts só para olhar os dados
dados/
  baixados/      arquivos de bases externas, não editar à mão
  manuais/       tabelas montadas ou revisadas pela equipe
resultados/
  1_plantas/  2_chembl/  3_rename/   saídas de cada etapa
```

Os scripts podem ser rodados de qualquer pasta, por exemplo `python scripts/1_plantas/cruzar_plantas.py`.

### Ordem de execução

1. Plantas: `cruzar_plantas.py` → `filtrar_ativos.py` → `finalizar_ativos.py` → `padronizar_compostos.py`
2. ChEMBL: `coleta_chembl.py`
3. RENAME: `cruzamento_DDI-RENAME.py` → `investigar_nao_encontrados.py` → `consolidar_resultados_DDI_RENAME.py` → `cruzamento_RENAME_enzimas.py` e `cruzamento_RENAME_FDA.py`

## Arquivos que precisam ser baixados

Alguns arquivos de dados são grandes demais para o repositório e não estão no git. Baixe-os antes de rodar os scripts e salve-os em `dados/baixados/`.

### `dados/baixados/coconut_csv-11-2026.csv` (COCONUT)

Usado por `cruzar_plantas.py` (cerca de 700 MB).

1. Acesse o site do COCONUT (coconut.naturalproducts.net) e abra a página **Download**.
2. Baixe a base completa em **CSV**. Não use a versão "lite" nem o arquivo de uma coleção só (por exemplo "Wikidata Natural Products"), pois faltam moléculas ou colunas.
3. Salve o arquivo em `dados/baixados/`. Se a versão baixada tiver outra data no nome, ajuste `ARQUIVO_COCONUT` em `cruzar_plantas.py`.

### `dados/baixados/relationships.tsv` (ClinPGx)

Usado por `cruzamento_RENAME_enzimas.py`.

1. Acesse o site do ClinPGx (clinpgx.org) e abra a página **Downloads**.
2. Baixe o arquivo **relationships.zip** (pode ser preciso aceitar a licença de uso).
3. Extraia o `relationships.tsv` e salve-o na pasta `dados/baixados/`.
