# fitointer-sus

## Organização

```
scripts/
  1_plantas/     plantas → compostos (COCONUT)
  2_chembl/      compostos → atividade em CYPs e P-gp (ChEMBL)
  3_rename/      fármacos da RENAME → enzimas (ClinPGx e FDA)
  4_interacoes/  planta × fármaco pela enzima
  5_banco/       banco SQLite com tudo
  6_analise/     tabelas-resumo da análise
  7_graficos/    figuras do relatório e página interativa
  utilitarios/   scripts só para olhar os dados
dados/
  baixados/      arquivos de bases externas, não editar à mão
  manuais/       tabelas montadas ou revisadas pela equipe
resultados/
  1_plantas/  2_chembl/  3_rename/  4_interacoes/  6_analise/  7_graficos/   saídas de cada etapa
  fitointer.db                                     banco SQLite
```

Os scripts podem ser rodados de qualquer pasta, por exemplo `python scripts/1_plantas/cruzar_plantas.py`.

### Ordem de execução

1. Plantas: `cruzar_plantas.py` → `filtrar_ativos.py` → `finalizar_ativos.py` → `padronizar_compostos.py`
2. ChEMBL: `coleta_chembl.py` → `padronizar_atividades.py` (dados da literatura entram por `dados/manuais/atividades_literatura.csv`)
3. RENAME: `cruzamento_DDI-RENAME.py` → `investigar_nao_encontrados.py` → `consolidar_resultados_DDI_RENAME.py` → `cruzamento_RENAME_enzimas.py` e `cruzamento_RENAME_FDA.py`
4. Interações: `cruzar_plantas_farmacos.py` → `camada_genetica.py` (variantes de CYP2C9 e CYP2D6)
5. Banco: `montar_banco.py` (gera `resultados/fitointer.db`; a view `vw_interacoes` junta planta, composto e fármaco, e `vw_interacoes_genetica` acrescenta a genética)
6. Análise: `analisar.py` (lê do banco; a análise principal usa só os ativos principais)
7. Gráficos: o notebook `graficos.ipynb` gera as figuras em PNG e SVG; `gerar_pagina.py` monta a página interativa `resultados/7_graficos/interacoes.html`, que abre no navegador sem servidor

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

### Anotações clínicas do ClinPGx

Usadas por `camada_genetica.py`. Na mesma página **Downloads** do ClinPGx:

1. Baixe **summaryAnnotations.zip** (seção *Summary and Variant Annotations Data*) e **guidelineAnnotations.json.zip** (seção *Clinical Guideline Annotations*). O **clinicalVariants.zip** (seção *Clinical Variant Data*) é opcional.
2. Extraia cada um numa pasta com o mesmo nome dentro de `dados/baixados/`: `summaryAnnotations/`, `guidelineAnnotations.json/` e `clinicalVariants/`.
