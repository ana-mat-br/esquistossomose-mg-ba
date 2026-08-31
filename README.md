# Esquistossomose na Região Endêmica, 2018–2025
### Norte de Minas Gerais e Bahia — análise espaço-temporal e determinantes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22216549.svg)](https://doi.org/10.5281/zenodo.22216549)
[![Dados: CC0 1.0](https://img.shields.io/badge/Dados-CC0%201.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
[![Código: MIT](https://img.shields.io/badge/C%C3%B3digo-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Pipeline em R para avaliar a dinâmica espaço-temporal da esquistossomose no
corredor epidemiológico MG–BA e sua associação com indicadores socioeconômicos
e de saneamento.

---

## Dados e licenças

Este repositório contém o código de análise, os bancos agregados por município
e ano e as saídas geradas — tabelas e figuras.

**Microdados individuais do SINAN não são redistribuídos.** As 24.024
notificações que originam as análises contêm ano de nascimento, sexo, raça/cor,
escolaridade, município de residência e os campos de nome da propriedade e da
coleção hídrica de infecção — combinação que permite reidentificação em
municípios de pequeno porte. O que está aqui é o agregado municipal, do qual
não se recupera o indivíduo. Os microdados são de acesso público no DATASUS:
`R/02_obtencao_dados.R` faz o download e `R/03_preparacao_painel.R` reproduz a
agregação, de modo que os arquivos de `dados/processados/` podem ser
regenerados a partir da fonte primária.

Malha municipal (IBGE/`geobr`), população (IBGE/SIDRA) e indicadores do Atlas
do Desenvolvimento Humano são baixados pelo pipeline e devem ser citados em
suas fontes originais.

Licenças: **CC0 1.0** para os dados, **MIT** para o código.

---

## Como rodar

```bash
cd ~/Esquistossomose

# 1. instalar dependências (uma vez): descomente a chamada em R/00_setup.R
R -e 'source("R/00_setup.R"); instalar_se_faltar(pacotes_cran)'

# 2. baixar as bases automáticas (malha + população)
Rscript R/02_obtencao_dados.R

# 3. colocar os arquivos manuais em dados/brutos/
#    -> instruções completas em dados/brutos/README_dados.md

# 4. pipeline completo
Rscript _run_all.R

# ou etapas isoladas
Rscript _run_all.R 5 6 7
```

---

## Estrutura

```
R/
  00_setup.R              pacotes, PARAMS, dicionário de covariáveis, temas
  01_funcoes_auxiliares.R vizinhança, EB, LISA, variância robusta, diagnósticos
  02_obtencao_dados.R     geobr, SIDRA, leitores das bases manuais
  03_preparacao_painel.R  limpeza do SINAN, harmonização, painel município-ano
  04_descritiva_temporal.R Tabela 1, Prais-Winsten, Mann-Kendall, esforço diagnóstico
  05_bayes_empirico.R     EB global e local (Marshall, 1991)
  06_moran_lisa.R         Moran global/EBI, correlograma, LISA, Gi*, bivariado
  07_satscan.R            varredura de Kulldorff (smerc) e SaTScan (rsatscan)
  08_poisson_robusta.R    Poisson modificada (Zou), GEE, sensibilidade
  09_espacial_sar_gwr.R   LM tests, SAR/SEM/SDM/SAC, impactos, ESF, GWR/GWPR
  10_bym2_inla.R          BYM2 espaço-temporal (Knorr-Held I e IV) — opcional
  11_mapas_tabelas.R      cartografia e tabela suplementar
  12_poder_simulacao.R    poder e precisão por Monte Carlo (roda sem os dados)
dados/brutos/             ← seus arquivos (ver README_dados.md)
dados/processados/        painel e agregados municipais
saidas/{figuras,tabelas}/ tabelas e figuras geradas
```

---

## Versões das bases: sempre a mais recente

O pipeline **não fixa anos** que possam sair de circulação. Em três pontos ele
descobre o que existe e usa o mais novo, relatando a escolha:

| Base | Comportamento |
|---|---|
| Malha municipal (geobr) | testa `PARAMS$anos_malha` (2024 → 2020) e usa o primeiro ano disponível **para as duas UFs**, evitando misturar vintages |
| IDHM (Atlas Brasil) | lê a coluna `ano` do CSV e usa o mais recente; force com `ler_idhm(ano_base = 2010)` para comparar com a literatura antiga |
| IVS (Ipea) | mesma regra |
| Recorte do Norte de MG | Região Geográfica Intermediária de Montes Claros (IBGE **2017**, divisão vigente) |

Os anos efetivamente usados ficam registrados em
`saidas/tabelas/00_fluxograma_dados.csv` e em
`dados/processados/metadados_recorte.rds`.

---

## Referências dos métodos

- Marshall RJ. Mapping disease and mortality rates using empirical Bayes
  estimators. *Appl Stat*. 1991;40(2):283–94.
- Assunção RM, Reis EA. A new proposal to adjust Moran's I for population
  density. *Stat Med*. 1999;18(16):2147–62.
- Anselin L. Local indicators of spatial association — LISA. *Geogr Anal*.
  1995;27(2):93–115.
- Kulldorff M. A spatial scan statistic. *Commun Stat Theory Methods*.
  1997;26(6):1481–96.
- Zou G. A modified Poisson regression approach to prospective studies with
  binary data. *Am J Epidemiol*. 2004;159(7):702–6.
- LeSage J, Pace RK. *Introduction to Spatial Econometrics*. CRC Press; 2009.
- Fotheringham AS, Brunsdon C, Charlton M. *Geographically Weighted Regression*.
  Wiley; 2002.
- Knorr-Held L. Bayesian modelling of inseparable space-time variation in
  disease risk. *Stat Med*. 2000;19(17-18):2555–67.
- Riebler A, Sørbye SH, Simpson D, Rue H. An intuitive Bayesian spatial model
  for disease mapping that accounts for scaling. *Stat Methods Med Res*.
  2016;25(4):1145–65.

---

## Reprodutibilidade

Semente fixa em `PARAMS$semente` (20250726). Registre o `sessionInfo()` junto
com os resultados:

```r
writeLines(capture.output(sessionInfo()), "saidas/sessionInfo.txt")
```
