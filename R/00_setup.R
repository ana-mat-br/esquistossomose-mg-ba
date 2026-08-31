# =============================================================================
# 00_setup.R — Ambiente, pacotes e parâmetros globais do estudo
#
# Projeto: Dinâmica espaço-temporal da esquistossomose no corredor endêmico
#          Norte de Minas Gerais + Bahia, 2018-2025
#
# Execute UMA VEZ para instalar dependências; depois é apenas `source()`ado
# pelos demais scripts.
# =============================================================================

# ---- 1. Instalação condicional de pacotes -----------------------------------

pacotes_cran <- c(
  # manipulação
  "tidyverse", "janitor", "data.table", "lubridate", "readxl", "writexl",
  "glue", "here", "jsonlite",
  # geoespacial
  "sf", "sp", "geobr", "spdep", "spatialreg", "tmap", "ggspatial", "classInt",
  # estatística / modelagem
  "sandwich", "lmtest", "MASS", "AER", "car", "geepack", "broom",
  "performance", "marginaleffects", "prais", "trend",
  # varredura espacial (alternativa 100% R ao SaTScan)
  "smerc", "SpatialEpi",
  # GWR / GWPR
  "GWmodel", "spgwr",
  # filtragem espacial (Moran eigenvector)
  "spmoran",
  # IBGE
  "sidrar",
  # apresentação
  "gtsummary", "gt", "patchwork", "scales", "knitr", "kableExtra"
)

instalar_se_faltar <- function(pkgs) {
  faltando <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(faltando)) {
    message("Instalando: ", paste(faltando, collapse = ", "))
    install.packages(faltando, dependencies = TRUE)
  }
  invisible(NULL)
}

# Descomente na primeira execução:
# instalar_se_faltar(pacotes_cran)

# Pacotes que NÃO estão no CRAN (instale manualmente se for usar):
#
#   # leitura de arquivos .dbc do DATASUS
#   remotes::install_github("danicat/read.dbc")
#
#   # download automatizado de microdados DATASUS (SINAN/SIM/SIH)
#   remotes::install_github("rfsaldanha/microdatasus")
#
#   # interface R -> SaTScan (exige o SaTScan instalado no sistema)
#   install.packages("rsatscan")
#
#   # INLA — modelos bayesianos BYM2 espaço-temporais (script 10)
#   install.packages("INLA",
#     repos = c(getOption("repos"),
#               INLA = "https://inla.r-inla-download.org/R/stable"),
#     dep = TRUE)


# ---- 2. Carregamento -------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(sf)
  library(spdep)
  library(glue)
})

options(
  stringsAsFactors = FALSE,
  scipen           = 6,     # evita notação científica em taxas, mas mantém
                            # p-valores muito pequenos legíveis

  digits           = 4,
  timeout          = 600          # downloads IBGE/DATASUS são lentos
)

sf::sf_use_s2(TRUE)


# ---- 3. Parâmetros globais do estudo ---------------------------------------

PARAMS <- list(

  ## Recorte temporal
  ano_ini = 2018L,
  ano_fim = 2025L,

  ## Recorte espacial
  ufs = c("MG", "BA"),
  # Em MG, restringir ao corredor endêmico do Norte:
  #   "imed"  = Região Geográfica Intermediária de Montes Claros — DIVISÃO
  #             REGIONAL VIGENTE DO IBGE (2017). É a definição atual e a que
  #             o IBGE usa hoje em todas as suas publicações.  <-- PADRÃO
  #   "meso"  = Mesorregião "Norte de Minas" — divisão de 1989, DESCONTINUADA
  #             pelo IBGE em 2017, mas ainda dominante na literatura de
  #             esquistossomose (útil para comparar com estudos anteriores)
  #   "lista" = arquivo dados/brutos/municipios_norte_mg.csv fornecido por você
  #
  # OS RECORTES NÃO COINCIDEM. Rode diagnostico_recortes() (script 02) para ver
  # exatamente quais municípios entram/saem antes de fechar a decisão.
  recorte_mg = "imed",

  ## Referências cartográficas
  # Anos candidatos da malha municipal, do mais recente para o mais antigo.
  # O script 02 testa em ordem e usa o primeiro disponível no geobr instalado —
  # assim o pipeline não quebra quando o IBGE/geobr publica um ano novo.
  anos_malha = c(2024L, 2023L, 2022L, 2021L, 2020L),

  ## Vizinhança espacial
  # "queen"  = contiguidade rainha (padrão em estudos de área)
  # "knn"    = k vizinhos mais próximos (garante conectividade; k abaixo)
  viz_tipo  = "queen",
  viz_k     = 5L,
  viz_style = "W",          # padronização por linha

  ## Simulações
  n_sim  = 9999L,           # permutações Monte Carlo (Moran, LISA)
  semente = 20250726L,

  ## Nível de significância
  alfa = 0.05,
  # Correção para múltiplas comparações no LISA:
  #   "fdr" (Benjamini-Hochberg) | "bonferroni" | "none"
  ajuste_p_lisa = "fdr",

  ## Multiplicador de taxa
  base_taxa = 1e5,          # casos por 100.000 hab.

  ## Varredura espacial
  satscan_ubpop   = 0.50,   # janela máx = 50% da população em risco
  satscan_nsim    = 999L,
  # Caminho do executável SaTScan (macOS). Ajuste se necessário.
  satscan_path    = "/Applications/SaTScan.app/Contents/app",

  ## Caminhos
  dir_brutos       = "dados/brutos",
  dir_processados  = "dados/processados",
  dir_figuras      = "saidas/figuras",
  dir_tabelas      = "saidas/tabelas",
  dir_modelos      = "saidas/modelos",
  dir_satscan      = "satscan"
)

set.seed(PARAMS$semente)

for (d in PARAMS[grepl("^dir_", names(PARAMS))]) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ---- 4. Dicionário de covariáveis ------------------------------------------
# Define, de forma explícita, quais variáveis entram nos modelos, sua natureza
# temporal (fixa vs. painel) e a transformação aplicada. Alterar AQUI propaga
# para os scripts 08-10.

COVARIAVEIS <- tibble::tribble(
  ~var,                  ~rotulo,                                  ~fonte,        ~temporal, ~transf,
  "idhm",                "IDHM (geral)",                           "Atlas/PNUD",  "fixa",    "z",
  "idhm_renda",          "IDHM Renda",                             "Atlas/PNUD",  "fixa",    "z",
  "idhm_educacao",       "IDHM Educação",                          "Atlas/PNUD",  "fixa",    "z",
  "idhm_longevidade",    "IDHM Longevidade",                       "Atlas/PNUD",  "fixa",    "z",
  "ivs",                 "Índice de Vulnerabilidade Social",       "Ipea",        "fixa",    "z",
  "gini",                "Índice de Gini da renda domiciliar",     "Atlas/PNUD",  "fixa",    "z",
  # --- vulnerabilidade social (substitutos municipais do IVS, Atlas DH 2010) ---
  # O IVS municipal não é redistribuído pelo Ipeadata; estas são as medidas de
  # privação que compõem o próprio IVS e estão disponíveis por município.
  "pct_pobres",          "% de pessoas pobres",                    "Atlas/PNUD",  "fixa",    "pp10",
  "pct_extrem_pobres",   "% de pessoas extremamente pobres",       "Atlas/PNUD",  "fixa",    "pp10",
  "pct_vulner_pobreza",  "% de vulneráveis à pobreza",             "Atlas/PNUD",  "fixa",    "pp10",
  "renda_per_capita",    "Renda domiciliar per capita (R$)",       "Atlas/PNUD",  "fixa",    "log",
  "taxa_analfabetismo",  "Taxa de analfabetismo (15 anos ou mais)","Atlas/PNUD",  "fixa",    "pp10",
  "pct_agua_rede",       "% domicílios com água de rede geral",    "Censo 2022",  "fixa",    "pp10",
  "pct_esgoto_adeq",     "% domicílios com esgotamento adequado",  "Censo 2022",  "fixa",    "pp10",
  "pct_lixo_coletado",   "% domicílios com coleta de lixo",        "Censo 2022",  "fixa",    "pp10",
  "pct_rural",           "% população rural",                      "Censo 2022",  "fixa",    "pp10",
  # --- variáveis de exposição hídrica específicas de esquistossomose ---
  # Medem contato e contaminação hídrica diretos, e não cobertura de serviço.
  # Distribuição fortemente assimétrica (mediana < 1%, máximo ~30% e ~89%),
  # por isso log1p em vez de escala linear.
  "pct_agua_superficial", "% domicílios abastecidos por rio/córrego/açude",
                                                                   "Censo 2022",  "fixa",    "log1p",
  "pct_esgoto_em_corpo_dagua", "% domicílios com esgoto lançado em rio/lago/córrego",
                                                                   "Censo 2022",  "fixa",    "log1p",
  "cob_agua_snis",       "Cobertura de abastecimento de água (%)", "SNIS",        "painel",  "pp10",
  "cob_esgoto_snis",     "Cobertura de coleta de esgoto (%)",      "SNIS",        "painel",  "pp10",
  "cob_esf",             "Cobertura da Estratégia Saúde Família",  "e-Gestor AB", "painel",  "pp10",
  "dens_demografica",    "Densidade demográfica (hab/km²)",        "IBGE",        "fixa",    "log"
)

# Transformações:
#   "z"    = padronização (média 0, dp 1)  -> RR por 1 desvio-padrão
#   "pp10" = divisão por 10                -> RR por 10 pontos percentuais
#   "log"  = log natural                   -> elasticidade
#   "id"   = sem transformação


# ---- 5. Tema gráfico -------------------------------------------------------

theme_mapa <- function(base_size = 11) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      legend.position   = "right",
      legend.key.height = grid::unit(1.1, "lines"),
      plot.title        = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle     = ggplot2::element_text(colour = "grey30"),
      plot.caption      = ggplot2::element_text(colour = "grey40", size = base_size - 2,
                                                hjust = 0)
    )
}

theme_grafico <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.caption     = ggplot2::element_text(colour = "grey40", size = base_size - 2,
                                               hjust = 0),
      strip.text       = ggplot2::element_text(face = "bold")
    )
}

# Paleta divergente para LISA / razões de risco
PAL_LISA <- c(
  "Alto-Alto"      = "#B2182B",
  "Baixo-Baixo"    = "#2166AC",
  "Alto-Baixo"     = "#EF8A62",
  "Baixo-Alto"     = "#67A9CF",
  "Não significativo" = "grey88"
)

message("Setup carregado. Período ", PARAMS$ano_ini, "-", PARAMS$ano_fim,
        " | UFs: ", paste(PARAMS$ufs, collapse = ", "))
