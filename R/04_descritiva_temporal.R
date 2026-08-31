# =============================================================================
# 04_descritiva_temporal.R — Descritiva e análise de tendência temporal
#
#   - Caracterização da área de estudo e dos desfechos
#   - Tendência por regressão de Prais-Winsten (corrige autocorrelação de
#     1ª ordem nos resíduos da série anual) -> Variação Percentual Anual (VPA)
#   - Mann-Kendall + declive de Sen (não paramétrico, robusto a n pequeno)
#   - Decomposição da mudança: quanto da variação da positividade é atribuível
#     a mudança de esforço diagnóstico (cobertura de exame) vs. transmissão
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

library(prais)
library(trend)

painel    <- readRDS(file.path(PARAMS$dir_processados, "painel.rds"))
agregado  <- readRDS(file.path(PARAMS$dir_processados, "agregado.rds"))
malha_est <- readRDS(file.path(PARAMS$dir_processados, "malha_est.rds"))


# =============================================================================
# 1. TABELA 1 — caracterização da área de estudo
# =============================================================================

# Covariáveis contextuais efetivamente presentes na base (o pipeline roda
# mesmo que IDHM/IVS/SNIS não tenham sido baixados)
vars_contexto <- intersect(
  c("idhm", "idhm_renda", "ivs", "gini", "pct_agua_rede", "pct_esgoto_adeq",
    "pct_agua_superficial", "pct_esgoto_em_corpo_dagua", "pct_lixo_coletado",
    "pct_rural", "dens_demografica"),
  names(agregado)
)
message("Covariáveis contextuais disponíveis: ",
        paste(vars_contexto, collapse = ", "))

tab1 <- agregado |>
  dplyr::group_by(regiao) |>
  dplyr::summarise(
    n_municipios     = dplyr::n(),
    pop_mediana      = median(pop_central, na.rm = TRUE),
    casos_totais     = sum(casos_sinan),
    inc_media_bruta  = sum(casos_sinan) / sum(pessoa_ano) * PARAMS$base_taxa,
    exames_totais    = sum(examinados, na.rm = TRUE),
    positivos_totais = sum(positivos, na.rm = TRUE),
    positividade_pct = dplyr::if_else(
      sum(examinados, na.rm = TRUE) > 0,
      sum(positivos, na.rm = TRUE) / sum(examinados, na.rm = TRUE) * 100,
      NA_real_),
    dplyr::across(dplyr::all_of(vars_contexto),
                  list(mediana = ~ median(.x, na.rm = TRUE),
                       q1 = ~ quantile(.x, .25, na.rm = TRUE),
                       q3 = ~ quantile(.x, .75, na.rm = TRUE)),
                  .names = "{.col}__{.fn}"),
    .groups = "drop"
  )

# Comparação entre regiões (Mann-Whitney para contínuas assimétricas)
comparar_regioes <- function(dados, vars) {
  purrr::map_dfr(vars, function(v) {
    if (!v %in% names(dados)) return(NULL)
    x <- dados[[v]][dados$regiao == "Norte de Minas"]
    y <- dados[[v]][dados$regiao == "Bahia"]
    if (sum(!is.na(x)) < 3 || sum(!is.na(y)) < 3) return(NULL)
    tt <- wilcox.test(x, y, exact = FALSE)
    tibble::tibble(
      variavel = v,
      mediana_NMG = median(x, na.rm = TRUE),
      mediana_BA  = median(y, na.rm = TRUE),
      W = unname(tt$statistic),
      p_valor = tt$p.value
    )
  })
}

vars_comp <- intersect(
  c(vars_contexto, "inc_sinan_bruta", "positividade"), names(agregado))

tab1_teste <- comparar_regioes(agregado, vars_comp) |>
  dplyr::mutate(p_fmt = fmt_p(p_valor))

salvar_tabela(tab1, "01_caracterizacao_area")
salvar_tabela(tab1_teste, "01_comparacao_regioes")
print(as.data.frame(tab1_teste))


# =============================================================================
# 2. SÉRIE TEMPORAL AGREGADA
# =============================================================================

serie <- painel |>
  dplyr::group_by(regiao, ano) |>
  dplyr::summarise(
    casos      = sum(casos_sinan, na.rm = TRUE),
    pop        = sum(populacao, na.rm = TRUE),
    examinados = sum(examinados, na.rm = TRUE),
    positivos  = sum(positivos, na.rm = TRUE),
    n_mun_com_pce = sum(!is.na(positividade)),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    incidencia   = casos / pop * PARAMS$base_taxa,
    positividade = dplyr::if_else(examinados > 0, positivos / examinados * 100,
                                  NA_real_),
    cob_exame    = examinados / pop * 100
  )

serie_total <- painel |>
  dplyr::group_by(ano) |>
  dplyr::summarise(
    casos = sum(casos_sinan, na.rm = TRUE),
    pop = sum(populacao, na.rm = TRUE),
    examinados = sum(examinados, na.rm = TRUE),
    positivos = sum(positivos, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(regiao = "Total",
                incidencia = casos / pop * PARAMS$base_taxa,
                positividade = dplyr::if_else(examinados > 0,
                                              positivos / examinados * 100, NA_real_),
                cob_exame = examinados / pop * 100)

serie_full <- dplyr::bind_rows(serie, serie_total)
salvar_tabela(serie_full, "02_serie_temporal")


# =============================================================================
# 3. TENDÊNCIA — REGRESSÃO DE PRAIS-WINSTEN
# =============================================================================
# Modelo: log10(taxa_t) = b0 + b1 * ano_t + e_t,  e_t ~ AR(1)
# VPA (%) = (10^b1 - 1) * 100
#
# Prais-Winsten é o padrão em séries epidemiológicas anuais brasileiras porque
# corrige a autocorrelação serial que infla a significância do OLS.
# Com apenas 8 pontos (2018-2025) o poder é baixo: reporte o IC, não só o p.

vpa_prais <- function(dados, var_taxa, rotulo) {

  d <- dados |>
    dplyr::filter(is.finite(.data[[var_taxa]]), .data[[var_taxa]] > 0) |>
    dplyr::mutate(log_taxa = log10(.data[[var_taxa]]),
                  ano_c = ano - min(ano)) |>
    dplyr::arrange(ano)

  if (nrow(d) < 4) {
    return(tibble::tibble(serie = rotulo, n_pontos = nrow(d),
                          VPA = NA_real_, VPA_li = NA_real_, VPA_ls = NA_real_,
                          p_valor = NA_real_, rho_AR1 = NA_real_,
                          classificacao = "insuficiente"))
  }

  pw <- prais::prais_winsten(log_taxa ~ ano_c, data = d, index = "ano_c")
  sm <- summary(pw)
  b1 <- coef(pw)["ano_c"]
  se <- sm$coefficients["ano_c", "Std. Error"]
  p  <- sm$coefficients["ano_c", ncol(sm$coefficients)]

  vpa    <- (10^b1 - 1) * 100
  vpa_li <- (10^(b1 - 1.96 * se) - 1) * 100
  vpa_ls <- (10^(b1 + 1.96 * se) - 1) * 100

  tibble::tibble(
    serie    = rotulo,
    n_pontos = nrow(d),
    b1       = unname(b1),
    ep       = unname(se),
    VPA      = unname(vpa),
    VPA_li   = unname(vpa_li),
    VPA_ls   = unname(vpa_ls),
    p_valor  = unname(p),
    rho_AR1  = if (!is.null(sm$rho)) sm$rho[length(sm$rho)] else NA_real_,
    classificacao = dplyr::case_when(
      p >= PARAMS$alfa ~ "Estacionária",
      vpa > 0          ~ "Crescente",
      TRUE             ~ "Decrescente"
    )
  )
}

#' Mann-Kendall + declive de Sen — alternativa não paramétrica
tendencia_mk <- function(dados, var_taxa, rotulo) {
  x <- dados |> dplyr::arrange(ano) |> dplyr::pull(.data[[var_taxa]])
  x <- x[is.finite(x)]
  if (length(x) < 4) {
    return(tibble::tibble(serie = rotulo, tau = NA_real_, p_MK = NA_real_,
                          sen_slope = NA_real_))
  }
  mk  <- trend::mk.test(x)
  sen <- trend::sens.slope(x)
  tibble::tibble(
    serie     = rotulo,
    tau       = unname(mk$estimates["tau"]),
    p_MK      = mk$p.value,
    sen_slope = unname(sen$estimates),
    sen_li    = sen$conf.int[1],
    sen_ls    = sen$conf.int[2]
  )
}

# Aplicar a todas as séries e desfechos
combos <- tidyr::expand_grid(
  reg = unique(serie_full$regiao),
  desfecho = c("incidencia", "positividade", "cob_exame")
)

tab_tendencia <- purrr::pmap_dfr(combos, function(reg, desfecho) {
  d <- dplyr::filter(serie_full, regiao == reg)
  rot <- glue::glue("{reg} — {desfecho}")
  pw <- vpa_prais(d, desfecho, rot)
  mk <- tendencia_mk(d, desfecho, rot)
  dplyr::left_join(pw, mk, by = "serie") |>
    dplyr::mutate(regiao = reg, desfecho = desfecho, .before = 1)
})

tab_tendencia <- tab_tendencia |>
  dplyr::mutate(
    VPA_fmt = ifelse(is.na(VPA), NA_character_,
                     sprintf("%.1f%% (%.1f a %.1f)", VPA, VPA_li, VPA_ls)),
    p_fmt   = fmt_p(p_valor),
    p_MK_fmt = fmt_p(p_MK)
  )

salvar_tabela(tab_tendencia, "03_tendencia_prais_winsten")
print(as.data.frame(dplyr::select(tab_tendencia, regiao, desfecho, VPA_fmt,
                                  p_fmt, classificacao, tau, p_MK_fmt)))


# =============================================================================
# 4. TENDÊNCIA MUNICIPAL (para mapear onde a doença cresce)
# =============================================================================
# Poisson por município com termo linear de ano e offset de população.
# Municípios com menos de 3 anos de dado ou < 5 casos no período são marcados
# como não estimáveis — evita RRs explosivos por instabilidade.

tendencia_municipal <- painel |>
  dplyr::group_by(code_muni) |>
  dplyr::group_modify(function(d, key) {
    if (sum(d$casos_sinan) < 5 || dplyr::n_distinct(d$ano) < 4) {
      return(tibble::tibble(RR_anual = NA_real_, li = NA_real_, ls = NA_real_,
                            p = NA_real_, motivo = "casos insuficientes"))
    }
    fit <- try(glm(casos_sinan ~ ano_c + offset(log(populacao)),
                   family = poisson(), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(tibble::tibble(RR_anual = NA_real_, li = NA_real_, ls = NA_real_,
                            p = NA_real_, motivo = "não convergiu"))
    }
    rr <- tidy_rr_robusto(fit, tipo_vcov = "HC0") |>
      dplyr::filter(termo == "ano_c")
    tibble::tibble(RR_anual = rr$RR, li = rr$RR_li, ls = rr$RR_ls,
                   p = rr$p_valor, motivo = NA_character_)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p_aj = p.adjust(p, method = "fdr"),
    classe_tendencia = dplyr::case_when(
      is.na(RR_anual)              ~ "Não estimável",
      p_aj > PARAMS$alfa           ~ "Estável",
      RR_anual > 1                 ~ "Crescente",
      TRUE                         ~ "Decrescente"
    )
  )

salvar_tabela(tendencia_municipal, "04_tendencia_municipal")
saveRDS(tendencia_municipal,
        file.path(PARAMS$dir_processados, "tendencia_municipal.rds"))

message("\nDistribuição das tendências municipais:")
print(table(tendencia_municipal$classe_tendencia))


# =============================================================================
# 5. ESFORÇO DIAGNÓSTICO — controle de confundimento crítico
# =============================================================================
# Em esquistossomose, queda de positividade pode refletir queda de busca ativa,
# não queda de transmissão. Este bloco quantifica a associação entre cobertura
# de exame e positividade — se forte, a análise de tendência DEVE ser ajustada.

dados_esforco <- painel |>
  dplyr::filter(!is.na(positividade), !is.na(cob_exame), cob_exame > 0)

if (nrow(dados_esforco) < 10) {
  message("\n>>> Sem dados de PCE: a checagem de confundimento por esforço ",
          "diagnóstico NÃO pôde ser feita.\n    Esta é uma lacuna relevante — ",
          "ver limitações do relatório.")
  corr_esforco <- tibble::tibble(ano = integer(), n = integer(),
                                 rho_spearman = numeric(), p = numeric())
} else {
corr_esforco <- dados_esforco |>
  dplyr::group_by(ano) |>
  dplyr::summarise(
    n = dplyr::n(),
    rho_spearman = cor(positividade, cob_exame, method = "spearman",
                       use = "complete.obs"),
    p = cor.test(positividade, cob_exame, method = "spearman",
                 exact = FALSE)$p.value,
    .groups = "drop"
  ) |>
  dplyr::mutate(p_fmt = fmt_p(p))
}

salvar_tabela(corr_esforco, "05_esforco_diagnostico")
message("\nCorrelação positividade x cobertura de exame, por ano:")
print(as.data.frame(corr_esforco))

if (any(abs(corr_esforco$rho_spearman) > 0.3, na.rm = TRUE)) {
  message(
    "\n>>> ALERTA: correlação |rho| > 0,3 entre positividade e cobertura de exame.\n",
    "    A variação temporal da positividade está parcialmente confundida pelo\n",
    "    esforço de busca. Os modelos dos scripts 08-10 já incluem cob_exame\n",
    "    como covariável de ajuste — mantenha e discuta essa limitação."
  )
}


# =============================================================================
# 6. FIGURAS
# =============================================================================

# Só entram no gráfico os indicadores que têm dado. Sem o PCE, painéis de
# positividade e cobertura de exame ficariam vazios — pior que ausentes.
inds <- c("incidencia", "positividade", "cob_exame")
inds <- inds[vapply(inds, function(v)
  any(is.finite(serie_full[[v]]) & serie_full[[v]] > 0), logical(1))]
message("Indicadores com dado para a figura: ", paste(inds, collapse = ", "))

g_serie <- serie_full |>
  tidyr::pivot_longer(dplyr::all_of(inds),
                      names_to = "indicador", values_to = "valor") |>
  dplyr::mutate(indicador = dplyr::recode(indicador,
    incidencia   = "Incidência SINAN (/100 mil hab.)",
    positividade = "Positividade PCE (%)",
    cob_exame    = "Cobertura de exame (% da população)")) |>
  ggplot2::ggplot(ggplot2::aes(ano, valor, colour = regiao, group = regiao)) +
  ggplot2::geom_line(linewidth = .8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~ indicador, scales = "free_y", ncol = 1) +
  ggplot2::scale_x_continuous(breaks = PARAMS$ano_ini:PARAMS$ano_fim) +
  ggplot2::scale_colour_manual(values = c("Norte de Minas" = "#1B7837",
                                          "Bahia" = "#762A83",
                                          "Total" = "grey30")) +
  ggplot2::labs(
    title = "Esquistossomose no corredor endêmico Norte de Minas–Bahia, 2018–2025",
    x = NULL, y = NULL, colour = NULL,
    caption = paste("Fontes: SINAN e PCE/SISPCE (Ministério da Saúde);",
                    "população: estimativas IBGE.")
  ) +
  theme_grafico()

salvar_figura(g_serie, "fig01_series_temporais",
              largura = 8, altura = max(3.2, 3 * length(inds)))

message("\nDescritiva concluída. Próximo: R/05_bayes_empirico.R")
