# =============================================================================
# 08_poisson_robusta.R — Regressão de Poisson com variância robusta
#
# ESTRUTURA:
#   8.1 Justificativa do estimador e definição dos modelos
#   8.2 Análise univariada (triagem, NÃO seleção automática)
#   8.3 Modelo multivariável — incidência SINAN (offset = log pessoa-ano)
#   8.4 Modelo multivariável — positividade PCE (offset = log examinados)
#   8.5 Modelos de painel com GEE (equações de estimação generalizadas)
#   8.6 Diagnósticos: dispersão, colinearidade, resíduos, influência
#   8.7 Autocorrelação espacial dos resíduos -> justifica o script 09
#   8.8 Análises de sensibilidade
#
# POR QUE POISSON MODIFICADA (Zou, 2004) E NÃO LOGÍSTICA/BINOMIAL-LOG:
#   - a positividade da esquistossomose em área endêmica NÃO é rara (>10%),
#     logo o odds ratio superestima o risco relativo;
#   - a regressão binomial com ligação log frequentemente não converge;
#   - o Poisson dá estimativa consistente do log(RR); o erro-padrão do Poisson
#     é conservador/incorreto para proporção, e a variância sanduíche corrige
#     isso sem custo de convergência.
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

library(sandwich)
library(lmtest)
library(MASS)
library(geepack)

agregado_eb <- readRDS(file.path(PARAMS$dir_processados, "agregado_eb.rds"))
painel      <- readRDS(file.path(PARAMS$dir_processados, "painel.rds"))
malha_est   <- readRDS(file.path(PARAMS$dir_processados, "malha_est.rds"))
viz         <- readRDS(file.path(PARAMS$dir_processados, "vizinhanca.rds"))
listw       <- viz$listw


# =============================================================================
# 8.1 PREPARO
# =============================================================================

# Covariáveis efetivamente disponíveis (as ausentes são descartadas com aviso)
cov_disp <- COVARIAVEIS |>
  dplyr::filter(var %in% names(agregado_eb)) |>
  dplyr::filter(purrr::map_lgl(var, ~ mean(!is.na(agregado_eb[[.x]])) > 0.70))

cov_faltantes <- setdiff(COVARIAVEIS$var, cov_disp$var)
if (length(cov_faltantes)) {
  message("Covariáveis excluídas (ausentes ou >30% de faltantes): ",
          paste(cov_faltantes, collapse = ", "))
}

dados_mod <- agregado_eb |>
  dplyr::mutate(cob_exame_periodo = examinados / pessoa_ano * 100) |>
  aplicar_transformacoes(dicionario = cov_disp)

painel_mod <- painel |>
  aplicar_transformacoes(dicionario = cov_disp)


# =============================================================================
# 8.2 ANÁLISE UNIVARIADA
# =============================================================================
# Serve para descrever a associação bruta de cada indicador — NÃO para
# selecionar variáveis por p-valor (prática que enviesa estimativas e
# subestima erros-padrão). A seleção do modelo final é por relevância teórica.

univariada <- function(dados, covs, desfecho, offset_var, escala = "HC0") {
  purrr::map_dfr(covs$var, function(v) {
    d <- dados |> dplyr::filter(!is.na(.data[[v]]), .data[[offset_var]] > 0)
    if (nrow(d) < 20) return(NULL)
    f <- as.formula(glue::glue("{desfecho} ~ {v} + offset(log({offset_var}))"))
    fit <- try(glm(f, family = poisson(), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    tidy_rr_robusto(fit, tipo_vcov = escala) |>
      dplyr::filter(termo == v) |>
      dplyr::mutate(
        covariavel = v,
        rotulo     = covs$rotulo[match(v, covs$var)],
        incremento = rotulo_incremento(v, covs),
        n          = nrow(d),
        .before = 1
      )
  })
}

uni_sinan <- univariada(dados_mod, cov_disp, "casos_sinan", "pessoa_ano")
uni_pce   <- univariada(dplyr::filter(dados_mod, examinados > 0),
                        cov_disp, "positivos", "examinados")

tab_uni <- dplyr::bind_rows(
  dplyr::mutate(uni_sinan, desfecho = "Incidência SINAN"),
  dplyr::mutate(uni_pce,   desfecho = "Positividade PCE")
) |>
  dplyr::mutate(
    RR_IC = fmt_rr(RR, RR_li, RR_ls),
    p_fmt = fmt_p(p_valor)
  ) |>
  dplyr::select(desfecho, covariavel, rotulo, incremento, n, RR_IC, p_fmt,
                RR, RR_li, RR_ls, p_valor)

salvar_tabela(tab_uni, "21_regressao_univariada")
print(as.data.frame(dplyr::select(tab_uni, desfecho, rotulo, incremento, RR_IC, p_fmt)))


# =============================================================================
# 8.3 MODELO MULTIVARIÁVEL — INCIDÊNCIA (SINAN)
# =============================================================================
# Bloco de covariáveis definido por plausibilidade causal, não por p-valor:
#   - condição de vida (IDHM ou IVS — NÃO os dois: colinearidade)
#   - saneamento (água e esgoto — vias de exposição e contaminação)
#   - contexto rural (proxy de contato com coleções hídricas)
#   - ajuste por esforço diagnóstico (cobertura de exame)
#   - efeito de região (MG vs. BA) para captar diferença de sistema de vigilância

montar_formula <- function(desfecho, offset_var, termos) {
  as.formula(glue::glue(
    "{desfecho} ~ {paste(termos, collapse = ' + ')} + offset(log({offset_var}))"
  ))
}

# Bloco a priori, definido por VIA DE TRANSMISSÃO e não por p-valor:
#   pct_agua_superficial      -> contato humano com coleção hídrica (exposição)
#   pct_esgoto_em_corpo_dagua -> contaminação fecal da coleção hídrica
#   pct_agua_rede             -> água encanada como fator de proteção
#   pct_esgoto_adeq           -> esgotamento adequado como fator de proteção
#   dens_demografica          -> proxy de ruralidade (pct_rural indisponível)
#   idhm                      -> condição de vida (se a base tiver sido obtida)
#
# pct_lixo_coletado fica FORA do modelo principal: é proxy geral de
# desenvolvimento, não via de transmissão da esquistossomose, e é fortemente
# correlacionado com os demais indicadores de cobertura. Permanece na análise
# univariada e na de sensibilidade.
termos_base <- intersect(
  c("idhm", "pct_agua_superficial", "pct_esgoto_em_corpo_dagua",
    "pct_agua_rede", "pct_esgoto_adeq", "pct_rural", "dens_demografica"),
  cov_disp$var
)
termos_sinan <- c(termos_base, "regiao")

d_sinan <- dados_mod |>
  dplyr::filter(pessoa_ano > 0) |>
  tidyr::drop_na(dplyr::all_of(termos_base))

message(glue::glue("\nModelo SINAN: n = {nrow(d_sinan)} municípios ",
                   "(de {nrow(dados_mod)}); casos = {sum(d_sinan$casos_sinan)}"))

m_pois_sinan <- glm(montar_formula("casos_sinan", "pessoa_ano", termos_sinan),
                    family = poisson(), data = d_sinan)

# Negativo binomial como comparador — se o alfa for grande, a superdispersão
# é estrutural (excesso de zeros/heterogeneidade não observada)
m_nb_sinan <- try(MASS::glm.nb(
  montar_formula("casos_sinan", "pessoa_ano", termos_sinan), data = d_sinan),
  silent = TRUE)

rr_sinan <- tidy_rr_robusto(m_pois_sinan, tipo_vcov = "HC0") |>
  dplyr::mutate(modelo = "Poisson robusto (HC0)")

if (!inherits(m_nb_sinan, "try-error")) {
  rr_sinan_nb <- broom::tidy(m_nb_sinan, conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::transmute(termo = term, RR = estimate, RR_li = conf.low,
                     RR_ls = conf.high, p_valor = p.value,
                     modelo = "Binomial negativo")
  rr_sinan <- dplyr::bind_rows(rr_sinan, rr_sinan_nb)
}


# =============================================================================
# 8.4 MODELO MULTIVARIÁVEL — POSITIVIDADE (PCE)
# =============================================================================
# Aqui a "regressão de Poisson modificada" no sentido estrito de Zou: desfecho
# é proporção (positivos/examinados), offset = log(examinados), RR interpretado
# como razão de prevalências.

termos_pce <- c(termos_base, "regiao")

d_pce <- dados_mod |>
  dplyr::filter(examinados > 0) |>
  tidyr::drop_na(dplyr::all_of(termos_base))

message(glue::glue("Modelo PCE: n = {nrow(d_pce)} municípios; ",
                   "{sum(d_pce$examinados)} exames, {sum(d_pce$positivos)} positivos"))

# O desfecho de positividade só existe se a base do PCE tiver sido obtida.
TEM_PCE <- nrow(d_pce) >= 30 && sum(d_pce$examinados, na.rm = TRUE) > 0
if (!TEM_PCE) {
  message("\n>>> SEM DADOS DE PCE: todos os modelos de POSITIVIDADE foram ",
          "pulados.\n    Os resultados abaixo se referem apenas à ",
          "incidência notificada (SINAN).")
}

m_pois_pce <- if (TEM_PCE) {
  glm(montar_formula("positivos", "examinados", termos_pce),
      family = poisson(), data = d_pce)
} else NULL

rr_pce <- if (TEM_PCE) {
  tidy_rr_robusto(m_pois_pce, tipo_vcov = "HC0") |>
    dplyr::mutate(modelo = "Poisson modificado (Zou) — RP")
} else NULL

# Comparador binomial com ligação log (o alvo teórico; frequentemente falha)
m_binlog <- if (!TEM_PCE) NULL else try(glm(cbind(positivos, examinados - positivos) ~ .,
                    family = binomial(link = "log"),
                    data = d_pce[, c("positivos", "examinados", termos_pce)],
                    start = c(log(mean(d_pce$positivos / d_pce$examinados)),
                              rep(0, length(termos_pce) +
                                    nlevels(d_pce$regiao) - 2))),
                silent = TRUE)
if (!is.null(m_binlog) && inherits(m_binlog, "try-error")) {
  message("Binomial-log não convergiu (esperado) — justifica o Poisson modificado.")
}


# =============================================================================
# 8.5 MODELOS DE PAINEL — GEE
# =============================================================================
# O painel município-ano viola independência (8 observações correlacionadas por
# município). GEE com estrutura de correlação de trabalho "exchangeable" e
# variância robusta (sanduíche) fornece estimativas populacionais-médias
# válidas mesmo se a estrutura de correlação estiver mal especificada.
#
# Alternativa equivalente: glm + vcovCL agrupado por município (calculado
# abaixo como checagem de consistência).

termos_painel <- c(termos_base, "regiao", "ano_c")

d_painel <- painel_mod |>
  dplyr::filter(populacao > 0) |>
  tidyr::drop_na(dplyr::all_of(termos_base)) |>
  dplyr::arrange(code_muni, ano) |>
  dplyr::mutate(id = as.integer(factor(code_muni)))

m_gee_sinan <- try(geepack::geeglm(
  montar_formula("casos_sinan", "populacao", termos_painel),
  family = poisson(), data = d_painel, id = id,
  corstr = "exchangeable", std.err = "san.se"
), silent = TRUE)

extrair_gee <- function(m, rotulo) {
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$coefficients
  tibble::tibble(
    termo   = rownames(s),
    RR      = exp(s[, "Estimate"]),
    RR_li   = exp(s[, "Estimate"] - 1.96 * s[, "Std.err"]),
    RR_ls   = exp(s[, "Estimate"] + 1.96 * s[, "Std.err"]),
    p_valor = s[, "Pr(>|W|)"],
    modelo  = rotulo
  )
}

rr_gee_sinan <- extrair_gee(m_gee_sinan, "GEE Poisson (exchangeable)")

# Painel do PCE
d_painel_pce <- painel_mod |>
  dplyr::filter(!is.na(examinados), examinados > 0) |>
  tidyr::drop_na(dplyr::all_of(termos_base)) |>
  dplyr::arrange(code_muni, ano) |>
  dplyr::mutate(id = as.integer(factor(code_muni)))

m_gee_pce <- if (!TEM_PCE || nrow(d_painel_pce) < 30) NULL else
  try(geepack::geeglm(
    montar_formula("positivos", "examinados", termos_painel),
    family = poisson(), data = d_painel_pce, id = id,
    corstr = "exchangeable", std.err = "san.se"
  ), silent = TRUE)

rr_gee_pce <- if (is.null(m_gee_pce)) NULL else
  extrair_gee(m_gee_pce, "GEE Poisson modificado — RP")

# Checagem: GLM + erro-padrão agrupado deve dar resultado próximo ao GEE
m_glm_cl <- glm(montar_formula("casos_sinan", "populacao", termos_painel),
                family = poisson(), data = d_painel)
rr_glm_cl <- tidy_rr_robusto(m_glm_cl, tipo_vcov = "cluster",
                             cluster_var = d_painel$code_muni) |>
  dplyr::mutate(modelo = "GLM + EP agrupado por município")


# =============================================================================
# TABELA CONSOLIDADA DE RESULTADOS
# =============================================================================

rotular_termo <- function(x) {
  dic <- setNames(cov_disp$rotulo, cov_disp$var)
  dplyr::case_when(
    x == "(Intercept)"    ~ "Intercepto",
    x == "regiaoBahia"    ~ "Bahia (ref.: Norte de Minas)",
    x == "ano_c"          ~ "Ano (tendência linear)",
    x %in% names(dic)     ~ unname(dic[x]),
    TRUE                  ~ x
  )
}

tab_multi <- dplyr::bind_rows(
  dplyr::mutate(rr_sinan,     desfecho = "Incidência SINAN"),
  dplyr::mutate(rr_gee_sinan, desfecho = "Incidência SINAN (painel)"),
  dplyr::mutate(rr_glm_cl,    desfecho = "Incidência SINAN (painel)"),
  if (!is.null(rr_pce))     dplyr::mutate(rr_pce,     desfecho = "Positividade PCE"),
  if (!is.null(rr_gee_pce)) dplyr::mutate(rr_gee_pce, desfecho = "Positividade PCE (painel)")
) |>
  dplyr::filter(termo != "(Intercept)") |>
  dplyr::mutate(
    rotulo = rotular_termo(termo),
    incremento = rotulo_incremento(termo, cov_disp),
    RR_IC  = fmt_rr(RR, RR_li, RR_ls),
    p_fmt  = fmt_p(p_valor)
  ) |>
  dplyr::select(desfecho, modelo, termo, rotulo, incremento, RR_IC, p_fmt,
                RR, RR_li, RR_ls, p_valor)

salvar_tabela(tab_multi, "22_regressao_multivariavel")
message("\n--- Modelos multivariáveis ---")
print(as.data.frame(dplyr::select(tab_multi, desfecho, modelo, rotulo, RR_IC, p_fmt)))


# =============================================================================
# 8.6 DIAGNÓSTICOS
# =============================================================================

diag <- list(
  dispersao_sinan = diag_dispersao(m_pois_sinan),
  dispersao_pce   = if (TEM_PCE) diag_dispersao(m_pois_pce) else NULL,
  vif_sinan       = diag_vif(m_pois_sinan),
  vif_pce         = if (TEM_PCE) diag_vif(m_pois_pce) else NULL
)

message("\nDispersão — modelo SINAN:")
print(as.data.frame(diag$dispersao_sinan))
if (TEM_PCE) { message("Dispersão — modelo PCE:")
  print(as.data.frame(diag$dispersao_pce)) }
message("\nVIF — modelo SINAN:")
print(as.data.frame(diag$vif_sinan))

salvar_tabela(dplyr::bind_rows(
  dplyr::mutate(diag$dispersao_sinan, modelo = "SINAN"),
  if (TEM_PCE) dplyr::mutate(diag$dispersao_pce, modelo = "PCE")
), "23_diagnostico_dispersao")
salvar_tabela(dplyr::bind_rows(
  dplyr::mutate(diag$vif_sinan, modelo = "SINAN"),
  if (TEM_PCE) dplyr::mutate(diag$vif_pce, modelo = "PCE")
), "24_diagnostico_vif")

if (any(diag$vif_sinan$VIF > 5, na.rm = TRUE)) {
  message("\n>>> VIF > 5 detectado. IDHM e saneamento são fortemente ",
          "correlacionados no Brasil. Considere: (a) manter apenas um bloco, ",
          "(b) usar componentes principais, ou (c) reportar modelos separados ",
          "por bloco (feito na sensibilidade, 8.8).")
}

# Observações influentes (distância de Cook)
cook <- cooks.distance(m_pois_sinan)
infl <- d_sinan |>
  dplyr::mutate(cook = cook) |>
  dplyr::filter(cook > 4 / nrow(d_sinan)) |>
  dplyr::select(code_muni, name_muni, abbrev_state, casos_sinan, pessoa_ano,
                inc_sinan_bruta, cook) |>
  dplyr::arrange(dplyr::desc(cook))
salvar_tabela(infl, "25_observacoes_influentes")
message("\n", nrow(infl), " município(s) com distância de Cook > 4/n.")


# =============================================================================
# 8.7 AUTOCORRELAÇÃO ESPACIAL DOS RESÍDUOS
# =============================================================================
# Se os resíduos do modelo não espacial ainda forem autocorrelacionados, os
# erros-padrão estão subestimados e há estrutura espacial não modelada —
# justificativa formal para o script 09 (SAR/GWR).

testar_residuos_espaciais <- function(modelo, dados, malha, rotulo) {
  res <- residuals(modelo, type = "pearson")
  idx <- match(dados$code_muni, malha$code_muni)
  ok  <- !is.na(idx)
  nb_sub <- spdep::subset.nb(viz$nb, malha$code_muni %in% dados$code_muni)
  lw_sub <- spdep::nb2listw(nb_sub, style = PARAMS$viz_style, zero.policy = TRUE)
  r_ord <- res[order(idx[ok])]
  mg <- moran_global(r_ord, lw_sub, nsim = PARAMS$n_sim)
  dplyr::mutate(mg, modelo = rotulo, .before = 1)
}

res_esp <- dplyr::bind_rows(
  testar_residuos_espaciais(m_pois_sinan, d_sinan, malha_est,
                            "Poisson robusto — SINAN"),
  if (TEM_PCE) testar_residuos_espaciais(m_pois_pce, d_pce, malha_est,
                                         "Poisson modificado — PCE")
) |>
  dplyr::mutate(p_fmt = fmt_p(p_valor),
                conclusao = dplyr::if_else(
                  p_valor <= PARAMS$alfa,
                  "Dependência espacial residual — usar modelo espacial (script 09)",
                  "Sem dependência espacial residual"))

salvar_tabela(res_esp, "26_moran_residuos")
message("\nI de Moran dos resíduos:")
print(as.data.frame(res_esp))


# =============================================================================
# 8.8 ANÁLISES DE SENSIBILIDADE
# =============================================================================

sensibilidade <- list()

# (a) HC3 no lugar de HC0 (melhor desempenho em amostras pequenas)
sensibilidade$hc3 <- tidy_rr_robusto(m_pois_sinan, tipo_vcov = "HC3") |>
  dplyr::mutate(cenario = "Variância HC3")

# (b) IVS no lugar do IDHM (medidas concorrentes de vulnerabilidade)
if ("ivs" %in% cov_disp$var) {
  termos_ivs <- gsub("^idhm$", "ivs", termos_sinan)
  d_ivs <- dados_mod |> tidyr::drop_na(dplyr::all_of(setdiff(termos_ivs, "regiao")))
  m_ivs <- glm(montar_formula("casos_sinan", "pessoa_ano", termos_ivs),
               family = poisson(), data = d_ivs)
  sensibilidade$ivs <- tidy_rr_robusto(m_ivs, "HC0") |>
    dplyr::mutate(cenario = "IVS no lugar do IDHM")
}

# (c) Excluindo capitais/polos assistenciais (Salvador, Feira de Santana,
#     Montes Claros) — testam se o resultado é dirigido por outliers urbanos
polos <- c("2927408", "2910800", "3143302")   # Salvador, Feira de Santana, Montes Claros
d_sem_polos <- dplyr::filter(d_sinan, !code_muni %in% polos)
m_sem_polos <- glm(montar_formula("casos_sinan", "pessoa_ano", termos_sinan),
                   family = poisson(), data = d_sem_polos)
sensibilidade$sem_polos <- tidy_rr_robusto(m_sem_polos, "HC0") |>
  dplyr::mutate(cenario = "Excluindo grandes polos assistenciais")

# (d) Ajuste adicional por cobertura de exame (esforço diagnóstico)
if ("cob_exame_periodo" %in% names(d_sinan)) {
  m_ajust <- glm(montar_formula("casos_sinan", "pessoa_ano",
                                c(termos_sinan, "cob_exame_periodo")),
                 family = poisson(),
                 data = dplyr::filter(d_sinan, is.finite(cob_exame_periodo)))
  sensibilidade$esforco <- tidy_rr_robusto(m_ajust, "HC0") |>
    dplyr::mutate(cenario = "Ajustado por cobertura de exame")
}

# (e) Restrito a municípios com pelo menos 1 caso (evita influência do excesso
#     de zeros estruturais por subnotificação)
d_com_caso <- dplyr::filter(d_sinan, casos_sinan > 0)
m_com_caso <- glm(montar_formula("casos_sinan", "pessoa_ano", termos_sinan),
                  family = poisson(), data = d_com_caso)
sensibilidade$com_caso <- tidy_rr_robusto(m_com_caso, "HC0") |>
  dplyr::mutate(cenario = glue::glue("Apenas municípios com >=1 caso (n={nrow(d_com_caso)})"))

tab_sens <- dplyr::bind_rows(sensibilidade) |>
  dplyr::filter(termo != "(Intercept)") |>
  dplyr::mutate(rotulo = rotular_termo(termo),
                RR_IC = fmt_rr(RR, RR_li, RR_ls),
                p_fmt = fmt_p(p_valor)) |>
  dplyr::select(cenario, rotulo, RR_IC, p_fmt, RR, RR_li, RR_ls, p_valor)

salvar_tabela(tab_sens, "27_analise_sensibilidade")
message("\n--- Análises de sensibilidade ---")
print(as.data.frame(dplyr::select(tab_sens, cenario, rotulo, RR_IC, p_fmt)))


# =============================================================================
# SALVAR MODELOS
# =============================================================================

saveRDS(list(
  poisson_sinan = m_pois_sinan,
  nb_sinan      = m_nb_sinan,
  poisson_pce   = m_pois_pce,
  gee_sinan     = m_gee_sinan,
  gee_pce       = m_gee_pce,
  dados_sinan   = d_sinan,
  dados_pce     = d_pce,
  dados_painel  = d_painel,
  termos_base   = termos_base,
  cov_disp      = cov_disp
), file.path(PARAMS$dir_modelos, "modelos_poisson.rds"))


# =============================================================================
# FIGURA — forest plot
# =============================================================================

g_forest <- tab_multi |>
  dplyr::filter(modelo %in% c("Poisson robusto (HC0)",
                              "Poisson modificado (Zou) — RP")) |>
  dplyr::mutate(rotulo = forcats::fct_reorder(rotulo, RR)) |>
  ggplot2::ggplot(ggplot2::aes(RR, rotulo, colour = desfecho)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  ggplot2::geom_pointrange(ggplot2::aes(xmin = RR_li, xmax = RR_ls),
                           position = ggplot2::position_dodge(width = .5)) +
  ggplot2::scale_x_continuous(trans = "log10") +
  ggplot2::scale_colour_manual(values = c("Incidência SINAN" = "#B2182B",
                                          "Positividade PCE" = "#2166AC")) +
  ggplot2::labs(
    title = "Determinantes socioeconômicos e de saneamento",
    subtitle = "Razões de risco/prevalência com IC95% por variância robusta",
    x = "RR / RP (escala logarítmica)", y = NULL, colour = NULL,
    caption = "Poisson com variância sanduíche (Zou, 2004). Contínuas padronizadas ou por 10 p.p."
  ) +
  theme_grafico()

salvar_figura(g_forest, "fig04_forest_plot", largura = 9, altura = 6)

message("\nRegressão concluída. Próximo: R/09_espacial_sar_gwr.R")
