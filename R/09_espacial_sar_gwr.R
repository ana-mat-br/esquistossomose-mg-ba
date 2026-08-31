# =============================================================================
# 09_espacial_sar_gwr.R — Modelos de regressão espacial
#
#   9.1  OLS de referência + testes do multiplicador de Lagrange
#   9.2  Modelos globais: SAR (lag), SEM (erro), SDM (Durbin), SAC
#   9.3  Impactos diretos, indiretos (transbordamento) e totais
#   9.4  Filtragem espacial por autovetores de Moran (ESF) aplicada ao Poisson
#        — alternativa que preserva a natureza de contagem do desfecho
#   9.5  GWR / GWPR: heterogeneidade espacial dos coeficientes
#   9.6  Teste de não estacionariedade e colinearidade local
#   9.7  Comparação e seleção de modelos
#
# LÓGICA DE ESCOLHA (siga esta ordem no artigo):
#   1. Resíduos do modelo não espacial são autocorrelacionados? (script 08)
#   2. Se sim, LM robusto indica lag ou erro?
#   3. Há transbordamento substantivo? -> reportar impactos do SAR/SDM
#   4. O efeito das covariáveis varia no espaço? -> GWR
#   Modelos globais (SAR/SEM) e locais (GWR) respondem perguntas DIFERENTES;
#   não são concorrentes e devem ser reportados juntos.
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

library(spatialreg)
library(GWmodel)

agregado_eb <- readRDS(file.path(PARAMS$dir_processados, "agregado_eb.rds"))
malha_est   <- readRDS(file.path(PARAMS$dir_processados, "malha_est.rds"))
viz         <- readRDS(file.path(PARAMS$dir_processados, "vizinhanca.rds"))
mods_pois   <- readRDS(file.path(PARAMS$dir_modelos, "modelos_poisson.rds"))

listw       <- viz$listw
termos_base <- mods_pois$termos_base
cov_disp    <- mods_pois$cov_disp

stopifnot(identical(agregado_eb$code_muni, malha_est$code_muni))


# =============================================================================
# 9.0 PREPARO
# =============================================================================
# Modelos SAR/SEM/GWR gaussianos exigem resposta contínua aproximadamente
# normal. Usamos log da taxa bayesiana empírica local — que já corrige o
# pequeno número (script 05) e estabiliza a variância.

dados_sp <- agregado_eb |>
  aplicar_transformacoes(dicionario = cov_disp) |>
  dplyr::mutate(
    y_inc   = log(inc_sinan_ebl + 0.5),
    y_posit = log(posit_ebl + 0.5),
    bahia   = as.integer(regiao == "Bahia")
  )

termos_sp <- c(termos_base, "bahia")
f_sp <- as.formula(paste("y_inc ~", paste(termos_sp, collapse = " + ")))

# Casos completos — SAR/SEM não aceitam NA e a subsetagem deve ser
# acompanhada da subsetagem da matriz de pesos.
completos <- complete.cases(dados_sp[, c("y_inc", termos_sp)])
message(glue::glue("Casos completos para regressão espacial: {sum(completos)} de {nrow(dados_sp)}"))

d_sp   <- dados_sp[completos, ]
sf_sp  <- malha_est[completos, ]
nb_sp  <- spdep::subset.nb(viz$nb, completos)
lw_sp  <- spdep::nb2listw(nb_sp, style = PARAMS$viz_style, zero.policy = TRUE)

# Normalidade da resposta transformada
sw <- shapiro.test(d_sp$y_inc)
message(glue::glue("Shapiro-Wilk de log(taxa EB): W = {round(sw$statistic,3)}, ",
                   "p = {signif(sw$p.value,3)}"))


# =============================================================================
# 9.1 OLS + TESTES DO MULTIPLICADOR DE LAGRANGE
# =============================================================================

m_ols <- lm(f_sp, data = d_sp)
summary(m_ols)

moran_res_ols <- spdep::lm.morantest(m_ols, lw_sp, zero.policy = TRUE)
message(glue::glue("\nI de Moran dos resíduos OLS: {round(moran_res_ols$estimate[1],4)} ",
                   "(p = {signif(moran_res_ols$p.value,4)})"))

lm_tests <- teste_ml_espacial(m_ols, lw_sp)
salvar_tabela(lm_tests, "28_testes_lagrange")
message("\nTestes do multiplicador de Lagrange:")
print(as.data.frame(lm_tests))

# Regra de decisão de Anselin (1988; Anselin & Rey, 2014)
decidir_modelo <- function(lm_tests) {
  g <- function(padroes) {
    i <- which(grepl(paste(padroes, collapse = "|"), lm_tests$teste, ignore.case = TRUE))
    if (length(i)) lm_tests$p_valor[i[1]] else NA_real_
  }
  p_err   <- g(c("^LMerr", "^RSerr"))
  p_lag   <- g(c("^LMlag", "^RSlag"))
  p_rerr  <- g(c("^RLMerr", "adjRSerr"))
  p_rlag  <- g(c("^RLMlag", "adjRSlag"))
  a <- PARAMS$alfa
  dplyr::case_when(
    p_err > a & p_lag > a                     ~ "OLS adequado (sem dependência espacial)",
    p_lag <= a & p_err > a                    ~ "SAR (lag espacial)",
    p_err <= a & p_lag > a                    ~ "SEM (erro espacial)",
    p_rlag <= a & (is.na(p_rerr) | p_rerr > a) ~ "SAR (lag) — pelo LM robusto",
    p_rerr <= a & (is.na(p_rlag) | p_rlag > a) ~ "SEM (erro) — pelo LM robusto",
    TRUE ~ "Ambos significativos — comparar SAR, SEM, SDM e SAC por AIC/LR"
  )
}
decisao <- decidir_modelo(lm_tests)
message("\n>>> Decisão pelos testes LM: ", decisao)


# =============================================================================
# 9.2 MODELOS ESPACIAIS GLOBAIS
# =============================================================================

m_sar <- spatialreg::lagsarlm(f_sp, data = d_sp, listw = lw_sp,
                              zero.policy = TRUE, method = "eigen")
m_sem <- spatialreg::errorsarlm(f_sp, data = d_sp, listw = lw_sp,
                                zero.policy = TRUE, method = "eigen")
m_sdm <- spatialreg::lagsarlm(f_sp, data = d_sp, listw = lw_sp, type = "mixed",
                              zero.policy = TRUE, method = "eigen")
m_sac <- try(spatialreg::sacsarlm(f_sp, data = d_sp, listw = lw_sp,
                                  zero.policy = TRUE), silent = TRUE)
m_slx <- spatialreg::lmSLX(f_sp, data = d_sp, listw = lw_sp, zero.policy = TRUE)

comparar_modelos <- function(lst) {
  purrr::imap_dfr(lst, function(m, nome) {
    if (inherits(m, "try-error")) return(NULL)
    ll <- as.numeric(logLik(m))
    tibble::tibble(
      modelo = nome,
      AIC = AIC(m), BIC = BIC(m), logLik = ll,
      rho    = if (!is.null(m$rho))    unname(m$rho)    else NA_real_,
      lambda = if (!is.null(m$lambda)) unname(m$lambda) else NA_real_
    )
  })
}

tab_comp <- comparar_modelos(list(
  OLS = m_ols, SLX = m_slx, SAR = m_sar, SEM = m_sem, SDM = m_sdm, SAC = m_sac
)) |>
  dplyr::arrange(AIC) |>
  dplyr::mutate(delta_AIC = AIC - min(AIC))

salvar_tabela(tab_comp, "29_comparacao_modelos_espaciais")
message("\nComparação de modelos espaciais:")
print(as.data.frame(tab_comp))

# Testes de razão de verossimilhança contra o OLS
lr <- function(m, nome) {
  if (inherits(m, "try-error")) return(NULL)
  t <- spatialreg::LR.Sarlm(m, m_ols)
  tibble::tibble(modelo = nome, LR = unname(t$statistic),
                 gl = unname(t$parameter), p_valor = t$p.value)
}
tab_lr <- dplyr::bind_rows(lr(m_sar, "SAR vs OLS"), lr(m_sem, "SEM vs OLS"),
                           lr(m_sdm, "SDM vs OLS")) |>
  dplyr::mutate(p_fmt = fmt_p(p_valor))
salvar_tabela(tab_lr, "30_teste_razao_verossimilhanca")

# Teste de Hausman: SEM consistente vs. OLS
hz <- try(spatialreg::Hausman.test(m_sem), silent = TRUE)
if (!inherits(hz, "try-error")) {
  message(glue::glue("Hausman (SEM): chi2 = {round(hz$statistic,2)}, ",
                     "p = {signif(hz$p.value,4)}"))
}

# Coeficientes do melhor modelo (menor AIC entre os espaciais)
melhor_nome <- tab_comp$modelo[tab_comp$modelo != "OLS"][1]
melhor <- list(SAR = m_sar, SEM = m_sem, SDM = m_sdm, SAC = m_sac,
               SLX = m_slx)[[melhor_nome]]
message("\n>>> Modelo espacial selecionado por AIC: ", melhor_nome)

extrair_coef_sarlm <- function(m, nome) {
  s <- summary(m)
  ct <- s$Coef
  tibble::tibble(
    modelo = nome,
    termo  = rownames(ct),
    beta   = ct[, 1], ep = ct[, 2], z = ct[, 3], p_valor = ct[, 4],
    # Em escala log-taxa, exp(beta) é razão de taxas
    razao_taxas = exp(ct[, 1]),
    rt_li = exp(ct[, 1] - 1.96 * ct[, 2]),
    rt_ls = exp(ct[, 1] + 1.96 * ct[, 2])
  )
}

tab_coef_sp <- dplyr::bind_rows(
  extrair_coef_sarlm(m_sar, "SAR"),
  extrair_coef_sarlm(m_sem, "SEM"),
  extrair_coef_sarlm(m_sdm, "SDM")
) |>
  dplyr::mutate(RT_IC = fmt_rr(razao_taxas, rt_li, rt_ls), p_fmt = fmt_p(p_valor))

salvar_tabela(tab_coef_sp, "31_coeficientes_modelos_espaciais")


# =============================================================================
# 9.3 IMPACTOS DIRETOS, INDIRETOS E TOTAIS
# =============================================================================
# Em modelos com defasagem da resposta (SAR/SDM/SAC), o coeficiente NÃO é o
# efeito marginal: parte do efeito se propaga pela vizinhança e retorna. A
# decomposição de LeSage & Pace (2009) é obrigatória para interpretação.
#   - direto   : efeito no próprio município (inclui retroalimentação)
#   - indireto : transbordamento para os vizinhos  <- achado epidemiológico
#                relevante: saneamento de um município protege os vizinhos
#   - total    : soma

if (melhor_nome %in% c("SAR", "SDM", "SAC", "SLX")) {

  # `listw =` usa a decomposição exata (adequada para N de algumas centenas).
  # Se o N crescer muito, troque por traços aproximados: trW(listw = lw_sp,
  # type = "MC").
  imp <- try(spatialreg::impacts(melhor, listw = lw_sp, R = 1000),
             silent = TRUE)
  if (inherits(imp, "try-error")) {
    W_tr <- spatialreg::trW(listw = lw_sp, type = "MC")
    imp <- spatialreg::impacts(melhor, tr = W_tr, R = 1000)
  }
  s_imp <- summary(imp, zstats = TRUE, short = TRUE)

  # `impacts` nomeia as linhas como "x1 dy/dx" — removemos o sufixo para
  # permitir a junção com o dicionário de covariáveis
  tab_impactos <- tibble::tibble(
    termo    = sub("\\s*dy/dx\\s*$", "", rownames(s_imp$direct_sum$statistics)),
    direto   = s_imp$direct_sum$statistics[, "Mean"],
    indireto = s_imp$indirect_sum$statistics[, "Mean"],
    total    = s_imp$total_sum$statistics[, "Mean"],
    p_direto   = s_imp$pzmat[, "Direct"],
    p_indireto = s_imp$pzmat[, "Indirect"],
    p_total    = s_imp$pzmat[, "Total"]
  ) |>
    dplyr::mutate(
      dplyr::across(c(p_direto, p_indireto, p_total), fmt_p, .names = "{.col}_fmt"),
      modelo = melhor_nome
    )

  salvar_tabela(tab_impactos, "32_impactos_espaciais")
  message("\nDecomposição de impactos (", melhor_nome, "):")
  print(as.data.frame(tab_impactos))

  sig_indireto <- tab_impactos$termo[tab_impactos$p_indireto <= PARAMS$alfa]
  if (length(sig_indireto)) {
    message("\n>>> Transbordamento significativo em: ",
            paste(sig_indireto, collapse = ", "),
            "\n    Interpretação: o indicador nos municípios vizinhos afeta o risco ",
            "local — evidência a favor de intervenção em escala regional, não municipal.")
  }
}


# =============================================================================
# 9.4 FILTRAGEM ESPACIAL (ESF) COM POISSON
# =============================================================================
# Alternativa metodologicamente mais coerente que SAR/SEM quando o desfecho é
# contagem: mantém-se a regressão de Poisson (com offset e variância robusta,
# como no script 08) e adiciona-se um conjunto de autovetores de Moran que
# absorvem a estrutura espacial residual. Os RRs continuam interpretáveis.

d_esf <- mods_pois$dados_sinan
idx   <- match(d_esf$code_muni, malha_est$code_muni)
nb_esf <- spdep::subset.nb(viz$nb, malha_est$code_muni %in% d_esf$code_muni)
lw_esf <- spdep::nb2listw(nb_esf, style = PARAMS$viz_style, zero.policy = TRUE)
d_esf  <- d_esf[order(idx), ]

f_esf <- as.formula(paste("casos_sinan ~",
                          paste(c(termos_base, "regiao"), collapse = " + "),
                          "+ offset(log(pessoa_ano))"))

me_fit <- try(spatialreg::ME(f_esf, data = d_esf, family = poisson(),
                             listw = lw_esf, alpha = 0.10, nsim = 199,
                             verbose = FALSE), silent = TRUE)

if (!inherits(me_fit, "try-error") && !is.null(me_fit$vectors) &&
    ncol(as.matrix(me_fit$vectors)) > 0) {

  ev <- as.data.frame(as.matrix(me_fit$vectors))
  names(ev) <- paste0("EV", seq_len(ncol(ev)))
  d_esf2 <- dplyr::bind_cols(d_esf, ev)

  f_esf2 <- as.formula(paste("casos_sinan ~",
                             paste(c(termos_base, "regiao", names(ev)), collapse = " + "),
                             "+ offset(log(pessoa_ano))"))
  m_esf <- glm(f_esf2, family = poisson(), data = d_esf2)

  rr_esf <- tidy_rr_robusto(m_esf, "HC0") |>
    dplyr::filter(!grepl("^EV|Intercept", termo)) |>
    dplyr::mutate(modelo = glue::glue("Poisson + ESF ({ncol(ev)} autovetores)"),
                  RR_IC = fmt_rr(RR, RR_li, RR_ls), p_fmt = fmt_p(p_valor))

  # Os resíduos ainda são autocorrelacionados?
  moran_esf <- moran_global(residuals(m_esf, type = "pearson"), lw_esf,
                            nsim = PARAMS$n_sim)
  message(glue::glue("\nESF: {ncol(ev)} autovetores selecionados. ",
                     "I de Moran residual = {round(moran_esf$I,4)} ",
                     "(p = {signif(moran_esf$p_valor,3)})"))

  salvar_tabela(rr_esf, "33_poisson_filtragem_espacial")
  print(as.data.frame(dplyr::select(rr_esf, termo, RR_IC, p_fmt)))
} else {
  message("Filtragem espacial: nenhum autovetor selecionado (ou ME falhou). ",
          "Isso indica que as covariáveis já absorvem a estrutura espacial.")
  rr_esf <- NULL
}


# =============================================================================
# 9.5 GWR / GWPR — HETEROGENEIDADE ESPACIAL DOS COEFICIENTES
# =============================================================================
# Pergunta respondida: o efeito do saneamento sobre a esquistossomose é o mesmo
# em todo o corredor endêmico, ou é mais forte em determinadas sub-regiões?
# Um coeficiente global médio pode mascarar sinais opostos entre MG e BA.

sp_gwr <- sf::as_Spatial(sf::st_transform(
  dplyr::bind_cols(sf::st_geometry(sf_sp) |> sf::st_sf(),
                   sf::st_drop_geometry(sf_sp)[, c("code_muni", "name_muni")],
                   d_sp[, c("y_inc", termos_sp, "casos_sinan", "pessoa_ano")]),
  5880))   # projeção métrica: bandas em metros/km, não graus

# A dummy de região (`bahia`) é EXCLUÍDA da GWR de propósito: dentro de janelas
# locais inteiramente contidas em uma UF ela é constante, tornando X'WX singular
# ("inv(): matrix is singular"). Além disso é conceitualmente redundante — a GWR
# já modela variação geográfica de forma contínua.
termos_gwr <- setdiff(termos_sp, "bahia")
f_gwr <- as.formula(paste("y_inc ~", paste(termos_gwr, collapse = " + ")))
message("GWR — covariáveis: ", paste(termos_gwr, collapse = ", "))

# Largura de banda ADAPTATIVA (nº de vizinhos): apropriada quando a densidade
# de municípios é heterogênea — caso da Bahia (municípios grandes no oeste,
# pequenos no Recôncavo).
# Matriz de distâncias pré-calculada: acelera as chamadas seguintes e é
# EXIGIDA por bw.ggwr/ggwr.basic em algumas versões do GWmodel (que não têm
# valor padrão para `dMat` e falham com "argumento 'dMat' ausente").
dMat_gwr <- GWmodel::gw.dist(dp.locat = sp::coordinates(sp_gwr))

bw_gwr <- GWmodel::bw.gwr(f_gwr, data = sp_gwr, approach = "AICc",
                          kernel = "bisquare", adaptive = TRUE, dMat = dMat_gwr)
message(glue::glue("\nLargura de banda GWR (adaptativa, AICc): {bw_gwr} vizinhos"))

m_gwr <- GWmodel::gwr.basic(f_gwr, data = sp_gwr, bw = bw_gwr,
                            kernel = "bisquare", adaptive = TRUE,
                            F123.test = TRUE, dMat = dMat_gwr)
print(m_gwr)

gwr_sdf <- as.data.frame(m_gwr$SDF)

# Variação espacial dos coeficientes
resumo_gwr <- purrr::map_dfr(termos_gwr, function(v) {
  if (!v %in% names(gwr_sdf)) return(NULL)
  b <- gwr_sdf[[v]]
  tse <- gwr_sdf[[paste0(v, "_TV")]]
  tibble::tibble(
    covariavel = v,
    rotulo = dplyr::coalesce(cov_disp$rotulo[match(v, cov_disp$var)], v),
    beta_min = min(b), beta_q1 = quantile(b, .25), beta_mediana = median(b),
    beta_q3 = quantile(b, .75), beta_max = max(b),
    amplitude_IQR = quantile(b, .75) - quantile(b, .25),
    pct_t_maior_196 = mean(abs(tse) > 1.96) * 100,
    troca_de_sinal = min(b) < 0 & max(b) > 0
  )
})
salvar_tabela(resumo_gwr, "34_coeficientes_gwr_resumo")
message("\nVariabilidade espacial dos coeficientes (GWR):")
print(as.data.frame(resumo_gwr))

# Superfície de coeficientes para mapeamento
gwr_sf <- sf_sp |>
  dplyr::select(code_muni, name_muni) |>
  dplyr::bind_cols(
    gwr_sdf[, intersect(c(termos_sp, paste0(termos_sp, "_TV"),
                          "Local_R2", "yhat", "residual"), names(gwr_sdf)),
            drop = FALSE]
  )
sf::st_write(gwr_sf, file.path(PARAMS$dir_processados, "gwr_superficie.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# --- GWPR (Poisson geograficamente ponderado) — tentativa ------------------
# Preserva a natureza de contagem. Nem toda versão do GWmodel aceita offset em
# ggwr.basic; se falhar, a GWR gaussiana sobre log-taxa EB permanece como
# análise principal e isso deve ser declarado como limitação.
sp_gwr$log_pop <- log(sp_gwr$pessoa_ano)
f_gwpr <- as.formula(paste("casos_sinan ~", paste(termos_sp, collapse = " + "),
                           "+ offset(log_pop)"))

m_gwpr <- try({
  bw_p <- GWmodel::bw.ggwr(f_gwpr, data = sp_gwr, family = "poisson",
                           approach = "AICc", kernel = "bisquare",
                           adaptive = TRUE, dMat = dMat_gwr)
  GWmodel::ggwr.basic(f_gwpr, data = sp_gwr, bw = bw_p, family = "poisson",
                      kernel = "bisquare", adaptive = TRUE, dMat = dMat_gwr)
}, silent = TRUE)

if (!inherits(m_gwpr, "try-error")) {
  message("GWPR estimado com sucesso.")
  gwpr_sdf <- as.data.frame(m_gwpr$SDF)
  gwpr_sf <- sf_sp |>
    dplyr::select(code_muni, name_muni) |>
    dplyr::bind_cols(gwpr_sdf[, intersect(termos_sp, names(gwpr_sdf)), drop = FALSE])
  sf::st_write(gwpr_sf, file.path(PARAMS$dir_processados, "gwpr_superficie.gpkg"),
               delete_dsn = TRUE, quiet = TRUE)
} else {
  message("GWPR não estimado (",
          substr(as.character(attr(m_gwpr, "condition")$message), 1, 120),
          "). Usando GWR gaussiana sobre log-taxa EB.")
}


# =============================================================================
# 9.6 NÃO ESTACIONARIEDADE E COLINEARIDADE LOCAL
# =============================================================================

# Teste de Monte Carlo: os coeficientes variam mais do que o esperado ao acaso?
mc_gwr <- try(GWmodel::gwr.montecarlo(f_gwr, data = sp_gwr, nsims = 99,
                                      kernel = "bisquare", adaptive = TRUE,
                                      bw = bw_gwr, dMat = dMat_gwr), silent = TRUE)
if (!inherits(mc_gwr, "try-error")) {
  # gwr.montecarlo devolve uma MATRIZ (termos nas linhas), não um vetor nomeado
  tab_mc <- (if (is.null(dim(mc_gwr))) {
    tibble::tibble(covariavel = names(mc_gwr), p_valor = as.numeric(mc_gwr))
  } else {
    tibble::tibble(covariavel = rownames(mc_gwr), p_valor = as.numeric(mc_gwr[, 1]))
  }) |>
    dplyr::mutate(p_fmt = fmt_p(p_valor),
                  conclusao = dplyr::if_else(p_valor <= PARAMS$alfa,
                                             "Efeito não estacionário no espaço",
                                             "Efeito estacionário"))
  salvar_tabela(tab_mc, "35_teste_nao_estacionariedade")
  message("\nTeste de não estacionariedade (Monte Carlo, 99 simulações):")
  print(as.data.frame(tab_mc))
}

# Colinearidade LOCAL: a GWR pode ser instável mesmo quando o VIF global é
# aceitável, porque em janelas pequenas as covariáveis podem ser quase
# colineares. Número de condição local > 30 é sinal de alerta.
collin <- try(GWmodel::gwr.collin.diagno(f_gwr, data = sp_gwr, bw = bw_gwr,
                                         kernel = "bisquare", adaptive = TRUE,
                                         dMat = dMat_gwr),
              silent = TRUE)
if (!inherits(collin, "try-error")) {
  cn <- collin$SDF$local_CN
  message(glue::glue("Número de condição local: mediana = {round(median(cn),1)}, ",
                     "máx = {round(max(cn),1)}; ",
                     "{round(mean(cn > 30)*100,1)}% dos municípios > 30"))
  if (mean(cn > 30) > 0.10) {
    message(">>> Colinearidade local relevante. Considere reduzir o número de ",
            "covariáveis na GWR ou usar GWR ridge (GWmodel::gwr.lcr).")
  }
}


# =============================================================================
# 9.7 COMPARAÇÃO FINAL
# =============================================================================

comparacao_final <- tibble::tibble(
  modelo = c("OLS", "SAR", "SEM", "SDM", "GWR"),
  AIC = c(AIC(m_ols), AIC(m_sar), AIC(m_sem), AIC(m_sdm), m_gwr$GW.diagnostic$AICc),
  R2_ajustado = c(summary(m_ols)$adj.r.squared, NA, NA, NA,
                  m_gwr$GW.diagnostic$gwR2.adj)
) |>
  dplyr::mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE)) |>
  dplyr::arrange(AIC)

salvar_tabela(comparacao_final, "36_comparacao_final_modelos")
message("\n--- Comparação final ---")
print(as.data.frame(comparacao_final))

saveRDS(list(
  ols = m_ols, sar = m_sar, sem = m_sem, sdm = m_sdm, sac = m_sac, slx = m_slx,
  gwr = m_gwr, gwpr = m_gwpr, bw_gwr = bw_gwr,
  esf = if (exists("m_esf")) m_esf else NULL,
  dados = d_sp, malha = sf_sp, listw = lw_sp, decisao_lm = decisao
), file.path(PARAMS$dir_modelos, "modelos_espaciais.rds"))

message("\nModelagem espacial concluída. Próximo: R/10_bym2_inla.R (opcional) ",
        "ou R/11_mapas_tabelas.R")
