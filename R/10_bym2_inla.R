# =============================================================================
# 10_bym2_inla.R — Modelo bayesiano hierárquico espaço-temporal (OPCIONAL)
#
# POR QUE ESTE SCRIPT EXISTE:
# Os estimadores bayesianos empíricos (script 05) suavizam taxas, mas NÃO
# quantificam incerteza nem permitem ajuste por covariáveis simultaneamente.
# Os modelos SAR/GWR (script 09) tratam a estrutura espacial, mas não a
# temporal. O modelo BYM2 com interação espaço-tempo resolve os dois de uma vez:
#
#   log(RR_it) = alfa + X_it*beta + (u_i + v_i) + gama_t + fi_t + delta_it
#
#   u_i     campo espacial estruturado (ICAR) — risco compartilhado com vizinhos
#   v_i     heterogeneidade não estruturada específica do município
#   gama_t  tendência temporal estruturada (RW1)
#   fi_t    ruído temporal não estruturado
#   delta_it interação espaço-tempo (Knorr-Held, tipos I-IV)
#
# A parametrização BYM2 (Riebler et al., 2016) separa a variância total de sua
# repartição entre componente estruturado e não estruturado (parâmetro phi),
# com priors PC (penalized complexity) que evitam a superajuste típico dos
# priors gamma "não informativos" do BYM clássico.
#
# REQUISITO: pacote INLA (não está no CRAN — ver R/00_setup.R).
# Se o INLA não estiver instalado, este script é ignorado sem quebrar o pipeline.
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

if (!requireNamespace("INLA", quietly = TRUE)) {
  message("INLA não instalado — script 10 ignorado.\n",
          "Para instalar:\n",
          '  install.packages("INLA", repos = c(getOption("repos"),\n',
          '    INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)')
  knitr_exit <- TRUE
} else {
  knitr_exit <- FALSE
  library(INLA)
}

if (!knitr_exit) {

painel      <- readRDS(file.path(PARAMS$dir_processados, "painel.rds"))
malha_est   <- readRDS(file.path(PARAMS$dir_processados, "malha_est.rds"))
viz         <- readRDS(file.path(PARAMS$dir_processados, "vizinhanca.rds"))
mods_pois   <- readRDS(file.path(PARAMS$dir_modelos, "modelos_poisson.rds"))
termos_base <- mods_pois$termos_base
cov_disp    <- mods_pois$cov_disp


# =============================================================================
# 10.1 GRAFO DE VIZINHANÇA
# =============================================================================

arq_grafo <- file.path(PARAMS$dir_processados, "mapa.graph")
spdep::nb2INLA(arq_grafo, viz$nb)
grafo <- INLA::inla.read.graph(filename = arq_grafo)
message("Grafo INLA: ", grafo$n, " nós, ", sum(grafo$nnbs), " arestas.")


# =============================================================================
# 10.2 DADOS
# =============================================================================
# Padronização indireta: E_it = pop_it * taxa global do período.
# O modelo estima o RR relativo a esse padrão.

d <- painel |>
  aplicar_transformacoes(dicionario = cov_disp) |>
  dplyr::arrange(ano, code_muni)

taxa_global <- sum(d$casos_sinan, na.rm = TRUE) / sum(d$populacao, na.rm = TRUE)

d <- d |>
  dplyr::mutate(
    esperado = populacao * taxa_global,
    id_esp   = match(code_muni, malha_est$code_muni),   # 1..N, alinhado ao grafo
    id_esp2  = id_esp,
    id_tempo = as.integer(factor(ano)),
    id_tempo2 = id_tempo,
    id_st    = seq_len(dplyr::n()),
    bahia    = as.integer(regiao == "Bahia")
  ) |>
  dplyr::filter(!is.na(esperado), esperado > 0)

stopifnot(max(d$id_esp) <= grafo$n)


# =============================================================================
# 10.3 PRIORS PC
# =============================================================================
# PC prior para precisão: P(sigma > U) = a. Com U=1 e a=0.01 assume-se que é
# improvável que o desvio-padrão do efeito aleatório exceda 1 na escala log —
# ou seja, RR além de ~e^2 = 7,4 é considerado implausível a priori. Isso é
# fracamente informativo e adequado para dados de doença.
# PC prior para phi (BYM2): P(phi < 0.5) = 2/3 favorece levemente a
# heterogeneidade não estruturada, evitando forçar suavização espacial.

pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc", param = c(0.5, 2/3)))

controle <- list(
  compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE,
                 return.marginals.predictor = TRUE),
  predictor = list(compute = TRUE, link = 1)
)


# =============================================================================
# 10.4 MODELOS
# =============================================================================
# Sequência de complexidade crescente. Selecionar por WAIC/DIC, penalizando
# complexidade; se dois modelos empatarem, prefira o mais simples.

termos_fixos <- paste(c(termos_base, "bahia"), collapse = " + ")

f_espacial <- as.formula(glue::glue(
  "casos_sinan ~ {termos_fixos} +
   f(id_esp, model = 'bym2', graph = grafo, scale.model = TRUE,
     constr = TRUE, hyper = pc_bym2)"
))

f_st_aditivo <- as.formula(glue::glue(
  "casos_sinan ~ {termos_fixos} +
   f(id_esp, model = 'bym2', graph = grafo, scale.model = TRUE,
     constr = TRUE, hyper = pc_bym2) +
   f(id_tempo, model = 'rw1', scale.model = TRUE, hyper = pc_prec) +
   f(id_tempo2, model = 'iid', hyper = pc_prec)"
))

# Interação tipo I (Knorr-Held): espaço e tempo não estruturados —
# capta discrepâncias município-ano isoladas (surtos pontuais)
f_st_int1 <- as.formula(glue::glue(
  "casos_sinan ~ {termos_fixos} +
   f(id_esp, model = 'bym2', graph = grafo, scale.model = TRUE,
     constr = TRUE, hyper = pc_bym2) +
   f(id_tempo, model = 'rw1', scale.model = TRUE, hyper = pc_prec) +
   f(id_tempo2, model = 'iid', hyper = pc_prec) +
   f(id_st, model = 'iid', hyper = pc_prec)"
))

# Interação tipo IV: espaço estruturado x tempo estruturado — permite que cada
# município tenha sua PRÓPRIA tendência temporal, correlacionada com a dos
# vizinhos. É o modelo mais informativo para "dinâmica espaço-temporal", mas o
# mais caro computacionalmente.
n_esp   <- grafo$n
n_tempo <- dplyr::n_distinct(d$id_tempo)

# Matriz de estrutura espacial (ICAR): R = diag(grau) - W
W_biny  <- spdep::nb2mat(viz$nb, style = "B", zero.policy = TRUE)
R_esp   <- diag(rowSums(W_biny)) - W_biny

# Matriz de estrutura temporal (RW1): R = D'D, com D a matriz de 1ª diferença
D_rw1   <- diff(diag(n_tempo))
R_tempo <- crossprod(D_rw1)

# Escalonamento (Sørbye & Rue, 2014): torna a variância generalizada unitária,
# permitindo que o mesmo prior PC signifique a mesma coisa nos dois componentes
R_esp   <- INLA::inla.scale.model(R_esp,
             constr = list(A = matrix(1, 1, n_esp), e = 0))
R_tempo <- INLA::inla.scale.model(R_tempo,
             constr = list(A = matrix(1, 1, n_tempo), e = 0))

R_int4 <- kronecker(R_esp, R_tempo)

# Restrições de soma zero por área e por tempo (Goicoa et al., 2018).
# São n_esp + n_tempo linhas com posto n_esp + n_tempo - 1 — a redundância é
# esperada. O INLA emite avisos "Matrix AA^t is numerical singular, remove
# singularity and move on"; são BENIGNOS e o ajuste prossegue corretamente.
restricoes_t4 <- list(
  A = rbind(kronecker(matrix(1, 1, n_esp), diag(n_tempo)),
            kronecker(diag(n_esp), matrix(1, 1, n_tempo))),
  e = rep(0, n_esp + n_tempo)
)

d$id_st4 <- (d$id_esp - 1) * n_tempo + d$id_tempo

f_st_int4 <- as.formula(glue::glue(
  "casos_sinan ~ {termos_fixos} +
   f(id_esp, model = 'bym2', graph = grafo, scale.model = TRUE,
     constr = TRUE, hyper = pc_bym2) +
   f(id_tempo, model = 'rw1', scale.model = TRUE, hyper = pc_prec) +
   f(id_tempo2, model = 'iid', hyper = pc_prec) +
   f(id_st4, model = 'generic0', Cmatrix = R_int4, rankdef = {n_esp + n_tempo - 1},
     constr = TRUE, extraconstr = restricoes_t4, hyper = pc_prec)"
))

ajustar <- function(f, rotulo) {
  message("\nAjustando: ", rotulo, " ...")
  t0 <- proc.time()[3]
  m <- try(INLA::inla(f, family = "poisson", data = d, E = esperado,
                      control.compute = controle$compute,
                      control.predictor = controle$predictor,
                      control.inla = list(strategy = "adaptive"),
                      verbose = FALSE), silent = TRUE)
  if (inherits(m, "try-error")) {
    message("  falhou: ", substr(as.character(m), 1, 150))
    return(NULL)
  }
  message(glue::glue("  ok em {round(proc.time()[3]-t0,1)}s | ",
                     "DIC = {round(m$dic$dic,1)} | WAIC = {round(m$waic$waic,1)}"))
  m
}

modelos <- list(
  `Espacial (BYM2)`              = ajustar(f_espacial,   "Espacial BYM2"),
  `Espaço-temporal aditivo`      = ajustar(f_st_aditivo, "ST aditivo"),
  `Interação tipo I`             = ajustar(f_st_int1,    "ST interação I"),
  `Interação tipo IV`            = ajustar(f_st_int4,    "ST interação IV")
)
modelos <- purrr::compact(modelos)

tab_ajuste <- purrr::imap_dfr(modelos, ~ tibble::tibble(
  modelo = .y,
  DIC = .x$dic$dic, pD = .x$dic$p.eff,
  WAIC = .x$waic$waic, pW = .x$waic$p.eff,
  LCPO = -mean(log(.x$cpo$cpo), na.rm = TRUE),
  n_falhas_cpo = sum(.x$cpo$failure > 0, na.rm = TRUE)
)) |>
  dplyr::arrange(WAIC) |>
  dplyr::mutate(delta_WAIC = WAIC - min(WAIC))

salvar_tabela(tab_ajuste, "37_ajuste_modelos_bayesianos")
message("\nComparação de modelos bayesianos:")
print(as.data.frame(tab_ajuste))

melhor_nome <- tab_ajuste$modelo[1]
m_final <- modelos[[melhor_nome]]
message("\n>>> Modelo bayesiano selecionado: ", melhor_nome)


# =============================================================================
# 10.5 EFEITOS FIXOS — razões de risco com intervalos de credibilidade
# =============================================================================

fixos <- m_final$summary.fixed |>
  tibble::rownames_to_column("termo") |>
  janitor::clean_names() |>
  dplyr::transmute(
    termo,
    RR      = exp(mean),
    RR_li   = exp(x0_025quant),
    RR_ls   = exp(x0_975quant),
    RR_IC   = fmt_rr(RR, RR_li, RR_ls),
    # Probabilidade posterior de efeito de risco (RR > 1) — substitui o p-valor
    prob_RR_maior_1 = purrr::map_dbl(
      m_final$marginals.fixed,
      ~ 1 - INLA::inla.pmarginal(0, .x)
    )
  ) |>
  dplyr::filter(termo != "(Intercept)")

salvar_tabela(fixos, "38_efeitos_fixos_bayesianos")
message("\nEfeitos fixos (RR e IC95% de credibilidade):")
print(as.data.frame(fixos))


# =============================================================================
# 10.6 HIPERPARÂMETROS — repartição da variância
# =============================================================================
# phi ~ 1  => risco quase todo espacialmente estruturado (contágio/ambiente
#             compartilhado); phi ~ 0 => heterogeneidade idiossincrática
#             municipal (mais compatível com variação de vigilância).

hiper <- m_final$summary.hyperpar |>
  tibble::rownames_to_column("parametro") |>
  janitor::clean_names()
salvar_tabela(hiper, "39_hiperparametros_bayesianos")
print(as.data.frame(hiper))

lin_phi <- grep("Phi for id_esp", rownames(m_final$summary.hyperpar))
if (length(lin_phi)) {
  phi <- m_final$summary.hyperpar[lin_phi, "mean"]
  message(glue::glue(
    "\nphi = {round(phi,3)} -> {round(phi*100,1)}% da variância do efeito ",
    "aleatório espacial é estruturada (compartilhada com vizinhos)."
  ))
}


# =============================================================================
# 10.7 RISCO RELATIVO SUAVIZADO E PROBABILIDADE DE EXCESSO
# =============================================================================
# A probabilidade posterior P(RR > 1 | dados) é a saída mais útil para gestão:
# identifica municípios com risco elevado com incerteza explicitada.
# Convenção usual: P > 0,80 = "área provavelmente de risco elevado";
#                  P > 0,95 = "evidência forte".

# As marginais só existem se `return.marginals.predictor = TRUE` tiver sido
# aceito pela versão do INLA instalada; sem elas, a probabilidade de excesso
# é aproximada pela normal na escala log.
marg <- m_final$marginals.fitted.values
prob_exc <- if (length(marg) == nrow(d)) {
  purrr::map_dbl(marg, ~ 1 - INLA::inla.pmarginal(1, .x))
} else {
  message("Marginais da preditora indisponíveis — usando aproximação normal.")
  mu <- m_final$summary.fitted.values$mean
  sd_log <- (log(m_final$summary.fitted.values$`0.975quant`) -
             log(m_final$summary.fitted.values$`0.025quant`)) / (2 * 1.96)
  pnorm(log(mu) / sd_log)
}

rr_it <- tibble::tibble(
  code_muni = d$code_muni,
  ano       = d$ano,
  observado = d$casos_sinan,
  esperado  = d$esperado,
  RR        = m_final$summary.fitted.values$mean,
  RR_li     = m_final$summary.fitted.values$`0.025quant`,
  RR_ls     = m_final$summary.fitted.values$`0.975quant`,
  prob_RR_maior_1 = prob_exc
) |>
  dplyr::mutate(
    classificacao = dplyr::case_when(
      prob_RR_maior_1 > 0.95 ~ "Risco elevado — evidência forte",
      prob_RR_maior_1 > 0.80 ~ "Risco elevado — evidência moderada",
      prob_RR_maior_1 < 0.20 ~ "Risco reduzido",
      TRUE                   ~ "Indeterminado"
    )
  )

salvar_tabela(rr_it, "40_rr_bayesiano_municipio_ano")
saveRDS(rr_it, file.path(PARAMS$dir_processados, "rr_bayesiano.rds"))

message("\nClassificação de risco (município-ano):")
print(table(rr_it$classificacao))

# Municípios com risco elevado persistente ao longo de todo o período
persistentes <- rr_it |>
  dplyr::group_by(code_muni) |>
  dplyr::summarise(n_anos_risco = sum(prob_RR_maior_1 > 0.80),
                   RR_medio = mean(RR), .groups = "drop") |>
  dplyr::filter(n_anos_risco >= 6) |>
  dplyr::left_join(sf::st_drop_geometry(malha_est)[, c("code_muni", "name_muni",
                                                       "abbrev_state")],
                   by = "code_muni") |>
  dplyr::arrange(dplyr::desc(RR_medio))

salvar_tabela(persistentes, "41_municipios_risco_persistente")
message("\n", nrow(persistentes),
        " município(s) com P(RR>1) > 0,80 em pelo menos 6 dos 8 anos.")


# =============================================================================
# 10.8 SALVAR
# =============================================================================

saveRDS(list(modelos = modelos, melhor = melhor_nome, dados = d,
             grafo = arq_grafo),
        file.path(PARAMS$dir_modelos, "modelos_inla.rds"))

message("\nModelagem bayesiana concluída. Próximo: R/11_mapas_tabelas.R")

}  # fim do if (!knitr_exit)
