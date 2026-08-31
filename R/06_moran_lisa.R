# =============================================================================
# 06_moran_lisa.R — Autocorrelação espacial global e local
#
#   6.1 I de Moran global (bruta, EB global, EB local)
#   6.2 Índice de Moran bayesiano empírico (Assunção & Reis, 1999) — corrige o
#       viés do I de Moran quando as populações municipais são muito desiguais,
#       situação exata do corredor MG/BA (Salvador vs. municípios de 3 mil hab.)
#   6.3 Correlograma espacial (I por ordem de vizinhança)
#   6.4 LISA / Moran local com permutação condicional + correção FDR
#   6.5 Gi* de Getis-Ord (hot/cold spots)
#   6.6 LISA ano a ano (persistência espaço-temporal dos clusters)
#   6.7 Moran bivariado: desfecho x indicadores socioeconômicos
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

agregado_eb <- readRDS(file.path(PARAMS$dir_processados, "agregado_eb.rds"))
painel_eb   <- readRDS(file.path(PARAMS$dir_processados, "painel_eb.rds"))
malha_est   <- readRDS(file.path(PARAMS$dir_processados, "malha_est.rds"))
viz         <- readRDS(file.path(PARAMS$dir_processados, "vizinhanca.rds"))

stopifnot(identical(agregado_eb$code_muni, malha_est$code_muni))
listw <- viz$listw
nb    <- viz$nb

set.seed(PARAMS$semente)


# =============================================================================
# 6.1 I DE MORAN GLOBAL
# =============================================================================
# Reportamos o I para os três estimadores. A comparação é informativa: se o I
# da taxa bruta for muito menor que o da taxa EB, o ruído de pequeno número
# estava mascarando a estrutura espacial real.

vars_moran <- c(
  "Incidência bruta"                = "inc_sinan_bruta",
  "Incidência EB global"            = "inc_sinan_ebg",
  "Incidência EB local"             = "inc_sinan_ebl",
  "RME (razão de morbidade padr.)"  = "rme_sinan",
  "Positividade bruta"              = "posit_bruta",
  "Positividade EB local"           = "posit_ebl"
)

moran_tab <- purrr::imap_dfr(vars_moran, function(v, rot) {
  if (!v %in% names(agregado_eb)) return(NULL)
  x <- agregado_eb[[v]]
  if (sum(is.finite(x)) < 20) return(NULL)
  # Municípios com NA (ex.: sem PCE) são excluídos com sua vizinhança
  ok <- is.finite(x)
  lw <- if (all(ok)) listw else
    spdep::nb2listw(spdep::subset.nb(nb, ok), style = PARAMS$viz_style,
                    zero.policy = TRUE)
  res <- moran_global(x[ok], lw, nsim = PARAMS$n_sim)
  dplyr::mutate(res, indicador = rot, variavel = v, n = sum(ok), .before = 1)
})


# =============================================================================
# 6.2 ÍNDICE DE MORAN BAYESIANO EMPÍRICO (Assunção & Reis, 1999)
# =============================================================================
# O I de Moran padrão aplicado a taxas assume variância homogênea. Com
# populações heterogêneas, a variância da taxa é maior em municípios pequenos,
# inflando ou deflacionando o I. O EBI corrige isso na própria estatística.
# É o teste global que deve ser reportado como principal em dados de contagem.

ebi_sinan <- spdep::EBImoran.mc(
  n = agregado_eb$casos_sinan,
  x = agregado_eb$pessoa_ano,
  listw = listw, nsim = PARAMS$n_sim, zero.policy = TRUE
)

moran_tab <- dplyr::bind_rows(
  moran_tab,
  tibble::tibble(
    indicador = "Incidência SINAN — I de Moran EB (Assunção-Reis)",
    variavel  = "EBI_sinan",
    n         = nrow(agregado_eb),
    I         = unname(ebi_sinan$statistic),
    I_esperado = -1 / (nrow(agregado_eb) - 1),
    p_valor   = ebi_sinan$p.value,
    n_sim     = PARAMS$n_sim
  )
)

tem_pce <- agregado_eb$examinados > 0
if (sum(tem_pce) > 20) {
  vp <- readRDS(file.path(PARAMS$dir_processados, "viz_pce.rds"))
  ebi_pce <- spdep::EBImoran.mc(
    n = agregado_eb$positivos[tem_pce],
    x = agregado_eb$examinados[tem_pce],
    listw = vp$viz$listw, nsim = PARAMS$n_sim, zero.policy = TRUE
  )
  moran_tab <- dplyr::bind_rows(moran_tab, tibble::tibble(
    indicador = "Positividade PCE — I de Moran EB (Assunção-Reis)",
    variavel = "EBI_pce", n = sum(tem_pce),
    I = unname(ebi_pce$statistic), I_esperado = -1 / (sum(tem_pce) - 1),
    p_valor = ebi_pce$p.value, n_sim = PARAMS$n_sim
  ))
}

moran_tab <- moran_tab |>
  dplyr::mutate(
    p_fmt = fmt_p(p_valor),
    interpretacao = dplyr::case_when(
      p_valor > PARAMS$alfa ~ "Distribuição aleatória",
      I > 0                 ~ "Autocorrelação positiva (clusterização)",
      TRUE                  ~ "Autocorrelação negativa (dispersão)"
    )
  )

salvar_tabela(moran_tab, "08_moran_global")
print(as.data.frame(dplyr::select(moran_tab, indicador, n, I, p_fmt, interpretacao)))


# =============================================================================
# 6.3 CORRELOGRAMA ESPACIAL
# =============================================================================
# I de Moran por ordem de vizinhança (1ª, 2ª, ... ordem). Mostra a ESCALA da
# dependência espacial: até quantas ordens de vizinhos o risco permanece
# correlacionado. Informa a escolha do raio máximo no SaTScan.

correlograma <- spdep::sp.correlogram(
  nb, agregado_eb$inc_sinan_ebl, order = 6, method = "I",
  style = PARAMS$viz_style, zero.policy = TRUE, randomisation = TRUE
)

corr_df <- as.data.frame(print(correlograma, p.adj.method = "holm")) |>
  tibble::rownames_to_column("ordem") |>
  janitor::clean_names()
salvar_tabela(corr_df, "09_correlograma_espacial")

png(file.path(PARAMS$dir_figuras, "figS02_correlograma.png"),
    width = 1800, height = 1400, res = 240)
plot(correlograma, main = "Correlograma espacial — incidência EB local",
     xlab = "Ordem de vizinhança", ylab = "I de Moran")
dev.off()


# =============================================================================
# 6.4 LISA — MORAN LOCAL
# =============================================================================

lisa_sinan <- lisa(agregado_eb$inc_sinan_ebl, listw,
                   nsim = PARAMS$n_sim, alfa = PARAMS$alfa,
                   ajuste_p = PARAMS$ajuste_p_lisa)

# Covariáveis contextuais efetivamente disponíveis
vars_ctx <- intersect(c("idhm", "ivs", "pct_agua_rede", "pct_esgoto_adeq",
                        "pct_agua_superficial", "pct_esgoto_em_corpo_dagua",
                        "pct_lixo_coletado", "pct_rural"), names(agregado_eb))

resultado_lisa <- agregado_eb |>
  dplyr::select(code_muni, name_muni, uf = abbrev_state, regiao,
                casos_sinan, pessoa_ano, inc_sinan_bruta, inc_sinan_ebl,
                rme_sinan, dplyr::all_of(vars_ctx)) |>
  dplyr::bind_cols(lisa_sinan)

# LISA para positividade (subconjunto com PCE)
if (sum(tem_pce) > 20) {
  vp <- readRDS(file.path(PARAMS$dir_processados, "viz_pce.rds"))
  lisa_pce <- lisa(agregado_eb$posit_ebl[tem_pce], vp$viz$listw,
                   nsim = PARAMS$n_sim, alfa = PARAMS$alfa,
                   ajuste_p = PARAMS$ajuste_p_lisa)
  resultado_lisa$cluster_pce <- factor(NA, levels = names(PAL_LISA))
  resultado_lisa$p_lisa_pce  <- NA_real_
  resultado_lisa$cluster_pce[tem_pce][match(vp$malha$code_muni,
                                            agregado_eb$code_muni[tem_pce])] <-
    lisa_pce$cluster
  resultado_lisa$p_lisa_pce[tem_pce][match(vp$malha$code_muni,
                                           agregado_eb$code_muni[tem_pce])] <-
    lisa_pce$p_lisa_aj
}

salvar_tabela(resultado_lisa, "10_lisa_municipal")
saveRDS(resultado_lisa, file.path(PARAMS$dir_processados, "lisa.rds"))

message("\nClusters LISA — incidência SINAN (EB local, p ajustado por ",
        PARAMS$ajuste_p_lisa, "):")
print(table(resultado_lisa$cluster, resultado_lisa$regiao))

# Perfil socioeconômico dos clusters — resultado central do artigo
perfil_cluster <- resultado_lisa |>
  dplyr::filter(cluster != "Não significativo") |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    n = dplyr::n(),
    inc_mediana = median(inc_sinan_ebl, na.rm = TRUE),
    dplyr::across(dplyr::all_of(vars_ctx), ~ median(.x, na.rm = TRUE),
                  .names = "{.col}_mediana"),
    .groups = "drop"
  )
salvar_tabela(perfil_cluster, "11_perfil_socioeconomico_clusters")
print(as.data.frame(perfil_cluster))


# =============================================================================
# 6.5 Gi* DE GETIS-ORD
# =============================================================================
# Complementar ao LISA: o LISA identifica associação local (inclusive
# outliers Alto-Baixo); o Gi* identifica concentrações de valores altos ou
# baixos. Reportar ambos fortalece a evidência de cluster.

gstar <- getis_ord(agregado_eb$inc_sinan_ebl, nb,
                   alfa = PARAMS$alfa, ajuste_p = PARAMS$ajuste_p_lisa)
resultado_lisa <- dplyr::bind_cols(resultado_lisa, gstar)

concordancia <- table(
  LISA = resultado_lisa$cluster,
  `Gi*` = resultado_lisa$Gi_cluster
)
message("\nConcordância LISA x Gi*:")
print(concordancia)
salvar_tabela(as.data.frame(concordancia), "12_concordancia_lisa_gistar")
saveRDS(resultado_lisa, file.path(PARAMS$dir_processados, "lisa.rds"))


# =============================================================================
# 6.6 LISA ANO A ANO — persistência dos clusters
# =============================================================================
# Um cluster que aparece em 1 de 8 anos é ruído; um que persiste em 6+ anos é
# área prioritária de intervenção. Esta é a leitura espaço-temporal exigida
# pelo objetivo do estudo.

lisa_anual <- painel_eb |>
  dplyr::arrange(ano, code_muni) |>
  dplyr::group_by(ano) |>
  dplyr::group_modify(function(d, key) {
    d <- d[match(malha_est$code_muni, d$code_muni), ]
    if (sum(d$casos_sinan, na.rm = TRUE) < 20) {
      d$cluster_ano <- factor("Não significativo", levels = names(PAL_LISA))
      d$I_ano <- NA_real_
      return(d)
    }
    l <- lisa(d$inc_ebl, listw, nsim = 999, alfa = PARAMS$alfa,
              ajuste_p = PARAMS$ajuste_p_lisa)
    d$cluster_ano <- l$cluster
    d$I_ano <- l$Ii
    d
  }) |>
  dplyr::ungroup()

persistencia <- lisa_anual |>
  dplyr::group_by(code_muni) |>
  dplyr::summarise(
    n_anos_AA = sum(cluster_ano == "Alto-Alto", na.rm = TRUE),
    n_anos_BB = sum(cluster_ano == "Baixo-Baixo", na.rm = TRUE),
    n_anos_sig = sum(cluster_ano != "Não significativo", na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    classe_persistencia = dplyr::case_when(
      n_anos_AA >= 6 ~ "Cluster de alto risco persistente (≥6 anos)",
      n_anos_AA >= 3 ~ "Cluster de alto risco intermitente (3-5 anos)",
      n_anos_AA >= 1 ~ "Cluster de alto risco esporádico (1-2 anos)",
      n_anos_BB >= 6 ~ "Cluster de baixo risco persistente",
      TRUE           ~ "Sem cluster relevante"
    )
  ) |>
  dplyr::left_join(
    dplyr::select(agregado_eb, code_muni, name_muni, abbrev_state, regiao,
                  inc_sinan_ebl, dplyr::all_of(vars_ctx)),
    by = "code_muni"
  )

salvar_tabela(persistencia, "13_persistencia_clusters")
saveRDS(persistencia, file.path(PARAMS$dir_processados, "persistencia.rds"))

message("\nPersistência dos clusters de alto risco:")
print(table(persistencia$classe_persistencia))

# Série do I de Moran global ano a ano
moran_anual <- painel_eb |>
  dplyr::arrange(ano, code_muni) |>
  dplyr::group_by(ano) |>
  dplyr::group_modify(function(d, key) {
    d <- d[match(malha_est$code_muni, d$code_muni), ]
    moran_global(d$inc_ebl, listw, nsim = 999)
  }) |>
  dplyr::ungroup()
salvar_tabela(moran_anual, "14_moran_global_anual")


# =============================================================================
# 6.7 MORAN BIVARIADO — desfecho x determinantes
# =============================================================================
# Mede se a incidência de um município se associa ao valor do indicador
# socioeconômico NOS VIZINHOS — evidência de efeito de contexto/transbordamento
# (spillover), que a regressão não espacial não capta.
#
# ATENÇÃO: o Moran bivariado NÃO controla a autocorrelação de cada variável
# isoladamente; é exploratório. A inferência causal fica com os modelos 08-10.

moran_bivariado <- function(x, y, listw, nsim = 999) {
  # I_B = correlação entre x padronizado e a defasagem espacial de y padronizado
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20) return(NULL)
  lw <- if (all(ok)) listw else
    spdep::nb2listw(spdep::subset.nb(nb, ok), style = PARAMS$viz_style,
                    zero.policy = TRUE)
  xs <- as.numeric(scale(x[ok]))
  ys <- as.numeric(scale(y[ok]))
  lag_y <- spdep::lag.listw(lw, ys, zero.policy = TRUE)
  I_obs <- sum(xs * lag_y) / sum(xs^2)
  # Inferência por permutação de x (mantém a estrutura espacial de y)
  I_sim <- replicate(nsim, {
    xp <- sample(xs)
    sum(xp * lag_y) / sum(xp^2)
  })
  p <- (sum(abs(I_sim) >= abs(I_obs)) + 1) / (nsim + 1)
  tibble::tibble(I_bivariado = I_obs, p_valor = p, n = sum(ok))
}

vars_socio <- intersect(
  c("idhm", "idhm_renda", "idhm_educacao", "ivs", "gini",
    "pct_agua_rede", "pct_esgoto_adeq", "pct_lixo_coletado", "pct_rural",
    "cob_agua_snis", "cob_esgoto_snis"),
  names(agregado_eb)
)

tab_bivariado <- purrr::map_dfr(vars_socio, function(v) {
  r <- moran_bivariado(agregado_eb$inc_sinan_ebl, agregado_eb[[v]],
                       listw, nsim = 999)
  if (is.null(r)) return(NULL)
  dplyr::mutate(r, covariavel = v,
                rotulo = COVARIAVEIS$rotulo[match(v, COVARIAVEIS$var)],
                .before = 1)
}) |>
  dplyr::mutate(p_aj = p.adjust(p_valor, "fdr"),
                p_fmt = fmt_p(p_aj)) |>
  dplyr::arrange(dplyr::desc(abs(I_bivariado)))

salvar_tabela(tab_bivariado, "15_moran_bivariado")
message("\nMoran bivariado (incidência EB x determinantes nos vizinhos):")
print(as.data.frame(tab_bivariado))


# =============================================================================
# 6.8 DIAGRAMA DE ESPALHAMENTO DE MORAN
# =============================================================================

z  <- as.numeric(scale(agregado_eb$inc_sinan_ebl))
wz <- spdep::lag.listw(listw, z, zero.policy = TRUE)

g_moran <- tibble::tibble(z = z, wz = wz,
                          cluster = resultado_lisa$cluster,
                          nome = agregado_eb$name_muni) |>
  ggplot2::ggplot(ggplot2::aes(z, wz)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
  ggplot2::geom_point(ggplot2::aes(colour = cluster), size = 2, alpha = .8) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "black",
                       linewidth = .7, formula = y ~ x) +
  ggplot2::scale_colour_manual(values = PAL_LISA, name = "Cluster LISA") +
  ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
                    label = sprintf("I de Moran = %.3f (p %s)",
                                    moran_tab$I[moran_tab$variavel == "inc_sinan_ebl"],
                                    fmt_p(moran_tab$p_valor[moran_tab$variavel == "inc_sinan_ebl"]))) +
  ggplot2::labs(
    title = "Diagrama de espalhamento de Moran",
    subtitle = "Incidência de esquistossomose (EB local), 2018–2025",
    x = "Incidência padronizada (z)",
    y = "Defasagem espacial da incidência (Wz)"
  ) +
  theme_grafico()

salvar_figura(g_moran, "fig03_moran_scatter", largura = 8, altura = 6.5)

message("\nAnálise de autocorrelação concluída. Próximo: R/07_satscan.R")
