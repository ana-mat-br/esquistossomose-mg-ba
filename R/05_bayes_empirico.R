# =============================================================================
# 05_bayes_empirico.R — Estimadores Bayesianos Empíricos (Marshall, 1991)
#
# PROBLEMA: taxas brutas em municípios de população pequena são instáveis —
# um único caso em um município de 3.000 habitantes gera taxa de 33/100 mil,
# que aparece no mapa como "área de alto risco" sem qualquer significado
# epidemiológico. É o problema clássico do "pequeno número".
#
# SOLUÇÃO: encolhimento (shrinkage) da taxa bruta em direção a uma média,
# com peso proporcional à confiabilidade da estimativa local:
#   - EB GLOBAL: encolhe para a média de TODA a área de estudo
#   - EB LOCAL : encolhe para a média da VIZINHANÇA — preferível aqui, porque
#                preserva o gradiente regional do corredor endêmico em vez de
#                achatá-lo contra uma média única MG+BA
#
# Saída: dados/processados/agregado_eb.rds e painel_eb.rds
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

agregado  <- readRDS(file.path(PARAMS$dir_processados, "agregado.rds"))
painel    <- readRDS(file.path(PARAMS$dir_processados, "painel.rds"))
malha_est <- readRDS(file.path(PARAMS$dir_processados, "malha_est.rds"))

# Garantir alinhamento entre a malha e os dados — pré-requisito absoluto de
# TODA análise espacial. Um desalinhamento silencioso aqui invalida Moran,
# LISA, SAR e GWR sem gerar erro.
stopifnot(identical(malha_est$code_muni, sort(malha_est$code_muni)))
agregado <- agregado[match(malha_est$code_muni, agregado$code_muni), ]
stopifnot(identical(agregado$code_muni, malha_est$code_muni))


# =============================================================================
# 1. VIZINHANÇA
# =============================================================================

viz <- construir_vizinhanca(malha_est,
                            tipo  = PARAMS$viz_tipo,
                            k     = PARAMS$viz_k,
                            style = PARAMS$viz_style)

message(glue::glue(
  "Vizinhança '{PARAMS$viz_tipo}': média de {round(mean(spdep::card(viz$nb)),1)} ",
  "vizinhos por município (mín {min(spdep::card(viz$nb))}, ",
  "máx {max(spdep::card(viz$nb))})."
))

saveRDS(viz, file.path(PARAMS$dir_processados, "vizinhanca.rds"))

# Figura diagnóstica da rede de vizinhança (vai para material suplementar)
png(file.path(PARAMS$dir_figuras, "figS01_rede_vizinhanca.png"),
    width = 2000, height = 2200, res = 240)
plot(sf::st_geometry(malha_est), border = "grey80",
     main = glue::glue("Matriz de vizinhança — contiguidade {PARAMS$viz_tipo}"))
plot(viz$nb, viz$coords, add = TRUE, col = "#B2182B", lwd = .5, points = FALSE)
dev.off()


# =============================================================================
# 2. EB PARA O AGREGADO DO PERÍODO (2018-2025)
# =============================================================================

# --- Desfecho 1: incidência SINAN (denominador = pessoa-ano) ----------------
eb_g_sinan <- eb_global(agregado$casos_sinan, agregado$pessoa_ano,
                        base = PARAMS$base_taxa)
eb_l_sinan <- eb_local(agregado$casos_sinan, agregado$pessoa_ano,
                       nb = viz$nb, base = PARAMS$base_taxa)

# --- Desfecho 2: positividade PCE (denominador = examinados) ---------------
# Municípios sem exame no período (examinados = 0) não têm positividade
# definida. Passamos NA em vez de zero — zero seria afirmar "0% positivo",
# o que é falso: não houve medição.
tem_pce <- agregado$examinados > 0

eb_l_pce <- tibble::tibble(taxa_bruta = NA_real_, taxa_ebl = NA_real_)[
  rep(1, nrow(agregado)), ]

if (sum(tem_pce) > 10) {
  # EB local exige vizinhança consistente; restringimos o subconjunto e
  # reconstruímos a vizinhança APENAS entre municípios com dado de PCE.
  malha_pce <- malha_est[tem_pce, ]
  viz_pce   <- construir_vizinhanca(malha_pce, tipo = PARAMS$viz_tipo,
                                    k = PARAMS$viz_k, style = PARAMS$viz_style)
  eb_sub <- eb_local(agregado$positivos[tem_pce],
                     agregado$examinados[tem_pce],
                     nb = viz_pce$nb, base = 100)   # base 100 => percentual
  eb_l_pce$taxa_bruta[tem_pce] <- eb_sub$taxa_bruta
  eb_l_pce$taxa_ebl[tem_pce]   <- eb_sub$taxa_ebl
  saveRDS(list(malha = malha_pce, viz = viz_pce),
          file.path(PARAMS$dir_processados, "viz_pce.rds"))
} else {
  warning("Menos de 10 municípios com dado de PCE — EB local do PCE não estimado.")
}

agregado_eb <- agregado |>
  dplyr::mutate(
    inc_sinan_bruta = eb_g_sinan$taxa_bruta,
    inc_sinan_ebg   = eb_g_sinan$taxa_ebg,
    inc_sinan_ebl   = eb_l_sinan$taxa_ebl,
    posit_bruta     = eb_l_pce$taxa_bruta,
    posit_ebl       = eb_l_pce$taxa_ebl,
    # log das taxas EB — variável-resposta dos modelos gaussianos espaciais
    # (script 09). +0,5/base evita log(0) sem distorcer taxas não nulas.
    log_inc_ebl     = log(inc_sinan_ebl + 0.5),
    log_posit_ebl   = log(posit_ebl + 0.5)
  )


# =============================================================================
# 3. DIAGNÓSTICO DO ENCOLHIMENTO
# =============================================================================
# Verifica se o EB fez o que deveria: reduzir variabilidade onde a população é
# pequena e quase não alterar onde é grande. Se o encolhimento for uniforme,
# algo está errado.

diag_eb <- agregado_eb |>
  dplyr::mutate(
    quintil_pop = dplyr::ntile(pessoa_ano, 5),
    razao_encolhimento = dplyr::if_else(inc_sinan_bruta > 0,
                                        inc_sinan_ebl / inc_sinan_bruta, NA_real_)
  ) |>
  dplyr::group_by(quintil_pop) |>
  dplyr::summarise(
    n = dplyr::n(),
    pop_mediana      = median(pessoa_ano),
    cv_bruta         = sd(inc_sinan_bruta, na.rm = TRUE) /
                       mean(inc_sinan_bruta, na.rm = TRUE),
    cv_ebl           = sd(inc_sinan_ebl, na.rm = TRUE) /
                       mean(inc_sinan_ebl, na.rm = TRUE),
    reducao_cv_pct   = (1 - cv_ebl / cv_bruta) * 100,
    razao_mediana    = median(razao_encolhimento, na.rm = TRUE),
    .groups = "drop"
  )

salvar_tabela(diag_eb, "06_diagnostico_encolhimento_eb")
message("\nEncolhimento por quintil de população (1 = menor):")
print(as.data.frame(diag_eb))

if (with(diag_eb, reducao_cv_pct[quintil_pop == 1] < reducao_cv_pct[quintil_pop == 5])) {
  warning("O encolhimento foi MENOR nos municípios pequenos — comportamento ",
          "inesperado. Verifique os denominadores.")
}

# Correlação entre estimadores (esperada alta; divergências grandes indicam
# municípios cuja taxa bruta era artefato de pequeno número)
message(glue::glue(
  "\nCorrelação de Spearman bruta x EB local: ",
  "{round(cor(agregado_eb$inc_sinan_bruta, agregado_eb$inc_sinan_ebl, ",
  "method='spearman', use='complete.obs'), 3)}"
))


# =============================================================================
# 4. EB POR ANO (para os mapas da série)
# =============================================================================
# Aplicado ano a ano com a MESMA estrutura de vizinhança, permitindo
# comparabilidade entre os painéis do mapa multitemporal.

painel_eb <- painel |>
  dplyr::arrange(ano, code_muni) |>
  dplyr::group_by(ano) |>
  dplyr::group_modify(function(d, key) {
    d <- d[match(malha_est$code_muni, d$code_muni), ]
    e <- eb_local(d$casos_sinan, d$populacao, nb = viz$nb, base = PARAMS$base_taxa)
    d$inc_bruta <- e$taxa_bruta
    d$inc_ebl   <- e$taxa_ebl
    d
  }) |>
  dplyr::ungroup()

# Positividade anual só onde houve exame naquele ano
painel_eb <- painel_eb |>
  dplyr::group_by(ano) |>
  dplyr::group_modify(function(d, key) {
    ok <- !is.na(d$examinados) & d$examinados > 0
    d$posit_ebl <- NA_real_
    if (sum(ok) > 10) {
      mp  <- malha_est[malha_est$code_muni %in% d$code_muni[ok], ]
      vzp <- construir_vizinhanca(mp, tipo = PARAMS$viz_tipo, k = PARAMS$viz_k)
      dd  <- d[ok, ][match(mp$code_muni, d$code_muni[ok]), ]
      e   <- eb_local(dd$positivos, dd$examinados, nb = vzp$nb, base = 100)
      d$posit_ebl[match(mp$code_muni, d$code_muni)] <- e$taxa_ebl
    }
    d
  }) |>
  dplyr::ungroup()


# =============================================================================
# 5. SALVAR
# =============================================================================

saveRDS(agregado_eb, file.path(PARAMS$dir_processados, "agregado_eb.rds"))
saveRDS(painel_eb,   file.path(PARAMS$dir_processados, "painel_eb.rds"))

salvar_tabela(
  agregado_eb |>
    dplyr::select(code_muni, name_muni, uf = abbrev_state, regiao,
                  casos_sinan, pessoa_ano,
                  inc_sinan_bruta, inc_sinan_ebg, inc_sinan_ebl,
                  examinados, positivos, posit_bruta, posit_ebl,
                  rme_sinan, rme_pce) |>
    dplyr::arrange(dplyr::desc(inc_sinan_ebl)),
  "07_taxas_municipais_eb"
)


# =============================================================================
# 6. FIGURA — efeito do encolhimento
# =============================================================================

g_eb <- agregado_eb |>
  ggplot2::ggplot(ggplot2::aes(inc_sinan_bruta, inc_sinan_ebl)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_point(ggplot2::aes(size = pessoa_ano, colour = regiao), alpha = .6) +
  ggplot2::scale_size_continuous(
    name = "População\n(pessoa-ano)", labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  ggplot2::scale_colour_manual(values = c("Norte de Minas" = "#1B7837",
                                          "Bahia" = "#762A83"), name = NULL) +
  ggplot2::scale_x_continuous(trans = "log1p") +
  ggplot2::scale_y_continuous(trans = "log1p") +
  ggplot2::labs(
    title = "Efeito do estimador bayesiano empírico local",
    subtitle = "Pontos abaixo da diagonal: taxa bruta superestimada por pequeno número",
    x = "Incidência bruta (/100 mil pessoa-ano, escala log)",
    y = "Incidência bayesiana empírica local (escala log)",
    caption = "Marshall (1991). O encolhimento é maior quanto menor a população."
  ) +
  theme_grafico()

salvar_figura(g_eb, "fig02_efeito_bayes_empirico", largura = 8, altura = 6)

message("\nEB concluído. Próximo: R/06_moran_lisa.R")
