# =============================================================================
# 11_mapas_tabelas.R — Cartografia e tabelas finais
#
# Convenções cartográficas adotadas (declare-as na legenda das figuras):
#   - Quebras por quantis para taxas (comparabilidade entre painéis temporais)
#     e por valores fixos para RR/LISA (interpretação absoluta)
#   - Escala sequencial multi-hue (viridis/YlOrRd) para intensidade;
#     divergente (RdBu invertido) para razões em torno de 1
#   - Norte, escala gráfica e CRS declarados
#   - Municípios sem dado hachurados/cinza, nunca em branco (branco confunde
#     com "valor zero")
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

library(ggspatial)
library(patchwork)
library(classInt)

agregado_eb <- readRDS(file.path(PARAMS$dir_processados, "agregado_eb.rds"))
painel_eb   <- readRDS(file.path(PARAMS$dir_processados, "painel_eb.rds"))
malha_est   <- readRDS(file.path(PARAMS$dir_processados, "malha_est.rds"))
lisa_res    <- readRDS(file.path(PARAMS$dir_processados, "lisa.rds"))
persist     <- readRDS(file.path(PARAMS$dir_processados, "persistencia.rds"))

# Contornos estaduais para referência
uf_contorno <- malha_est |>
  dplyr::group_by(abbrev_state) |>
  dplyr::summarise(.groups = "drop")

CAPTION_FONTE <- paste(
  "Fontes: SINAN e PCE/SISPCE (MS); população: IBGE; malha: IBGE 2022.",
  "Projeção: SIRGAS 2000 (EPSG:4674)."
)

#' Quebras por quantis, com rótulos legíveis
quebras_quantil <- function(x, n = 5, dig = 1) {
  x <- x[is.finite(x)]
  br <- unique(quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE))
  if (length(br) < 3) br <- unique(pretty(x, n))
  br
}

mapa_base <- function(sf_dados, var, titulo, subtitulo, rotulo_legenda,
                      paleta = "YlOrRd", n_classes = 5, direcao = 1) {
  x <- sf_dados[[var]]
  br <- quebras_quantil(x, n_classes)
  sf_dados$.classe <- cut(x, breaks = br, include.lowest = TRUE, dig.lab = 4)
  ggplot2::ggplot(sf_dados) +
    ggplot2::geom_sf(ggplot2::aes(fill = .classe), colour = "white", linewidth = .05) +
    ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey20", linewidth = .4) +
    ggplot2::scale_fill_brewer(palette = paleta, direction = direcao,
                               na.value = "grey85", name = rotulo_legenda,
                               drop = FALSE) +
    ggspatial::annotation_scale(location = "bl", width_hint = .25,
                                text_cex = .6) +
    ggspatial::annotation_north_arrow(location = "tr",
                                      style = ggspatial::north_arrow_minimal(),
                                      height = grid::unit(.8, "cm"),
                                      width = grid::unit(.8, "cm")) +
    ggplot2::labs(title = titulo, subtitle = subtitulo, caption = CAPTION_FONTE) +
    theme_mapa()
}


# =============================================================================
# FIG 5 — Taxa bruta vs. bayesiana empírica local
# =============================================================================

sf_ag <- malha_est |>
  dplyr::left_join(sf::st_drop_geometry(
    dplyr::select(agregado_eb, code_muni, inc_sinan_bruta, inc_sinan_ebl,
                  posit_bruta, posit_ebl, rme_sinan)) |>
    tibble::as_tibble(), by = "code_muni")

# (o left_join acima aceita agregado_eb como data.frame simples)
if (!"inc_sinan_ebl" %in% names(sf_ag)) {
  sf_ag <- malha_est |>
    dplyr::left_join(dplyr::select(as.data.frame(agregado_eb),
                                   code_muni, inc_sinan_bruta, inc_sinan_ebl,
                                   posit_bruta, posit_ebl, rme_sinan),
                     by = "code_muni")
}

m_bruta <- mapa_base(sf_ag, "inc_sinan_bruta",
                     "A. Taxa bruta", "Incidência por 100 mil pessoa-ano",
                     "Casos/100 mil")
m_ebl <- mapa_base(sf_ag, "inc_sinan_ebl",
                   "B. Bayesiana empírica local", "Taxa suavizada (Marshall, 1991)",
                   "Casos/100 mil")

fig05 <- (m_bruta | m_ebl) +
  patchwork::plot_annotation(
    title = "Esquistossomose no corredor Norte de Minas–Bahia, 2018–2025",
    subtitle = "Efeito da suavização bayesiana sobre a instabilidade de pequenos números",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  )
salvar_figura(fig05, "fig05_taxa_bruta_vs_eb", largura = 13, altura = 7)


# =============================================================================
# FIG 6 — Clusters LISA
# =============================================================================

sf_lisa <- malha_est |>
  dplyr::left_join(dplyr::select(as.data.frame(lisa_res),
                                 code_muni, cluster, p_lisa_aj,
                                 dplyr::any_of(c("cluster_pce", "Gi_cluster"))),
                   by = "code_muni")

fig06 <- ggplot2::ggplot(sf_lisa) +
  ggplot2::geom_sf(ggplot2::aes(fill = cluster), colour = "white", linewidth = .05) +
  ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey20", linewidth = .4) +
  ggplot2::scale_fill_manual(values = PAL_LISA, name = "Cluster LISA",
                             na.value = "grey85", drop = FALSE) +
  ggspatial::annotation_scale(location = "bl", width_hint = .25, text_cex = .6) +
  ggspatial::annotation_north_arrow(location = "tr",
                                    style = ggspatial::north_arrow_minimal(),
                                    height = grid::unit(.8, "cm"),
                                    width = grid::unit(.8, "cm")) +
  ggplot2::labs(
    title = "Clusters espaciais de esquistossomose (LISA)",
    subtitle = glue::glue(
      "Incidência bayesiana empírica local, 2018–2025 | ",
      "{PARAMS$n_sim} permutações, p ajustado por {toupper(PARAMS$ajuste_p_lisa)}"),
    caption = CAPTION_FONTE
  ) +
  theme_mapa()

salvar_figura(fig06, "fig06_lisa", largura = 8, altura = 9)

# Gi* de Getis-Ord
if ("Gi_cluster" %in% names(sf_lisa)) {
  fig06b <- ggplot2::ggplot(sf_lisa) +
    ggplot2::geom_sf(ggplot2::aes(fill = Gi_cluster), colour = "white", linewidth = .05) +
    ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey20", linewidth = .4) +
    ggplot2::scale_fill_manual(
      values = c("Hot spot" = "#B2182B", "Cold spot" = "#2166AC",
                 "Não significativo" = "grey88"),
      name = "Gi* de Getis-Ord", na.value = "grey85", drop = FALSE) +
    ggplot2::labs(title = "Concentrações de risco (Gi* de Getis-Ord)",
                  subtitle = "Análise complementar ao LISA",
                  caption = CAPTION_FONTE) +
    theme_mapa()
  salvar_figura(fig06b, "fig07_getis_ord", largura = 8, altura = 9)
}


# =============================================================================
# FIG 8 — Clusters de varredura (Kulldorff / SaTScan)
# =============================================================================

arq_clusters <- file.path(PARAMS$dir_processados, "clusters.gpkg")
if (file.exists(arq_clusters)) {
  clusters_sf <- sf::st_read(arq_clusters, quiet = TRUE) |>
    dplyr::filter(desfecho == "Incidência SINAN")

  fig08 <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf_ag, ggplot2::aes(fill = inc_sinan_ebl),
                     colour = "white", linewidth = .05) +
    ggplot2::scale_fill_viridis_c(option = "rocket", direction = -1,
                                  trans = "log1p", na.value = "grey85",
                                  name = "Incidência EB\n(/100 mil)") +
    ggplot2::geom_sf(data = clusters_sf, fill = NA, colour = "#00429d",
                     linewidth = 1) +
    ggplot2::geom_sf_text(data = clusters_sf,
                          ggplot2::aes(label = paste0("C", cluster)),
                          colour = "#00429d", fontface = "bold", size = 3.5) +
    ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey20",
                     linewidth = .4) +
    ggspatial::annotation_scale(location = "bl", width_hint = .25, text_cex = .6) +
    ggplot2::labs(
      title = "Clusters espaciais de alto risco (varredura de Kulldorff)",
      subtitle = glue::glue("Modelo Poisson, janela máxima = ",
                            "{PARAMS$satscan_ubpop*100}% da população em risco"),
      caption = CAPTION_FONTE) +
    theme_mapa()

  salvar_figura(fig08, "fig08_clusters_varredura", largura = 8, altura = 9)
}


# =============================================================================
# FIG 9 — Painel multitemporal
# =============================================================================

sf_anual <- painel_eb |>
  dplyr::select(code_muni, ano, inc_ebl) |>
  dplyr::left_join(dplyr::select(malha_est, code_muni), by = "code_muni") |>
  sf::st_as_sf()

# Quebras COMUNS a todos os anos — essencial para comparabilidade visual
br_comum <- quebras_quantil(sf_anual$inc_ebl, 5)
sf_anual$.classe <- cut(sf_anual$inc_ebl, breaks = br_comum,
                        include.lowest = TRUE, dig.lab = 4)

fig09 <- ggplot2::ggplot(sf_anual) +
  ggplot2::geom_sf(ggplot2::aes(fill = .classe), colour = NA) +
  ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey30", linewidth = .25) +
  ggplot2::facet_wrap(~ ano, nrow = 2) +
  ggplot2::scale_fill_brewer(palette = "YlOrRd", na.value = "grey85",
                             name = "Incidência EB\n(/100 mil hab.)", drop = FALSE) +
  ggplot2::labs(
    title = "Distribuição espaço-temporal da esquistossomose, 2018–2025",
    subtitle = "Taxas bayesianas empíricas locais; classes idênticas em todos os painéis",
    caption = CAPTION_FONTE) +
  theme_mapa(base_size = 10)

salvar_figura(fig09, "fig09_painel_multitemporal", largura = 14, altura = 8)


# =============================================================================
# FIG 10 — Persistência dos clusters
# =============================================================================

sf_persist <- malha_est |>
  dplyr::left_join(dplyr::select(as.data.frame(persist), code_muni,
                                 n_anos_AA, classe_persistencia),
                   by = "code_muni") |>
  dplyr::mutate(classe_persistencia = factor(
    classe_persistencia,
    levels = c("Cluster de alto risco persistente (≥6 anos)",
               "Cluster de alto risco intermitente (3-5 anos)",
               "Cluster de alto risco esporádico (1-2 anos)",
               "Cluster de baixo risco persistente",
               "Sem cluster relevante")))

fig10 <- ggplot2::ggplot(sf_persist) +
  ggplot2::geom_sf(ggplot2::aes(fill = classe_persistencia),
                   colour = "white", linewidth = .05) +
  ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey20", linewidth = .4) +
  ggplot2::scale_fill_manual(
    values = c("#67001F", "#D6604D", "#F4A582", "#4393C3", "grey88"),
    name = "Persistência do cluster", na.value = "grey85", drop = FALSE) +
  ggplot2::labs(
    title = "Persistência temporal dos clusters de alto risco",
    subtitle = "Número de anos (2018–2025) classificados como Alto-Alto pelo LISA",
    caption = paste(CAPTION_FONTE,
                    "\nÁreas em vermelho escuro são prioritárias para intervenção.")) +
  theme_mapa()

salvar_figura(fig10, "fig10_persistencia_clusters", largura = 8, altura = 9)


# =============================================================================
# FIG 11 — Superfície de coeficientes da GWR
# =============================================================================

arq_gwr <- file.path(PARAMS$dir_processados, "gwr_superficie.gpkg")
if (file.exists(arq_gwr)) {
  gwr_sf <- sf::st_read(arq_gwr, quiet = TRUE)

  vars_gwr <- intersect(c("pct_agua_rede", "pct_esgoto_adeq", "idhm", "pct_rural"),
                        names(gwr_sf))

  mapas_gwr <- purrr::map(vars_gwr, function(v) {
    rot <- COVARIAVEIS$rotulo[match(v, COVARIAVEIS$var)]
    lim <- max(abs(gwr_sf[[v]]), na.rm = TRUE)
    tv  <- paste0(v, "_TV")
    g <- ggplot2::ggplot(gwr_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data[[v]]), colour = NA) +
      ggplot2::scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim),
                                    name = "Coeficiente\nlocal") +
      ggplot2::labs(title = rot) +
      theme_mapa(base_size = 9)
    # Hachura nos municípios onde o coeficiente não é localmente significativo
    if (tv %in% names(gwr_sf)) {
      nao_sig <- dplyr::filter(gwr_sf, abs(.data[[tv]]) <= 1.96)
      g <- g + ggplot2::geom_sf(data = nao_sig, fill = "grey70", alpha = .55,
                                colour = NA)
    }
    g
  })

  fig11 <- patchwork::wrap_plots(mapas_gwr, ncol = 2) +
    patchwork::plot_annotation(
      title = "Variação espacial dos efeitos (GWR)",
      subtitle = paste("Coeficientes locais sobre log da taxa bayesiana empírica.",
                       "Áreas acinzentadas: |t| ≤ 1,96 (não significativo localmente)."),
      caption = CAPTION_FONTE,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))

  salvar_figura(fig11, "fig11_gwr_coeficientes", largura = 11, altura = 10)

  # R² local — onde o modelo explica bem e onde falha
  if ("Local_R2" %in% names(gwr_sf)) {
    fig11b <- ggplot2::ggplot(gwr_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = Local_R2), colour = NA) +
      ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey20",
                       linewidth = .4) +
      ggplot2::scale_fill_viridis_c(option = "mako", direction = -1,
                                    name = expression(R^2~"local")) +
      ggplot2::labs(title = "Poder explicativo local do modelo (GWR)",
                    subtitle = "Baixo R² local indica determinantes não incluídos no modelo",
                    caption = CAPTION_FONTE) +
      theme_mapa()
    salvar_figura(fig11b, "fig12_gwr_r2_local", largura = 8, altura = 9)
  }
}


# =============================================================================
# FIG 13 — Risco relativo bayesiano e probabilidade de excesso
# =============================================================================

arq_rr <- file.path(PARAMS$dir_processados, "rr_bayesiano.rds")
if (file.exists(arq_rr)) {
  rr_bayes <- readRDS(arq_rr)

  rr_medio <- rr_bayes |>
    dplyr::group_by(code_muni) |>
    dplyr::summarise(RR = mean(RR),
                     prob = mean(prob_RR_maior_1), .groups = "drop")

  sf_rr <- malha_est |> dplyr::left_join(rr_medio, by = "code_muni")

  m_rr <- ggplot2::ggplot(sf_rr) +
    ggplot2::geom_sf(ggplot2::aes(fill = RR), colour = NA) +
    ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey20", linewidth = .4) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                                  midpoint = 1, trans = "log10",
                                  name = "RR posterior", na.value = "grey85") +
    ggplot2::labs(title = "A. Risco relativo suavizado (BYM2)") +
    theme_mapa()

  m_prob <- ggplot2::ggplot(sf_rr) +
    ggplot2::geom_sf(ggplot2::aes(fill = cut(prob,
                        breaks = c(0, .2, .8, .95, 1),
                        labels = c("<0,20 (risco reduzido)", "0,20–0,80",
                                   "0,80–0,95", ">0,95 (risco elevado)"),
                        include.lowest = TRUE)), colour = NA) +
    ggplot2::geom_sf(data = uf_contorno, fill = NA, colour = "grey20", linewidth = .4) +
    ggplot2::scale_fill_brewer(palette = "RdYlBu", direction = -1,
                               name = "P(RR > 1 | dados)", na.value = "grey85") +
    ggplot2::labs(title = "B. Probabilidade posterior de risco excedente") +
    theme_mapa()

  fig13 <- (m_rr | m_prob) +
    patchwork::plot_annotation(
      title = "Modelo bayesiano hierárquico espaço-temporal (BYM2)",
      subtitle = "Média do período 2018–2025, ajustada por indicadores socioeconômicos e de saneamento",
      caption = CAPTION_FONTE,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))

  salvar_figura(fig13, "fig13_bayesiano_rr_prob", largura = 13, altura = 7)
}


# =============================================================================
# FIG 14 — Determinantes: mapas dos indicadores socioeconômicos
# =============================================================================

vars_mapa <- intersect(c("idhm", "pct_agua_rede", "pct_esgoto_adeq", "pct_rural"),
                       names(agregado_eb))

if (length(vars_mapa) >= 2) {
  sf_det <- malha_est |>
    dplyr::left_join(dplyr::select(as.data.frame(agregado_eb),
                                   code_muni, dplyr::all_of(vars_mapa)),
                     by = "code_muni")

  mapas_det <- purrr::map(vars_mapa, function(v) {
    rot <- COVARIAVEIS$rotulo[match(v, COVARIAVEIS$var)]
    ggplot2::ggplot(sf_det) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data[[v]]), colour = NA) +
      ggplot2::scale_fill_viridis_c(option = "viridis", na.value = "grey85",
                                    name = NULL) +
      ggplot2::labs(title = rot) +
      theme_mapa(base_size = 9)
  })

  fig14 <- patchwork::wrap_plots(mapas_det, ncol = 2) +
    patchwork::plot_annotation(
      title = "Determinantes socioeconômicos e de saneamento",
      caption = "Fontes: Atlas Brasil (PNUD/Ipea/FJP) e Censo Demográfico 2022 (IBGE).",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
  salvar_figura(fig14, "fig14_determinantes", largura = 11, altura = 10)
}


# =============================================================================
# TABELA CONSOLIDADA PARA O ARTIGO
# =============================================================================
# Junta os principais resultados por município — vai para material suplementar.

suplementar <- agregado_eb |>
  as.data.frame() |>
  dplyr::select(code_muni, name_muni, uf = abbrev_state, regiao,
                pop_central, casos_sinan, pessoa_ano,
                inc_sinan_bruta, inc_sinan_ebl, rme_sinan,
                examinados, positivos, posit_ebl,
                dplyr::any_of(c("idhm", "ivs", "pct_agua_rede",
                                "pct_esgoto_adeq", "pct_rural"))) |>
  dplyr::left_join(dplyr::select(as.data.frame(lisa_res), code_muni, cluster,
                                 p_lisa_aj, dplyr::any_of("Gi_cluster")),
                   by = "code_muni") |>
  dplyr::left_join(dplyr::select(as.data.frame(persist), code_muni,
                                 n_anos_AA, classe_persistencia),
                   by = "code_muni") |>
  dplyr::arrange(dplyr::desc(inc_sinan_ebl))

if (file.exists(file.path(PARAMS$dir_processados, "tendencia_municipal.rds"))) {
  tend <- readRDS(file.path(PARAMS$dir_processados, "tendencia_municipal.rds"))
  suplementar <- dplyr::left_join(
    suplementar,
    dplyr::select(tend, code_muni, RR_anual, classe_tendencia),
    by = "code_muni")
}

salvar_tabela(suplementar, "S1_resultados_municipais_completos")

message("\n=== Pipeline concluído ===")
message("Figuras: ", PARAMS$dir_figuras)
message("Tabelas: ", PARAMS$dir_tabelas)
