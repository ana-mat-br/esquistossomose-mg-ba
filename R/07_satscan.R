# =============================================================================
# 07_satscan.R — Estatística de varredura espacial e espaço-temporal
#
# DUAS ROTAS, controladas por ROTA_VARREDURA:
#
#   "smerc"   (padrão) — 100% em R, sem software externo. Varredura circular de
#                        Kulldorff, modelo Poisson, inferência Monte Carlo.
#                        Cobre a varredura PURAMENTE ESPACIAL.
#
#   "satscan" — chama o SaTScan™ via rsatscan. Necessário para a varredura
#               ESPAÇO-TEMPORAL retrospectiva (detecta clusters que emergem ou
#               desaparecem dentro de 2018-2025) e para a análise prospectiva
#               de vigilância. Requer o SaTScan instalado.
#
# A varredura de Kulldorff é complementar ao LISA, não redundante: o LISA testa
# associação local par a par; a varredura testa janelas de tamanho variável
# contra a hipótese de risco constante, com correção embutida para múltiplos
# testes — por isso é o método de escolha para DELIMITAR áreas de intervenção.
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

ROTA_VARREDURA <- "smerc"   # "smerc" | "satscan" | "ambas"

agregado_eb <- readRDS(file.path(PARAMS$dir_processados, "agregado_eb.rds"))
painel      <- readRDS(file.path(PARAMS$dir_processados, "painel.rds"))
malha_est   <- readRDS(file.path(PARAMS$dir_processados, "malha_est.rds"))

stopifnot(identical(agregado_eb$code_muni, malha_est$code_muni))

# Centroides em coordenadas geográficas (SaTScan usa lat/long; smerc aceita
# ambos, com longlat = TRUE para distância em km sobre a esfera)
suppressWarnings({
  cent <- sf::st_coordinates(
    sf::st_point_on_surface(sf::st_transform(sf::st_geometry(malha_est), 4326))
  )
})
colnames(cent) <- c("lon", "lat")


# =============================================================================
# ROTA A — smerc (varredura puramente espacial)
# =============================================================================

varredura_smerc <- function(casos, pop, coords, rotulo,
                            ubpop = PARAMS$satscan_ubpop,
                            nsim  = PARAMS$satscan_nsim) {

  if (!requireNamespace("smerc", quietly = TRUE)) {
    stop("Instale o pacote 'smerc'.")
  }
  ok <- is.finite(casos) & is.finite(pop) & pop > 0
  if (sum(ok) < 20) {
    warning("Dados insuficientes para varredura: ", rotulo)
    return(NULL)
  }

  message(glue::glue("\n--- Varredura de Kulldorff ({rotulo}) ---"))
  esperado <- sum(casos[ok]) / sum(pop[ok]) * pop[ok]

  out <- smerc::scan.test(
    coords  = coords[ok, , drop = FALSE],
    cases   = as.integer(round(casos[ok])),
    pop     = pop[ok],
    ex      = esperado,
    nsim    = nsim,
    alpha   = PARAMS$alfa,
    ubpop   = ubpop,
    longlat = TRUE      # distância geodésica em km
  )

  idx_orig <- which(ok)

  # Campos devolvidos por smerc::scan.test em cada cluster:
  #   locids, centroid, r, max_dist, population, cases, expected, smr, rr,
  #   loglikrat, test_statistic, pvalue, w
  cl <- purrr::imap_dfr(out$clusters, function(cc, i) {
    tibble::tibble(
      cluster    = i,
      tipo       = if (i == 1) "Primário" else "Secundário",
      n_municipios = length(cc$locids),
      centroide_idx = idx_orig[cc$locids[1]],
      raio_km    = cc$max_dist,
      casos_obs  = cc$cases,
      casos_esp  = cc$expected,
      populacao  = cc$population,
      RME        = cc$smr,        # observado/esperado dentro do cluster
      RR         = cc$rr,         # risco dentro vs. fora do cluster
      LLR        = cc$loglikrat,
      p_valor    = cc$pvalue
    )
  }) |>
    dplyr::filter(p_valor <= PARAMS$alfa)

  if (nrow(cl) == 0) {
    message("Nenhum cluster significativo (", rotulo, ").")
    return(list(tabela = cl, membros = tibble::tibble(), obj = out))
  }

  membros <- purrr::imap_dfr(out$clusters, function(cc, i) {
    if (cc$pvalue > PARAMS$alfa) return(NULL)
    tibble::tibble(cluster = i, idx = idx_orig[cc$locids])
  }) |>
    dplyr::mutate(
      code_muni = malha_est$code_muni[idx],
      name_muni = malha_est$name_muni[idx],
      uf        = malha_est$abbrev_state[idx]
    )

  cl <- cl |>
    dplyr::mutate(
      municipio_centro = malha_est$name_muni[centroide_idx],
      uf_centro        = malha_est$abbrev_state[centroide_idx],
      desfecho         = rotulo
    )

  message(glue::glue("{nrow(cl)} cluster(s) significativo(s); ",
                     "primário: RR = {round(cl$RR[1],2)}, p = {cl$p_valor[1]}"))
  list(tabela = cl, membros = membros, obj = out)
}

if (ROTA_VARREDURA %in% c("smerc", "ambas")) {

  set.seed(PARAMS$semente)

  # Desfecho 1 — incidência SINAN (denominador: pessoa-ano)
  vs_sinan <- varredura_smerc(agregado_eb$casos_sinan, agregado_eb$pessoa_ano,
                              cent, "Incidência SINAN")

  # Desfecho 2 — positividade PCE (denominador: examinados)
  tem_pce <- agregado_eb$examinados > 0
  vs_pce <- varredura_smerc(
    ifelse(tem_pce, agregado_eb$positivos, NA),
    ifelse(tem_pce, agregado_eb$examinados, NA),
    cent, "Positividade PCE"
  )

  # varredura_smerc() devolve NULL quando o desfecho não tem dados (ex.: PCE
  # ausente). Guardar contra isso mantém o pipeline rodando com um só desfecho.
  marcar <- function(res, rotulo) {
    if (is.null(res) || is.null(res$membros) || nrow(res$membros) == 0) return(NULL)
    dplyr::mutate(res$membros, desfecho = rotulo)
  }
  tabela_clusters <- dplyr::bind_rows(
    if (!is.null(vs_sinan)) vs_sinan$tabela,
    if (!is.null(vs_pce))   vs_pce$tabela
  )
  membros_clusters <- dplyr::bind_rows(
    marcar(vs_sinan, "Incidência SINAN"),
    marcar(vs_pce,   "Positividade PCE")
  )
  if (is.null(membros_clusters)) membros_clusters <- tibble::tibble()

  if (nrow(tabela_clusters) > 0) {
    tabela_clusters <- tabela_clusters |>
      dplyr::mutate(
        RR_fmt = sprintf("%.2f", RR),
        p_fmt  = fmt_p(p_valor)
      ) |>
      dplyr::select(desfecho, cluster, tipo, municipio_centro, uf_centro,
                    n_municipios, raio_km, populacao, casos_obs, casos_esp,
                    RR_fmt, LLR, p_fmt, dplyr::everything())
    salvar_tabela(tabela_clusters,  "16_clusters_varredura")
    salvar_tabela(membros_clusters, "17_municipios_por_cluster")
  }

  saveRDS(list(tabela = tabela_clusters, membros = membros_clusters),
          file.path(PARAMS$dir_processados, "clusters_varredura.rds"))

  # --- Perfil socioeconômico: dentro vs. fora do cluster primário -----------
  vars_ctx_cl <- intersect(
    c("idhm", "ivs", "pct_agua_rede", "pct_esgoto_adeq", "pct_agua_superficial",
      "pct_esgoto_em_corpo_dagua", "pct_lixo_coletado", "pct_rural"),
    names(agregado_eb))

  if (nrow(membros_clusters) > 0) {
    dentro <- membros_clusters |>
      dplyr::filter(desfecho == "Incidência SINAN", cluster == 1) |>
      dplyr::pull(code_muni)

    perfil <- agregado_eb |>
      dplyr::mutate(no_cluster = dplyr::if_else(code_muni %in% dentro,
                                                "Dentro do cluster primário",
                                                "Fora")) |>
      dplyr::group_by(no_cluster) |>
      dplyr::summarise(
        n = dplyr::n(),
        inc_mediana = median(inc_sinan_ebl, na.rm = TRUE),
        dplyr::across(dplyr::all_of(vars_ctx_cl), ~ median(.x, na.rm = TRUE)),
        .groups = "drop"
      )
    salvar_tabela(perfil, "18_perfil_cluster_primario")
    print(as.data.frame(perfil))

    # Teste formal de diferença
    teste_dentro_fora <- purrr::map_dfr(
      vars_ctx_cl,
      function(v) {
        g <- agregado_eb$code_muni %in% dentro
        if (sum(!is.na(agregado_eb[[v]][g])) < 3) return(NULL)
        tt <- wilcox.test(agregado_eb[[v]][g], agregado_eb[[v]][!g], exact = FALSE)
        tibble::tibble(variavel = v, p_valor = tt$p.value)
      }) |>
      dplyr::mutate(p_aj = p.adjust(p_valor, "fdr"), p_fmt = fmt_p(p_aj))
    salvar_tabela(teste_dentro_fora, "18b_teste_cluster_vs_fora")
  }
}


# =============================================================================
# ROTA B — SaTScan (espacial + espaço-temporal retrospectivo)
# =============================================================================

preparar_arquivos_satscan <- function(painel, malha, cent,
                                      dir = PARAMS$dir_satscan) {

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  # .cas — id, casos, ano
  cas <- painel |>
    dplyr::transmute(locationid = code_muni,
                     numcases   = as.integer(casos_sinan),
                     year       = as.integer(ano)) |>
    dplyr::filter(!is.na(numcases))

  # .pop — id, ano, população
  pop <- painel |>
    dplyr::transmute(locationid = code_muni,
                     year       = as.integer(ano),
                     population = round(populacao)) |>
    dplyr::filter(!is.na(population))

  # .geo — id, latitude, longitude
  geo <- tibble::tibble(locationid = malha$code_muni,
                        latitude   = cent[, "lat"],
                        longitude  = cent[, "lon"])

  list(cas = as.data.frame(cas), pop = as.data.frame(pop),
       geo = as.data.frame(geo), dir = dir)
}

rodar_satscan <- function(tipo_analise = c("espacial", "espaco_temporal",
                                           "prospectiva"),
                          arquivos,
                          nome = "esquisto") {

  tipo_analise <- match.arg(tipo_analise)
  if (!requireNamespace("rsatscan", quietly = TRUE)) {
    stop("Instale: install.packages('rsatscan') e o software SaTScan ",
         "(https://www.satscan.org).")
  }
  if (!dir.exists(PARAMS$satscan_path)) {
    stop(glue::glue("SaTScan não encontrado em {PARAMS$satscan_path}. ",
                    "Ajuste PARAMS$satscan_path em R/00_setup.R."))
  }

  # Códigos SaTScan:
  #   AnalysisType : 1=puramente espacial, 3=espaço-temporal retrospectivo,
  #                  4=espaço-temporal prospectivo
  #   ModelType    : 0=Poisson discreto
  #   ScanAreas    : 1=taxas altas
  #   PrecisionCaseTimes / TimeAggregationUnits : 1=ano
  cod_analise <- switch(tipo_analise, espacial = 1, espaco_temporal = 3,
                        prospectiva = 4)

  rsatscan::ss.options(reset = TRUE)
  rsatscan::ss.options(list(
    CaseFile             = paste0(nome, ".cas"),
    PopulationFile       = paste0(nome, ".pop"),
    CoordinatesFile      = paste0(nome, ".geo"),
    PrecisionCaseTimes   = 1,
    StartDate            = as.character(PARAMS$ano_ini),
    EndDate              = as.character(PARAMS$ano_fim),
    CoordinatesType      = 1,          # latitude/longitude
    AnalysisType         = cod_analise,
    ModelType            = 0,
    ScanAreas            = 1,
    TimeAggregationUnits = 1,
    TimeAggregationLength = 1,
    MaxSpatialSizeInPopulationAtRisk = PARAMS$satscan_ubpop * 100,
    MaxTemporalSizeInterpretation = 1,
    MaxTemporalSize      = 50,
    MonteCarloReps       = PARAMS$satscan_nsim,
    ReportGiniClusters   = "n",
    NonCompactnessPenalty = 1,         # penaliza clusters muito alongados
    SpatialWindowShapeType = 0,        # 0=circular; 1=elíptico
    ResultsFile          = file.path(arquivos$dir, paste0(nome, "_", tipo_analise, ".txt")),
    OutputGoogleEarthKML = "n",
    MostLikelyClusterEachCentroidASCII = "y",
    MostLikelyClusterCaseInfoEachCentASCII = "y",
    CensusAreasReportedClustersASCII = "y"
  ))

  rsatscan::write.cas(arquivos$cas, arquivos$dir, nome)
  rsatscan::write.pop(arquivos$pop, arquivos$dir, nome)
  rsatscan::write.geo(arquivos$geo, arquivos$dir, nome)
  rsatscan::write.ss.prm(arquivos$dir, nome)

  message("Executando SaTScan — ", tipo_analise, " ...")
  res <- rsatscan::satscan(
    prmlocation   = arquivos$dir,
    prmfilename   = nome,
    sslocation    = PARAMS$satscan_path,
    ssbatchfilename = "satscan",
    verbose = FALSE
  )
  res
}

if (ROTA_VARREDURA %in% c("satscan", "ambas")) {

  arq <- preparar_arquivos_satscan(painel, malha_est, cent)

  res_esp <- try(rodar_satscan("espacial", arq, "esquisto_esp"), silent = TRUE)
  res_st  <- try(rodar_satscan("espaco_temporal", arq, "esquisto_st"), silent = TRUE)

  extrair_clusters_ss <- function(res, rotulo) {
    if (inherits(res, "try-error") || is.null(res$col)) return(NULL)
    res$col |>
      janitor::clean_names() |>
      dplyr::mutate(analise = rotulo, .before = 1)
  }

  clusters_ss <- dplyr::bind_rows(
    extrair_clusters_ss(res_esp, "Puramente espacial"),
    extrair_clusters_ss(res_st,  "Espaço-temporal retrospectiva")
  )

  if (!is.null(clusters_ss) && nrow(clusters_ss) > 0) {
    salvar_tabela(clusters_ss, "19_clusters_satscan")
    # .gis traz o município-a-cluster
    gis <- dplyr::bind_rows(
      if (!inherits(res_esp, "try-error") && !is.null(res_esp$gis))
        dplyr::mutate(janitor::clean_names(res_esp$gis), analise = "Puramente espacial"),
      if (!inherits(res_st, "try-error") && !is.null(res_st$gis))
        dplyr::mutate(janitor::clean_names(res_st$gis), analise = "Espaço-temporal")
    )
    if (nrow(gis) > 0) salvar_tabela(gis, "20_municipios_cluster_satscan")
    saveRDS(list(esp = res_esp, st = res_st),
            file.path(PARAMS$dir_modelos, "satscan.rds"))
  } else {
    message("SaTScan não retornou clusters (ou falhou). Verifique ",
            file.path(arq$dir, "*.txt"))
  }
}


# =============================================================================
# GEOMETRIA DOS CLUSTERS (para os mapas do script 11)
# =============================================================================

if (exists("membros_clusters") && nrow(membros_clusters) > 0) {
  clusters_sf <- membros_clusters |>
    dplyr::left_join(dplyr::select(malha_est, code_muni), by = "code_muni") |>
    sf::st_as_sf() |>
    dplyr::group_by(desfecho, cluster) |>
    dplyr::summarise(n_mun = dplyr::n(), .groups = "drop")
  sf::st_write(clusters_sf,
               file.path(PARAMS$dir_processados, "clusters.gpkg"),
               delete_dsn = TRUE, quiet = TRUE)
  message("Geometria dos clusters salva em dados/processados/clusters.gpkg")
}

message("\nVarredura concluída. Próximo: R/08_poisson_robusta.R")
