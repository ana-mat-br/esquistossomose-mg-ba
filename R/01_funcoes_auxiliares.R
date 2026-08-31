# =============================================================================
# 01_funcoes_auxiliares.R — Funções reutilizadas por todos os scripts
# =============================================================================

# `%||%` existe no R base a partir da 4.4.0; definimos para compatibilidade
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a


# ---- 1. Códigos de município ------------------------------------------------

#' Converte código IBGE de 6 para 7 dígitos (recupera o dígito verificador)
#'
#' O SINAN grava o município de residência com 6 dígitos; IBGE/geobr usam 7.
#' A conversão é feita por *lookup* na tabela de municípios (não por cálculo do
#' DV), evitando erros em municípios extintos/criados no período.
cod6_para_cod7 <- function(cod6, tabela_muni) {
  stopifnot(all(c("code_muni") %in% names(tabela_muni)))
  chave <- tibble::tibble(
    cod6 = substr(as.character(tabela_muni$code_muni), 1, 6),
    cod7 = as.character(tabela_muni$code_muni)
  ) |> dplyr::distinct()

  out <- chave$cod7[match(as.character(cod6), chave$cod6)]
  n_na <- sum(is.na(out) & !is.na(cod6))
  if (n_na > 0) {
    warning(glue::glue("{n_na} códigos de 6 dígitos sem correspondência na malha."))
  }
  out
}

#' Padroniza qualquer coluna de código municipal para character de 7 dígitos
pad_cod_muni <- function(x) {
  x <- gsub("[^0-9]", "", as.character(x))
  ifelse(nchar(x) == 6, NA_character_, sprintf("%07s", x)) |>
    (\(z) gsub(" ", "0", z))()
}


# ---- 2. Vizinhança espacial -------------------------------------------------

#' Constrói lista de vizinhança e matriz de pesos espaciais
#'
#' Trata explicitamente ilhas (municípios sem vizinho por contiguidade),
#' que quebram Moran/SAR. Se `tipo = "queen"` e houver ilhas, elas são ligadas
#' ao vizinho mais próximo (relatado na saída).
#'
#' @return lista com `nb`, `listw` e `n_ilhas`
construir_vizinhanca <- function(sf_obj,
                                 tipo  = c("queen", "rook", "knn"),
                                 k     = 5,
                                 style = "W") {
  tipo <- match.arg(tipo)
  stopifnot(inherits(sf_obj, "sf"))

  # Coordenadas dos centroides (para knn e para religar ilhas).
  # st_point_on_surface garante ponto interno em polígonos côncavos.
  suppressWarnings({
    coords <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(sf_obj)))
  })

  if (tipo == "knn") {
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k), sym = TRUE)
    n_ilhas <- 0L
  } else {
    nb <- spdep::poly2nb(sf_obj, queen = (tipo == "queen"))
    ilhas <- which(spdep::card(nb) == 0)
    n_ilhas <- length(ilhas)
    if (n_ilhas > 0) {
      message(glue::glue(
        "{n_ilhas} município(s) sem vizinho por contiguidade -> ligado(s) ao vizinho mais próximo."
      ))
      nb_knn1 <- spdep::knn2nb(spdep::knearneigh(coords, k = 1))
      for (i in ilhas) {
        viz <- nb_knn1[[i]]
        nb[[i]] <- as.integer(viz)
        nb[[viz]] <- sort(unique(c(nb[[viz]][nb[[viz]] != 0L], as.integer(i))))
      }
      nb <- spdep::make.sym.nb(nb)
    }
  }

  listw <- spdep::nb2listw(nb, style = style, zero.policy = TRUE)

  list(nb = nb, listw = listw, n_ilhas = n_ilhas, coords = coords)
}


# ---- 3. Estimadores bayesianos empíricos ------------------------------------

#' Estimador bayesiano empírico GLOBAL (Marshall, 1991)
#'
#' Encolhe a taxa bruta municipal em direção à média do conjunto todo.
#' @param casos numerador; @param pop denominador; @param base multiplicador
eb_global <- function(casos, pop, base = 1e5) {
  eb <- spdep::EBest(n = casos, x = pop)
  tibble::tibble(
    taxa_bruta = eb$raw * base,
    taxa_ebg   = eb$estmm * base
  )
}

#' Estimador bayesiano empírico LOCAL (Marshall, 1991)
#'
#' Encolhe a taxa bruta em direção à média da VIZINHANÇA, corrigindo a
#' instabilidade de taxas em municípios de população pequena sem apagar
#' heterogeneidade regional. É o estimador indicado quando há autocorrelação
#' espacial — exatamente o caso do corredor endêmico MG/BA.
#'
#' @param nb objeto nb (de construir_vizinhanca)
eb_local <- function(casos, pop, nb, base = 1e5) {
  eb <- spdep::EBlocal(ri = casos, ni = pop, nb = nb, zero.policy = TRUE)
  tibble::tibble(
    taxa_bruta = eb$raw * base,
    taxa_ebl   = eb$est * base
  )
}


# ---- 4. Moran / LISA --------------------------------------------------------

#' I de Moran global com inferência por permutação
moran_global <- function(x, listw, nsim = 9999) {
  ok <- is.finite(x)
  if (any(!ok)) warning(sum(!ok), " valor(es) não finito(s) removido(s) do Moran global.")
  mc <- spdep::moran.mc(x, listw, nsim = nsim, zero.policy = TRUE, na.action = na.exclude)
  tibble::tibble(
    I         = unname(mc$statistic),
    I_esperado = -1 / (length(x) - 1),
    p_valor   = mc$p.value,
    n_sim     = nsim
  )
}

#' LISA (Moran local) com permutação condicional e classificação em quadrantes
#'
#' Usa `localmoran_perm`, cujo p-valor por permutação é preferível ao analítico
#' (o teste analítico assume normalidade, violada por taxas de doença).
#' A classificação de quadrantes vem do atributo `quadr` (comparação com a
#' MEDIANA — mais robusta a assimetria do que a média).
#'
#' @param ajuste_p "fdr", "bonferroni" ou "none"
lisa <- function(x, listw, nsim = 9999, alfa = 0.05, ajuste_p = "fdr") {

  lm_perm <- spdep::localmoran_perm(x, listw, nsim = nsim, zero.policy = TRUE,
                                    na.action = na.exclude)
  quad <- attr(lm_perm, "quadr")

  # Coluna de p-valor por permutação ("folded", 2 caudas) é a última nomeada
  # com "Sim"; usamos busca por nome para não depender da posição.
  nm <- colnames(lm_perm)
  col_p <- nm[grepl("Pr\\(folded\\)", nm)]
  if (length(col_p) == 0) col_p <- nm[grepl("^Pr", nm)][1]
  p_bruto <- as.numeric(lm_perm[, col_p])

  p_aj <- if (ajuste_p == "none") p_bruto else p.adjust(p_bruto, method = ajuste_p)

  cluster <- as.character(quad$median)
  cluster <- dplyr::recode(cluster,
    "High-High" = "Alto-Alto", "Low-Low" = "Baixo-Baixo",
    "High-Low"  = "Alto-Baixo", "Low-High" = "Baixo-Alto",
    .default = "Não significativo"
  )
  cluster[is.na(p_aj) | p_aj > alfa] <- "Não significativo"

  tibble::tibble(
    Ii       = as.numeric(lm_perm[, "Ii"]),
    Z_Ii     = as.numeric(lm_perm[, "Z.Ii"]),
    p_lisa   = p_bruto,
    p_lisa_aj = p_aj,
    cluster  = factor(cluster, levels = names(PAL_LISA))
  )
}

#' Gi* de Getis-Ord (hot/cold spots), complementar ao LISA
getis_ord <- function(x, nb, listw_style = "B", alfa = 0.05, ajuste_p = "fdr") {
  # Gi* inclui o próprio polígono na vizinhança
  nb_self <- spdep::include.self(nb)
  lw_self <- spdep::nb2listw(nb_self, style = listw_style, zero.policy = TRUE)
  g <- spdep::localG(x, lw_self, zero.policy = TRUE)
  z <- as.numeric(g)
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  p_aj <- if (ajuste_p == "none") p else p.adjust(p, method = ajuste_p)
  cls <- dplyr::case_when(
    p_aj <= alfa & z > 0 ~ "Hot spot",
    p_aj <= alfa & z < 0 ~ "Cold spot",
    TRUE                 ~ "Não significativo"
  )
  tibble::tibble(Gi_z = z, Gi_p = p, Gi_p_aj = p_aj,
                 Gi_cluster = factor(cls, levels = c("Hot spot", "Cold spot",
                                                     "Não significativo")))
}


# ---- 5. Modelos com variância robusta ---------------------------------------

#' Extrai RR (IC95%) de um GLM Poisson/binomial-log com variância robusta
#'
#' Implementa a "regressão de Poisson modificada" (Zou, 2004): estimativa
#' pontual pelo Poisson (consistente para o log da razão) + erro-padrão
#' sanduíche, que corrige a superestimação da variância do Poisson quando o
#' desfecho é proporção, e a subestimação quando há superdispersão.
#'
#' @param tipo_vcov "HC0" (Zou clássico), "HC3" (melhor em n pequeno) ou
#'   "cluster" (agrupado por `cluster_var`, obrigatório em dados de painel)
#' @param cluster_var vetor de identificadores de cluster (ex.: código municipal)
tidy_rr_robusto <- function(modelo,
                            tipo_vcov   = c("HC0", "HC3", "cluster"),
                            cluster_var = NULL,
                            conf_level  = 0.95,
                            exponenciar = TRUE) {
  tipo_vcov <- match.arg(tipo_vcov)

  V <- switch(tipo_vcov,
    "HC0"     = sandwich::vcovHC(modelo, type = "HC0"),
    "HC3"     = sandwich::vcovHC(modelo, type = "HC3"),
    "cluster" = {
      if (is.null(cluster_var)) stop("cluster_var é obrigatório para tipo_vcov='cluster'.")
      sandwich::vcovCL(modelo, cluster = cluster_var, type = "HC0")
    }
  )

  ct <- lmtest::coeftest(modelo, vcov. = V)
  z  <- qnorm(1 - (1 - conf_level) / 2)

  res <- tibble::tibble(
    termo    = rownames(ct),
    estimativa = ct[, 1],
    ep_robusto = ct[, 2],
    estat_z    = ct[, 3],
    p_valor    = ct[, 4],
    li         = ct[, 1] - z * ct[, 2],
    ls         = ct[, 1] + z * ct[, 2]
  )

  if (exponenciar) {
    res <- res |>
      dplyr::mutate(dplyr::across(c(estimativa, li, ls), exp)) |>
      dplyr::rename(RR = estimativa, RR_li = li, RR_ls = ls)
  }
  res
}

#' Formata RR (IC95%) para tabela de artigo
fmt_rr <- function(rr, li, ls, dig = 2) {
  sprintf(paste0("%.", dig, "f (%.", dig, "f–%.", dig, "f)"), rr, li, ls)
}

#' Formata p-valor no padrão de periódico
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ NA_character_,
    p < 0.001  ~ "<0,001",
    TRUE       ~ sub("\\.", ",", sprintf("%.3f", p))
  )
}


# ---- 6. Diagnósticos --------------------------------------------------------

#' Teste de superdispersão (Cameron & Trivedi) + razão deviance/gl
diag_dispersao <- function(modelo_poisson) {
  disp_ratio <- sum(residuals(modelo_poisson, type = "pearson")^2) /
    df.residual(modelo_poisson)
  ct <- try(AER::dispersiontest(modelo_poisson, trafo = 1), silent = TRUE)
  tibble::tibble(
    razao_pearson_gl = disp_ratio,
    deviance_gl      = modelo_poisson$deviance / df.residual(modelo_poisson),
    alfa_CT          = if (inherits(ct, "try-error")) NA_real_ else unname(ct$estimate),
    p_CT             = if (inherits(ct, "try-error")) NA_real_ else ct$p.value,
    veredito = dplyr::case_when(
      disp_ratio > 1.5 ~ "Superdispersão relevante — usar NB ou variância robusta",
      disp_ratio < 0.7 ~ "Subdispersão — variância robusta recomendada",
      TRUE             ~ "Dispersão aceitável"
    )
  )
}

#' Multicolinearidade: VIF (>5 alerta, >10 crítico)
diag_vif <- function(modelo) {
  v <- try(car::vif(modelo), silent = TRUE)
  if (inherits(v, "try-error")) return(tibble::tibble(termo = NA_character_, VIF = NA_real_))
  if (is.matrix(v)) v <- v[, "GVIF^(1/(2*Df))"]^2
  tibble::tibble(termo = names(v), VIF = as.numeric(v)) |>
    dplyr::mutate(alerta = dplyr::case_when(VIF > 10 ~ "crítico",
                                            VIF > 5  ~ "atenção",
                                            TRUE     ~ "ok")) |>
    dplyr::arrange(dplyr::desc(VIF))
}

#' Testes de multiplicador de Lagrange para dependência espacial em resíduos OLS
#'
#' Compatível com spdep antigo (lm.LMtests) e >= 1.3 (lm.RStests).
teste_ml_espacial <- function(modelo_ols, listw) {
  fun <- if ("lm.RStests" %in% getNamespaceExports("spdep")) {
    getExportedValue("spdep", "lm.RStests")
  } else {
    getExportedValue("spdep", "lm.LMtests")
  }
  testes <- c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA")
  out <- try(fun(modelo_ols, listw, test = testes, zero.policy = TRUE), silent = TRUE)
  if (inherits(out, "try-error")) {
    # nomenclatura nova: RSerr, RSlag, adjRSerr, adjRSlag, SARMA
    out <- fun(modelo_ols, listw,
               test = c("RSerr", "RSlag", "adjRSerr", "adjRSlag", "SARMA"),
               zero.policy = TRUE)
  }
  purrr::map_dfr(out, ~ tibble::tibble(estatistica = unname(.x$statistic),
                                       gl = unname(.x$parameter),
                                       p_valor = .x$p.value),
                 .id = "teste")
}


# ---- 7. Transformações de covariáveis ---------------------------------------

#' Aplica as transformações declaradas em COVARIAVEIS (script 00)
aplicar_transformacoes <- function(dados, dicionario = COVARIAVEIS) {
  dic <- dicionario |> dplyr::filter(var %in% names(dados))
  for (i in seq_len(nrow(dic))) {
    v <- dic$var[i]
    dados[[v]] <- switch(dic$transf[i],
      "z"    = as.numeric(scale(dados[[v]])),
      "pp10" = dados[[v]] / 10,
      "log"   = log(dados[[v]] + 1),
      "log1p" = log1p(dados[[v]]),
      "id"   = dados[[v]],
      dados[[v]]
    )
  }
  attr(dados, "transformacoes") <- dic
  dados
}

#' Rótulo interpretativo do RR conforme a transformação aplicada
rotulo_incremento <- function(var, dicionario = COVARIAVEIS) {
  tr <- dicionario$transf[match(var, dicionario$var)]
  dplyr::case_when(
    tr == "z"    ~ "por 1 DP de aumento",
    tr == "pp10" ~ "por 10 p.p. de aumento",
    tr == "log"   ~ "por aumento de 1 unidade em log",
    tr == "log1p" ~ "por duplicação aproximada do percentual",
    TRUE         ~ "por 1 unidade"
  )
}


# ---- 8. Utilitários ---------------------------------------------------------

#' Salva tabela em CSV (padrão brasileiro) e XLSX
salvar_tabela <- function(x, nome, dir = PARAMS$dir_tabelas) {
  readr::write_csv2(x, file.path(dir, paste0(nome, ".csv")), na = "")
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(x, file.path(dir, paste0(nome, ".xlsx")))
  }
  message("Tabela salva: ", file.path(dir, nome))
  invisible(x)
}

#' Salva figura em PNG (alta resolução) e PDF vetorial
salvar_figura <- function(plot, nome, largura = 9, altura = 7,
                          dir = PARAMS$dir_figuras, dpi = 320) {
  ggplot2::ggsave(file.path(dir, paste0(nome, ".png")), plot,
                  width = largura, height = altura, dpi = dpi, bg = "white")
  # cairo_pdf exige X11/cairo, ausente em muitas instalações macOS.
  # Tenta cairo e cai para o device pdf padrão se indisponível.
  ok <- try(ggplot2::ggsave(file.path(dir, paste0(nome, ".pdf")), plot,
              width = largura, height = altura, device = grDevices::cairo_pdf),
            silent = TRUE)
  if (inherits(ok, "try-error")) {
    ggplot2::ggsave(file.path(dir, paste0(nome, ".pdf")), plot,
                    width = largura, height = altura)
  }
  message("Figura salva: ", file.path(dir, nome))
  invisible(plot)
}

#' Checagem defensiva do painel antes de modelar
checar_painel <- function(df, id = "code_muni", tempo = "ano") {
  n_mun  <- dplyr::n_distinct(df[[id]])
  n_ano  <- dplyr::n_distinct(df[[tempo]])
  message(glue::glue(
    "Painel: {n_mun} municípios x {n_ano} anos = {n_mun * n_ano} esperado; ",
    "{nrow(df)} observado."
  ))
  dup <- df |> dplyr::count(.data[[id]], .data[[tempo]]) |> dplyr::filter(n > 1)
  if (nrow(dup) > 0) warning(nrow(dup), " combinações município-ano duplicadas.")
  na_tab <- df |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.x)))) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "variavel",
                        values_to = "n_faltantes") |>
    dplyr::filter(n_faltantes > 0) |>
    dplyr::arrange(dplyr::desc(n_faltantes))
  if (nrow(na_tab) > 0) {
    message("Variáveis com faltantes:")
    print(as.data.frame(na_tab))
  }
  invisible(list(duplicados = dup, faltantes = na_tab))
}
