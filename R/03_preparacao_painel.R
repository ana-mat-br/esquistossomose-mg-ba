# =============================================================================
# 03_preparacao_painel.R — Limpeza, harmonização e montagem do painel
#                          município x ano (2018-2025)
#
# Produz:
#   dados/processados/painel.rds     — data.frame do painel (long)
#   dados/processados/malha_est.rds  — sf da área de estudo
#   saidas/tabelas/00_fluxograma_dados.csv — rastro de exclusões (para o artigo)
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")
source("R/02_obtencao_dados.R")

# Registro de exclusões — vai virar o fluxograma de dados do artigo
FLUXO <- tibble::tibble(etapa = character(), n = numeric(), obs = character())
registrar <- function(etapa, n, obs = "") {
  FLUXO <<- dplyr::add_row(FLUXO, etapa = etapa, n = n, obs = obs)
  message(glue::glue("[fluxo] {etapa}: {n}  {obs}"))
}


# =============================================================================
# 1. ÁREA DE ESTUDO
# =============================================================================

malha         <- obter_malha()
cods_norte_mg <- readRDS(file.path(PARAMS$dir_processados, "cods_norte_mg.rds"))

malha_est <- malha |>
  dplyr::filter(abbrev_state == "BA" | code_muni %in% cods_norte_mg) |>
  dplyr::mutate(
    regiao = dplyr::if_else(abbrev_state == "BA", "Bahia", "Norte de Minas"),
    area_km2 = as.numeric(sf::st_area(sf::st_transform(geom, 5880))) / 1e6
  ) |>
  dplyr::arrange(code_muni)

registrar("Municípios na área de estudo", nrow(malha_est),
          glue::glue("BA={sum(malha_est$regiao=='Bahia')}, ",
                     "N-MG={sum(malha_est$regiao=='Norte de Minas')}"))

cods_estudo <- malha_est$code_muni


# =============================================================================
# 2. SINAN — casos notificados
# =============================================================================

#' Limpa o SINAN e agrega em casos por município de residência x ano
#'
#' Decisões (todas explicitáveis no artigo):
#'   - unidade = município de RESIDÊNCIA (ID_MN_RESI), não de notificação
#'   - ano     = ano do primeiro sintoma se disponível; senão, ano de notificação
#'   - casos confirmados: CLASSI_FIN == 1 quando a variável existir
#'   - duplicatas removidas por (NU_NOTIFIC, ID_MUNICIP, DT_NOTIFIC)
preparar_sinan <- function(sinan_bruto, cods_estudo, malha) {

  nm <- names(sinan_bruto)
  pega <- function(...) {
    alvos <- tolower(c(...))
    hit <- nm[tolower(nm) %in% alvos]
    if (length(hit)) hit[1] else NA_character_
  }

  v_resi   <- pega("id_mn_resi", "id_mun_resi", "co_mun_res", "municipio_residencia")
  v_ano    <- pega("nu_ano", "ano", "ano_notif")
  v_dtnot  <- pega("dt_notific", "dt_notifica")
  v_dtsin  <- pega("dt_sin_pri", "dt_sintoma")
  v_class  <- pega("classi_fin", "classificacao_final")
  v_notif  <- pega("nu_notific", "nu_notifica")
  v_munnot <- pega("id_municip", "id_mn_notif")

  if (is.na(v_resi)) {
    stop("Coluna de município de residência não localizada no SINAN. ",
         "Colunas disponíveis: ", paste(head(nm, 40), collapse = ", "))
  }

  d <- sinan_bruto
  registrar("SINAN — registros brutos", nrow(d))

  # -- deduplicação ---------------------------------------------------------
  # O SINAN de esquistossomose NÃO traz NU_NOTIFIC. Deduplicar só por
  # (município, data de notificação) colapsaria casos genuínos ocorridos no
  # mesmo município e no mesmo dia — erro grave em município endêmico, onde
  # notificações em bloco após inquérito são a regra.
  #
  # Sem identificador de registro, usamos uma identidade individual aproximada
  # (local + data + atributos da pessoa). Se não houver campos suficientes para
  # discriminar indivíduos, NÃO deduplicamos e registramos a decisão.
  if (!is.na(v_notif)) {
    chaves <- na.omit(c(v_notif, v_munnot, v_dtnot))
  } else {
    candidatos <- c(v_dtnot, v_munnot, v_resi, v_dtsin,
                    pega("dt_nasc"), pega("ano_nasc"), pega("cs_sexo"),
                    pega("nu_idade_n"), pega("cs_raca"), pega("id_ocupa_n"))
    chaves <- unique(na.omit(candidatos))
  }

  if (length(chaves) >= 5) {
    n_antes <- nrow(d)
    d <- dplyr::distinct(d, dplyr::across(dplyr::all_of(chaves)), .keep_all = TRUE)
    registrar("SINAN — após remoção de duplicatas", nrow(d),
              glue::glue("{n_antes - nrow(d)} removidos; chaves: ",
                         "{paste(chaves, collapse='+')}"))
  } else {
    registrar("SINAN — deduplicação NÃO aplicada", nrow(d),
              glue::glue("apenas {length(chaves)} campo(s) discriminante(s) ",
                         "disponível(is) — risco de remover casos verdadeiros"))
  }

  # -- ano de referência ----------------------------------------------------
  data_ref <- if (!is.na(v_dtsin)) d[[v_dtsin]] else if (!is.na(v_dtnot)) d[[v_dtnot]] else NA
  ano_calc <- suppressWarnings(lubridate::year(lubridate::ymd(data_ref)))
  if (!is.na(v_ano)) {
    ano_decl <- suppressWarnings(as.integer(d[[v_ano]]))
    ano_final <- dplyr::coalesce(ano_calc, ano_decl)
  } else {
    ano_final <- ano_calc
  }
  d$.ano <- ano_final

  # -- município de residência (SINAN grava 6 dígitos) ----------------------
  cod_raw <- gsub("[^0-9]", "", as.character(d[[v_resi]]))
  d$.cod7 <- dplyr::case_when(
    nchar(cod_raw) == 7 ~ cod_raw,
    nchar(cod_raw) == 6 ~ cod6_para_cod7(cod_raw, malha),
    TRUE                ~ NA_character_
  )

  n_sem_cod <- sum(is.na(d$.cod7))
  if (n_sem_cod > 0) {
    registrar("SINAN — código de residência ausente/ignorado (não localizado na lista nacional do IBGE)",
              n_sem_cod)
  }
  registrar("SINAN — com município de residência válido (todo o Brasil)",
            sum(!is.na(d$.cod7)))

  # -- confirmação diagnóstica ---------------------------------------------
  if (!is.na(v_class)) {
    cf <- as.character(d[[v_class]])
    n_antes <- nrow(d)
    d <- d[cf %in% c("1", "01") | is.na(cf), , drop = FALSE]
    registrar("SINAN — após manter apenas casos confirmados",
              nrow(d), glue::glue("excluídos {n_antes - nrow(d)} descartados/inconclusivos"))
  } else {
    message("AVISO: CLASSI_FIN ausente — todos os registros tratados como confirmados.")
  }

  # -- recorte espaço-temporal ---------------------------------------------
  d <- d |>
    dplyr::filter(!is.na(.cod7), !is.na(.ano),
                  .ano >= PARAMS$ano_ini, .ano <= PARAMS$ano_fim)
  registrar("SINAN — dentro do período 2018-2025", nrow(d))

  d_est <- dplyr::filter(d, .cod7 %in% cods_estudo)
  registrar("SINAN — residentes na área de estudo", nrow(d_est))

  d_est |>
    dplyr::count(code_muni = .cod7, ano = .ano, name = "casos_sinan")
}

sinan_bruto <- ler_sinan(rota = "arquivo")

# A conversão 6->7 dígitos usa a lista NACIONAL de municípios. Com apenas a
# malha MG+BA, todo residente de outra UF seria contado como "código inválido",
# inflando as exclusões e corrompendo o fluxograma.
tab_muni_br <- tabela_municipios_brasil()
tabela_conversao <- if (is.null(tab_muni_br)) malha else tab_muni_br

casos_sinan <- preparar_sinan(sinan_bruto, cods_estudo, tabela_conversao)


# =============================================================================
# 3. PCE — examinados e positivos
# =============================================================================

pce <- ler_pce() |> dplyr::filter(code_muni %in% cods_estudo)
registrar("PCE — registros município-ano", nrow(pce),
          glue::glue("{sum(pce$examinados, na.rm=TRUE)} exames, ",
                     "{sum(pce$positivos, na.rm=TRUE)} positivos"))

# Consistência: positivos não podem exceder examinados
inconsist <- dplyr::filter(pce, positivos > examinados)
if (nrow(inconsist) > 0) {
  warning(nrow(inconsist), " registros do PCE com positivos > examinados. ",
          "Verifique a exportação; foram marcados como NA.")
  pce <- pce |>
    dplyr::mutate(dplyr::across(c(examinados, positivos),
                                ~ dplyr::if_else(positivos > examinados, NA_real_, .x)))
}


# =============================================================================
# 4. POPULAÇÃO
# =============================================================================

populacao <- obter_populacao() |> dplyr::filter(code_muni %in% cods_estudo)

falta_pop <- expand.grid(code_muni = cods_estudo,
                         ano = PARAMS$ano_ini:PARAMS$ano_fim,
                         stringsAsFactors = FALSE) |>
  dplyr::anti_join(populacao, by = c("code_muni", "ano"))

if (nrow(falta_pop) > 0) {
  message(glue::glue("{nrow(falta_pop)} pares município-ano sem população. ",
                     "Aplicando interpolação/extrapolação log-linear intramunicipal."))
  populacao <- dplyr::bind_rows(populacao, falta_pop |> dplyr::mutate(populacao = NA_real_)) |>
    dplyr::arrange(code_muni, ano) |>
    dplyr::group_by(code_muni) |>
    dplyr::mutate(
      populacao = if (sum(!is.na(populacao)) >= 2) {
        exp(approx(ano[!is.na(populacao)], log(populacao[!is.na(populacao)]),
                   xout = ano, rule = 2)$y)
      } else {
        zoo_fill <- populacao[!is.na(populacao)]
        if (length(zoo_fill)) rep(zoo_fill[1], dplyr::n()) else populacao
      }
    ) |>
    dplyr::ungroup()
}


# =============================================================================
# 5. COVARIÁVEIS ESTRUTURAIS (fixas no tempo)
# =============================================================================

# Bloco socioeconômico: Ipeadata (automático). Um arquivo manual do Atlas em
# dados/brutos/atlas_idhm.csv, se existir, tem precedência.
atlas <- obter_atlas_ipeadata()
idhm  <- if (file.exists(file.path(PARAMS$dir_brutos, "atlas_idhm.csv"))) {
  ler_idhm()
} else if (!is.null(atlas)) {
  atlas
} else {
  tibble::tibble(code_muni = character())
}

ivs       <- ler_ivs()   # IVS municipal só se a planilha do Ipea for fornecida
saneam    <- obter_saneamento_censo()

# Registra os anos-base efetivamente usados — vai para o texto do artigo
registrar("Ano-base do bloco socioeconômico (Atlas DH)",
          attr(idhm, "ano_base") %||% attr(atlas, "ano_base") %||% NA,
          glue::glue("{max(0, ncol(idhm)-1)} indicadores; covariáveis estruturais fixas"))
if (nrow(ivs) > 0) {
  registrar("Ano-base do IVS", attr(ivs, "ano_base") %||% NA,
            "covariável estrutural fixa")
}

estruturais <- malha_est |>
  sf::st_drop_geometry() |>
  dplyr::select(code_muni, name_muni, abbrev_state, regiao, area_km2) |>
  dplyr::left_join(idhm,   by = "code_muni") |>
  dplyr::left_join(ivs,    by = "code_muni") |>
  dplyr::left_join(saneam, by = "code_muni")

# Densidade demográfica com a população do ano central do período
pop_central <- populacao |>
  dplyr::filter(ano == round(mean(c(PARAMS$ano_ini, PARAMS$ano_fim)))) |>
  dplyr::select(code_muni, pop_central = populacao)

estruturais <- estruturais |>
  dplyr::left_join(pop_central, by = "code_muni") |>
  dplyr::mutate(dens_demografica = pop_central / area_km2)

# Diagnóstico de cobertura das covariáveis
cobertura_cov <- estruturais |>
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric),
                                 ~ sum(!is.na(.x)) / dplyr::n() * 100)) |>
  tidyr::pivot_longer(dplyr::everything(), names_to = "variavel",
                      values_to = "pct_preenchido") |>
  dplyr::arrange(pct_preenchido)
print(as.data.frame(cobertura_cov))
salvar_tabela(cobertura_cov, "00_cobertura_covariaveis")


# =============================================================================
# 6. COVARIÁVEIS DE PAINEL (variam no tempo)
# =============================================================================

snis <- ler_snis() |> dplyr::filter(code_muni %in% cods_estudo)

# O SNIS é declaratório: ausência de resposta NÃO é cobertura zero.
# Estratégia: manter NA e, para os modelos, usar (a) LOCF/imputação por
# última observação válida do próprio município quando houver >= 1 resposta,
# (b) NA quando o município nunca respondeu — a ser tratado como missing.
if (nrow(snis) > 0) {
  snis <- snis |>
    dplyr::arrange(code_muni, ano) |>
    dplyr::group_by(code_muni) |>
    dplyr::mutate(dplyr::across(
      c(cob_agua_snis, cob_esgoto_snis),
      ~ if (all(is.na(.x))) .x else approx(ano[!is.na(.x)], .x[!is.na(.x)],
                                           xout = ano, rule = 2)$y,
      .names = "{.col}"
    )) |>
    dplyr::ungroup()

  nunca_respondeu <- snis |>
    dplyr::group_by(code_muni) |>
    dplyr::summarise(sempre_na = all(is.na(cob_agua_snis))) |>
    dplyr::filter(sempre_na)
  registrar("SNIS — municípios sem nenhuma declaração no período",
            nrow(nunca_respondeu),
            "mantidos como NA (não imputar zero)")
}


# =============================================================================
# 7. MONTAGEM DO PAINEL
# =============================================================================

esqueleto <- tidyr::expand_grid(
  code_muni = cods_estudo,
  ano       = PARAMS$ano_ini:PARAMS$ano_fim
)

painel <- esqueleto |>
  dplyr::left_join(populacao,   by = c("code_muni", "ano")) |>
  dplyr::left_join(casos_sinan, by = c("code_muni", "ano")) |>
  dplyr::left_join(pce,         by = c("code_muni", "ano")) |>
  dplyr::left_join(estruturais, by = "code_muni") |>
  dplyr::left_join(snis,        by = c("code_muni", "ano")) |>
  dplyr::mutate(
    # Ausência de linha no SINAN = zero caso notificado (o SINAN só registra
    # ocorrência). Ausência no PCE = NÃO houve inquérito -> permanece NA,
    # pois positividade sem examinados é indefinida, não zero.
    casos_sinan = tidyr::replace_na(casos_sinan, 0),

    # --- desfechos -------------------------------------------------------
    inc_sinan   = casos_sinan / populacao * PARAMS$base_taxa,
    positividade = dplyr::if_else(examinados > 0, positivos / examinados * 100,
                                  NA_real_),
    cob_exame   = dplyr::if_else(populacao > 0, examinados / populacao * 100,
                                 NA_real_),

    # --- estrutura -------------------------------------------------------
    ano_c   = ano - min(ano),           # ano centrado, para termos de tendência
    ano_f   = factor(ano),
    uf      = abbrev_state,
    regiao  = factor(regiao, levels = c("Norte de Minas", "Bahia"))
  ) |>
  dplyr::arrange(code_muni, ano)

chk <- checar_painel(painel)

registrar("Painel final — linhas município-ano", nrow(painel))
registrar("Painel — com dado de PCE", sum(!is.na(painel$positividade)),
          glue::glue("{round(mean(!is.na(painel$positividade))*100,1)}% das células"))


# =============================================================================
# 8. AGREGADO DO PERÍODO (para análise espacial de corte único)
# =============================================================================
# Análises espaciais (Moran, LISA, SaTScan puramente espacial, SAR, GWR)
# ganham estabilidade ao agregar os 8 anos. Somamos numeradores e usamos
# população-ano acumulada (pessoa-ano) como denominador — não a média.

agregado <- painel |>
  dplyr::group_by(code_muni) |>
  dplyr::summarise(
    casos_sinan  = sum(casos_sinan, na.rm = TRUE),
    pessoa_ano   = sum(populacao, na.rm = TRUE),
    examinados   = sum(examinados, na.rm = TRUE),
    positivos    = sum(positivos, na.rm = TRUE),
    anos_com_pce = sum(!is.na(positividade)),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    inc_sinan_bruta = casos_sinan / pessoa_ano * PARAMS$base_taxa,
    positividade    = dplyr::if_else(examinados > 0, positivos / examinados * 100,
                                     NA_real_)
  ) |>
  dplyr::left_join(estruturais, by = "code_muni")

# Covariáveis de painel: média do período
if (nrow(snis) > 0) {
  agregado <- agregado |>
    dplyr::left_join(
      painel |>
        dplyr::group_by(code_muni) |>
        dplyr::summarise(dplyr::across(dplyr::any_of(c("cob_agua_snis", "cob_esgoto_snis")),
                                       ~ mean(.x, na.rm = TRUE)), .groups = "drop"),
      by = "code_muni"
    )
}

# Razão de incidência padronizada (SIR/RME): observado / esperado, com o
# esperado calculado por padronização indireta usando a taxa global do estudo.
taxa_global_sinan <- sum(agregado$casos_sinan) / sum(agregado$pessoa_ano)
prop_global_pce   <- sum(agregado$positivos, na.rm = TRUE) /
                     sum(agregado$examinados, na.rm = TRUE)

agregado <- agregado |>
  dplyr::mutate(
    esperado_sinan = pessoa_ano * taxa_global_sinan,
    rme_sinan      = casos_sinan / esperado_sinan,
    esperado_pce   = examinados * prop_global_pce,
    rme_pce        = dplyr::if_else(examinados > 0, positivos / esperado_pce, NA_real_)
  )

message(glue::glue(
  "\nTaxa global SINAN: {round(taxa_global_sinan*PARAMS$base_taxa, 2)}/100 mil pessoa-ano",
  "\nPositividade global PCE: {round(prop_global_pce*100, 2)}%"
))


# =============================================================================
# 9. SALVAR
# =============================================================================

saveRDS(painel,    file.path(PARAMS$dir_processados, "painel.rds"))
saveRDS(agregado,  file.path(PARAMS$dir_processados, "agregado.rds"))
saveRDS(malha_est, file.path(PARAMS$dir_processados, "malha_est.rds"))
salvar_tabela(FLUXO, "00_fluxograma_dados")

message("\nPainel pronto. Próximo: R/04_descritiva_temporal.R")
