# =============================================================================
# 02_obtencao_dados.R — Aquisição das bases
#
# ESTRATÉGIA:
#   (A) automatizado  -> malha municipal (geobr) e população (SIDRA/IBGE)
#   (B) semiautomático -> SINAN via FTP DATASUS (.dbc)
#   (C) manual        -> PCE/SISPCE, SNIS, Atlas Brasil (IDHM), Ipea (IVS),
#                        Censo 2022 (saneamento)
#
# As fontes em (C) NÃO têm API pública estável. Os arquivos devem ser baixados
# e colocados em dados/brutos/ com os nomes indicados abaixo. Cada leitor
# valida o arquivo e falha com mensagem explícita se ele não existir.
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

library(geobr)
library(sidrar)


# =============================================================================
# (A1) MALHA MUNICIPAL — geobr
# =============================================================================

#' Baixa a malha do geobr usando o ano MAIS RECENTE disponível
#'
#' Percorre PARAMS$anos_malha do mais novo para o mais antigo e retorna o
#' primeiro que o geobr conseguir entregar. Evita fixar um ano que ainda não
#' existe (ou que deixou de existir) na versão instalada do pacote.
#' Tabela de TODOS os municípios do Brasil (código de 7 dígitos)
#'
#' Necessária para converter os códigos de 6 dígitos do SINAN. Usar apenas a
#' malha de MG+BA faria com que todo residente de outra UF fosse classificado
#' como "código inválido", inflando artificialmente as exclusões e produzindo
#' um fluxograma incorreto.
#'
#' Fonte: API de Localidades do IBGE (leve, ~5.570 registros).
tabela_municipios_brasil <- function(
    cache = file.path(PARAMS$dir_processados, "municipios_brasil.rds")) {

  if (file.exists(cache)) return(readRDS(cache))

  message("Baixando lista nacional de municípios (IBGE Localidades)...")
  url <- "https://servicodados.ibge.gov.br/api/v1/localidades/municipios"
  js <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  if (is.null(js)) {
    warning("Falha ao obter a lista nacional; usando apenas a malha do estudo. ",
            "As exclusões por 'município inválido' ficarão superestimadas.")
    return(NULL)
  }
  tab <- tibble::tibble(
    code_muni = as.character(js$id),
    name_muni = js$nome,
    uf        = js$microrregiao$mesorregiao$UF$sigla
  )
  saveRDS(tab, cache)
  message("Lista nacional: ", nrow(tab), " municípios")
  tab
}

geobr_mais_recente <- function(fun, anos = PARAMS$anos_malha, ...) {
  for (a in anos) {
    res <- try(suppressWarnings(fun(year = a, showProgress = FALSE, ...)),
               silent = TRUE)
    if (!inherits(res, "try-error") && !is.null(res) && nrow(res) > 0) {
      attr(res, "ano_usado") <- a
      return(res)
    }
  }
  stop("Nenhum dos anos candidatos está disponível no geobr: ",
       paste(anos, collapse = ", "),
       "\nVerifique a conexão ou rode geobr::list_geobr() para ver os anos válidos.")
}

obter_malha <- function(ufs = PARAMS$ufs,
                        cache = file.path(PARAMS$dir_processados, "malha_mun.rds")) {

  if (file.exists(cache)) {
    m <- readRDS(cache)
    PARAMS$ano_malha <<- attr(m, "ano_malha")
    message("Malha lida do cache: ", cache, " (ano ", attr(m, "ano_malha"), ")")
    return(m)
  }

  message("Baixando malha municipal (geobr) — testando anos ",
          paste(PARAMS$anos_malha, collapse = " > "), " ...")

  # Descobre o ano mais recente que atende TODAS as UFs de uma vez, para não
  # misturar vintages de malha entre MG e BA.
  ano_final <- NULL
  for (a in PARAMS$anos_malha) {
    tentativa <- lapply(ufs, function(uf)
      try(suppressWarnings(
        geobr::read_municipality(code_muni = uf, year = a, showProgress = FALSE)),
        silent = TRUE))
    if (all(!vapply(tentativa, inherits, logical(1), "try-error"))) {
      ano_final <- a
      partes <- tentativa
      break
    }
  }
  if (is.null(ano_final)) {
    stop("Nenhum ano de malha disponível para todas as UFs: ",
         paste(PARAMS$anos_malha, collapse = ", "),
         "\nRode geobr::list_geobr() para ver os anos válidos.")
  }

  malha <- dplyr::bind_rows(partes) |>
    dplyr::mutate(code_muni = as.character(code_muni)) |>
    sf::st_make_valid() |>
    sf::st_transform(4674)   # SIRGAS 2000 geográfico

  attr(malha, "ano_malha") <- ano_final
  # Guarda o ano escolhido para os demais scripts (regiões, cálculos métricos)
  PARAMS$ano_malha <<- ano_final

  saveRDS(malha, cache)
  message("Malha ", ano_final, ": ", nrow(malha), " municípios em ",
          paste(ufs, collapse = "/"))
  malha
}


# =============================================================================
# (A2) DELIMITAÇÃO DO NORTE DE MINAS
# =============================================================================

#' Códigos dos municípios de um recorte regional de MG
#'
#' Estratégias (PARAMS$recorte_mg):
#'   "imed"  — Região Geográfica Intermediária de Montes Claros (IBGE 2017),
#'             a divisão regional VIGENTE  <-- padrão do projeto
#'   "meso"  — Mesorregião "Norte de Minas" (IBGE 1989), descontinuada mas
#'             ainda dominante na literatura de esquistossomose
#'   "lista" — CSV fornecido pela pesquisadora
#'
#' A definição escolhida DEVE ser declarada no artigo — os recortes não
#' coincidem e alteram o N de municípios.
municipios_norte_mg <- function(malha, estrategia = PARAMS$recorte_mg,
                                verboso = TRUE) {

  mg <- malha |> dplyr::filter(abbrev_state == "MG")
  ano <- attr(malha, "ano_malha") %||% PARAMS$anos_malha[1]

  if (estrategia == "lista") {
    arq <- file.path(PARAMS$dir_brutos, "municipios_norte_mg.csv")
    if (!file.exists(arq)) {
      stop(glue::glue("Arquivo não encontrado: {arq}\n",
                      "Deve conter uma coluna 'code_muni' (7 dígitos)."))
    }
    lst <- readr::read_csv2(arq, show_col_types = FALSE) |> janitor::clean_names()
    return(pad_cod_muni(lst$code_muni))
  }

  # As geografias regionais do geobr NÃO acompanham os anos da malha municipal
  # (read_intermediate_region existe só para 2017/2019/2020). Por isso cada uma
  # resolve seu próprio ano mais recente, independentemente do ano da malha.
  tentar_anos <- function(fun, anos, ...) {
    for (a in anos) {
      r <- try(suppressWarnings(fun(year = a, showProgress = FALSE, ...)),
               silent = TRUE)
      if (!inherits(r, "try-error") && !is.null(r) && nrow(r) > 0) {
        attr(r, "ano_usado") <- a
        return(r)
      }
    }
    stop("Geografia indisponível nos anos: ", paste(anos, collapse = ", "))
  }

  poligono <- if (estrategia == "meso") {
    m <- tentar_anos(geobr::read_meso_region,
                     c(2020, 2019, 2018, 2017, 2010), code_meso = "MG")
    if (verboso) message("Recorte: Mesorregião Norte de Minas (IBGE 1989, ",
                         "descontinuada) — geografia ", attr(m, "ano_usado"))
    dplyr::filter(m, grepl("Norte de Minas", name_meso, ignore.case = TRUE))
  } else {
    r <- tentar_anos(geobr::read_intermediate_region,
                     c(2020, 2019, 2017), code_intermediate = "MG")
    if (verboso) message("Recorte: Região Geográfica Intermediária de Montes ",
                         "Claros (IBGE 2017, divisão vigente) — geografia ",
                         attr(r, "ano_usado"))
    dplyr::filter(r, grepl("Montes Claros", name_intermediate, ignore.case = TRUE))
  }

  if (nrow(poligono) == 0) {
    stop("Região não encontrada na malha do ano ", ano,
         ". Verifique se o geobr oferece essa geografia nesse ano.")
  }

  poligono <- sf::st_transform(poligono, sf::st_crs(mg))

  # Junção por ponto interno ao polígono (robusto a diferenças de borda)
  suppressWarnings({ pts <- sf::st_point_on_surface(mg) })
  dentro <- sf::st_within(pts, sf::st_union(poligono), sparse = FALSE)[, 1]

  cods <- mg$code_muni[dentro]
  if (verboso) message("Municípios no recorte: ", length(cods))
  cods
}


#' Compara os dois recortes possíveis do Norte de MG
#'
#' Rode ANTES de fechar a decisão metodológica. Mostra quantos municípios cada
#' divisão inclui, quantos são comuns e — o que importa de verdade — QUAIS
#' entram ou saem. Salva a tabela para o material suplementar.
#'
#' Contexto epidemiológico: a área endêmica de esquistossomose no norte de MG
#' acompanha o vale do São Francisco e as bordas do Jequitinhonha. A divisão de
#' 2017 realocou municípios entre as regiões intermediárias de Montes Claros e
#' Teófilo Otoni; verifique se algum município historicamente endêmico ficou
#' de fora antes de adotar o recorte.
diagnostico_recortes <- function(malha) {

  imed <- municipios_norte_mg(malha, "imed", verboso = FALSE)
  meso <- municipios_norte_mg(malha, "meso", verboso = FALSE)

  info <- sf::st_drop_geometry(malha)[, c("code_muni", "name_muni")]
  rotular <- function(x) info$name_muni[match(x, info$code_muni)]

  so_imed <- setdiff(imed, meso)
  so_meso <- setdiff(meso, imed)
  comuns  <- intersect(imed, meso)

  message(glue::glue(
    "\n--- Comparação dos recortes do Norte de MG ---\n",
    "Reg. Interm. Montes Claros (2017): {length(imed)} municípios\n",
    "Mesorregião Norte de Minas (1989): {length(meso)} municípios\n",
    "Em ambos: {length(comuns)}\n",
    "Só na divisão 2017 (ENTRAM): {length(so_imed)}\n",
    "Só na divisão 1989 (SAEM):   {length(so_meso)}\n"
  ))
  if (length(so_meso)) {
    message("SAEM ao adotar 2017: ", paste(rotular(so_meso), collapse = ", "))
  }
  if (length(so_imed)) {
    message("ENTRAM ao adotar 2017: ", paste(rotular(so_imed), collapse = ", "))
  }

  tab <- dplyr::bind_rows(
    tibble::tibble(code_muni = comuns,  situacao = "Em ambos os recortes"),
    tibble::tibble(code_muni = so_imed, situacao = "Só na Reg. Interm. Montes Claros (2017)"),
    tibble::tibble(code_muni = so_meso, situacao = "Só na Mesorregião Norte de Minas (1989)")
  ) |>
    dplyr::mutate(name_muni = rotular(code_muni)) |>
    dplyr::arrange(situacao, name_muni)

  salvar_tabela(tab, "00_comparacao_recortes_norte_mg")
  invisible(list(imed = imed, meso = meso, tabela = tab))
}


# =============================================================================
# (A3) POPULAÇÃO RESIDENTE — SIDRA/IBGE
# =============================================================================

#' Baixa estimativas populacionais municipais por ano
#'
#' Tabela SIDRA 6579 = "População residente estimada".
#' Anos censitários (2022) podem vir da tabela 4714/9514 — o IBGE revisou a
#' série após o Censo 2022. A função tenta o SIDRA e, em caso de falha
#' (indisponibilidade, mudança de tabela), instrui o uso do CSV manual.
obter_populacao <- function(ufs = PARAMS$ufs,
                            anos = PARAMS$ano_ini:PARAMS$ano_fim,
                            cache = file.path(PARAMS$dir_processados, "populacao.rds")) {

  if (file.exists(cache)) {
    message("População lida do cache: ", cache)
    return(readRDS(cache))
  }

  cod_uf <- c(MG = 31, BA = 29)[ufs]

  # NOTA: o parâmetro se chama `ano_alvo`, não `ano`. Dentro de transmute(), o
  # nome `ano` resolveria para a COLUNA "Ano" devolvida pelo SIDRA (character),
  # e não para o argumento da função — produzindo incompatibilidade de tipo
  # silenciosa ao concatenar com outras fontes.
  puxar_ano <- function(ano_alvo) {
    message("SIDRA 6579 — ano ", ano_alvo)
    Sys.sleep(1)  # cortesia com a API
    tryCatch({
      sidrar::get_sidra(
        x = 6579,
        variable = 9324,
        period = as.character(ano_alvo),
        geo = "City",
        geo.filter = list("State" = unname(cod_uf))
      ) |>
        janitor::clean_names() |>
        dplyr::transmute(
          code_muni = pad_cod_muni(municipio_codigo),
          ano       = as.integer(ano_alvo),
          populacao = as.numeric(valor)
        )
    }, error = function(e) {
      warning(glue::glue("Falha SIDRA em {ano_alvo}: {conditionMessage(e)}"))
      NULL
    })
  }

  pop <- purrr::map_dfr(anos, puxar_ano)

  # A tabela 6579 (estimativas) NÃO cobre anos censitários. Em 2022 a população
  # vem do Censo Demográfico (tabela 4709), que é contagem, não estimativa —
  # fonte preferível quando disponível.
  if (2022 %in% setdiff(anos, unique(pop$ano))) {
    message("SIDRA 4709 — Censo 2022 (contagem)")
    censo <- tryCatch({
      sidrar::get_sidra(
        x = 4709, variable = 93, period = "2022", geo = "City",
        geo.filter = list("State" = unname(cod_uf))
      ) |>
        janitor::clean_names() |>
        dplyr::transmute(code_muni = pad_cod_muni(municipio_codigo),
                         ano = 2022L, populacao = as.numeric(valor))
    }, error = function(e) {
      warning("Falha no Censo 2022 (tabela 4709): ", conditionMessage(e)); NULL
    })
    if (!is.null(censo)) {
      pop <- dplyr::bind_rows(pop, censo)
      message("Censo 2022 incorporado: ", nrow(censo), " municípios")
    }
  }

  anos_faltando <- setdiff(anos, unique(pop$ano))
  if (length(anos_faltando) > 0) {
    arq <- file.path(PARAMS$dir_brutos, "populacao_manual.csv")
    message(glue::glue(
      "\nAnos sem dado no SIDRA: {paste(anos_faltando, collapse=', ')}\n",
      "-> Baixe em https://sidra.ibge.gov.br/tabela/6579 (ou Tabnet/IBGE Cidades)\n",
      "   e salve como {arq} com colunas: code_muni;ano;populacao"
    ))
    if (file.exists(arq)) {
      manual <- readr::read_csv2(arq, show_col_types = FALSE) |>
        janitor::clean_names() |>
        dplyr::transmute(code_muni = pad_cod_muni(code_muni),
                         ano = as.integer(ano),
                         populacao = as.numeric(populacao)) |>
        dplyr::filter(ano %in% anos_faltando)
      pop <- dplyr::bind_rows(pop, manual)
      message("Complementado com ", nrow(manual), " linhas de populacao_manual.csv")
    }
  }

  pop <- dplyr::distinct(pop, code_muni, ano, .keep_all = TRUE)
  saveRDS(pop, cache)
  pop
}


# =============================================================================
# (B) SINAN — ESQUISTOSSOMOSE
# =============================================================================

#' Lê notificações de esquistossomose do SINAN
#'
#' Duas rotas:
#'   rota = "microdatasus" — download automático dos .dbc (requer o pacote;
#'          confira o nome do sistema de informação na documentação da versão
#'          instalada, pois a lista de agravos suportados muda entre releases)
#'   rota = "arquivo"      — .dbc/.csv já baixados manualmente em dados/brutos/
#'
#' FTP DATASUS (rota manual):
#'   ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS/FINAIS/  (ESQUBRAA.dbc)
#'   ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS/PRELIM/  (anos recentes)
#'   AA = últimos 2 dígitos do ano.  Ex.: ESQUBR18.dbc = 2018.
#'
#' ATENÇÃO METODOLÓGICA: use SEMPRE o município de RESIDÊNCIA (ID_MN_RESI),
#' não o de notificação — caso contrário os polos assistenciais (Montes Claros,
#' Feira de Santana, Salvador) concentram artificialmente os casos.
ler_sinan <- function(anos = PARAMS$ano_ini:PARAMS$ano_fim,
                      rota = c("arquivo", "microdatasus"),
                      cache = file.path(PARAMS$dir_processados, "sinan.bruto.rds")) {

  rota <- match.arg(rota)
  if (file.exists(cache)) {
    message("SINAN lido do cache: ", cache)
    return(readRDS(cache))
  }

  if (rota == "microdatasus") {
    if (!requireNamespace("microdatasus", quietly = TRUE)) {
      stop("Instale: remotes::install_github('rfsaldanha/microdatasus')")
    }
    # O identificador do agravo varia entre versões do pacote. Verifique com
    # ?microdatasus::fetch_datasus antes de rodar.
    dados <- purrr::map_dfr(anos, function(a) {
      message("microdatasus — ", a)
      tryCatch(
        microdatasus::fetch_datasus(
          year_start = a, year_end = a,
          information_system = "SINAN-ESQUISTOSSOMOSE"
        ),
        error = function(e) {
          warning(glue::glue("Ano {a} indisponível: {conditionMessage(e)}"))
          NULL
        }
      )
    })
  } else {
    arqs <- list.files(PARAMS$dir_brutos, pattern = "^ESQUBR.*\\.(dbc|DBC|csv|CSV)$",
                       full.names = TRUE)
    if (length(arqs) == 0) {
      stop(glue::glue(
        "Nenhum arquivo SINAN em {PARAMS$dir_brutos}.\n",
        "Baixe ESQUBR18.dbc ... ESQUBR25.dbc do FTP do DATASUS\n",
        "(ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS/) ou exporte\n",
        "CSV do TabNet e salve com o mesmo prefixo."
      ))
    }
    ler_um <- function(f) {
      if (grepl("\\.dbc$", f, ignore.case = TRUE)) {
        if (!requireNamespace("read.dbc", quietly = TRUE)) {
          stop("Instale: remotes::install_github('danicat/read.dbc')")
        }
        read.dbc::read.dbc(f, as.is = TRUE)
      } else {
        readr::read_csv2(f, show_col_types = FALSE, guess_max = 1e5)
      }
    }
    dados <- purrr::map_dfr(arqs, ~ ler_um(.x) |> dplyr::mutate(dplyr::across(dplyr::everything(), as.character)))
  }

  dados <- janitor::clean_names(dados)
  saveRDS(dados, cache)
  message("SINAN bruto: ", nrow(dados), " registros, ", ncol(dados), " colunas")
  dados
}


# =============================================================================
# (C1) PCE / SISPCE — Programa de Controle da Esquistossomose
# =============================================================================

#' Lê dados do PCE (examinados e positivos por município-ano)
#'
#' O SISPCE não tem API. Fontes:
#'   - TabNet/DATASUS: http://tabnet.datasus.gov.br/cgi/deftohtm.exe?pce/cnv/pce*.def
#'   - Secretarias estaduais de saúde (SES-MG / SESAB) — normalmente a única
#'     via para 2023-2025
#'
#' Formato esperado (dados/brutos/pce.csv, separador ";"):
#'   code_muni;ano;examinados;positivos
#'   3143302;2018;1240;87
#'
#' Se as planilhas vierem por faixa etária/localidade, agregue antes ou
#' adapte o bloco de agregação abaixo.
ler_pce <- function(arquivo = file.path(PARAMS$dir_brutos, "pce.csv")) {

  if (!file.exists(arquivo)) {
    warning(glue::glue(
      "PCE NÃO ENCONTRADO ({arquivo}).\n",
      "As análises de POSITIVIDADE serão puladas; as de incidência (SINAN) ",
      "prosseguem normalmente.\n\n",
      "Exporte do TabNet (Programa de Controle da Esquistossomose) com:\n",
      "  Linha    = Município de residência\n",
      "  Coluna   = Ano\n",
      "  Conteúdo = Nº de examinados E Nº de positivos (duas exportações)\n",
      "e consolide no formato: code_muni;ano;examinados;positivos"
    ))
    return(tibble::tibble(code_muni = character(), ano = integer(),
                          examinados = numeric(), positivos = numeric()))
  }

  pce <- readr::read_csv2(arquivo, show_col_types = FALSE,
                          locale = readr::locale(encoding = "Latin1")) |>
    janitor::clean_names()

  obrig <- c("code_muni", "ano", "examinados", "positivos")
  falta <- setdiff(obrig, names(pce))
  if (length(falta)) stop("Colunas ausentes em pce.csv: ", paste(falta, collapse = ", "))

  pce |>
    dplyr::transmute(
      code_muni  = pad_cod_muni(code_muni),
      ano        = as.integer(ano),
      examinados = as.numeric(examinados),
      positivos  = as.numeric(positivos)
    ) |>
    dplyr::filter(!is.na(code_muni), ano >= PARAMS$ano_ini, ano <= PARAMS$ano_fim) |>
    dplyr::group_by(code_muni, ano) |>
    dplyr::summarise(examinados = sum(examinados, na.rm = TRUE),
                     positivos  = sum(positivos,  na.rm = TRUE),
                     .groups = "drop")
}


# =============================================================================
# (C2) INDICADORES SOCIOECONÔMICOS — Atlas Brasil (IDHM) e Ipea (IVS)
# =============================================================================

#' Seleciona o ano MAIS RECENTE presente num arquivo de indicadores
#'
#' Em vez de fixar 2010 (ou qualquer ano), lê o que existe no arquivo e usa o
#' mais novo, relatando a escolha. Assim, quando o Atlas publicar uma série
#' baseada no Censo 2022, basta rebaixar o CSV — nenhum código muda.
#'
#' @param ano_base "max" (mais recente), "min", ou um ano numérico específico
selecionar_ano <- function(dados, ano_base = "max", rotulo = "indicador") {
  if (!"ano" %in% names(dados) || all(is.na(dados$ano))) {
    message(glue::glue("{rotulo}: arquivo sem coluna 'ano' — usando todas as linhas."))
    return(dados)
  }
  anos <- sort(unique(dados$ano[!is.na(dados$ano)]))
  escolhido <- switch(as.character(ano_base),
    "max" = max(anos), "min" = min(anos), as.numeric(ano_base))
  if (!escolhido %in% anos) {
    stop(glue::glue("{rotulo}: ano {escolhido} não existe no arquivo. ",
                    "Disponíveis: {paste(anos, collapse=', ')}"))
  }
  message(glue::glue(
    "{rotulo}: anos no arquivo = {paste(anos, collapse=', ')} -> usando {escolhido}"
  ))
  out <- dplyr::filter(dados, ano == escolhido)
  attr(out, "ano_base") <- escolhido
  out
}

#' Baixa o bloco socioeconômico do Atlas do Desenvolvimento Humano via Ipeadata
#'
#' O Atlas Brasil não tem API, mas o Ipeadata redistribui as séries municipais
#' do Atlas DH e expõe uma API estável, acessível pelo pacote `ipeadatar`.
#'
#' NOTA SOBRE O `ipeadatar`: nas séries municipais a coluna `uname` volta como
#' NA. O nível territorial é identificado pelo COMPRIMENTO do `tcode`
#' (7 dígitos = município, 2 = UF, 1 = Brasil). Filtrar por `uname` descartaria
#' todos os municípios silenciosamente.
#'
#' IVS MUNICIPAL NÃO ESTÁ DISPONÍVEL. A série `AVS_IVS` do Ipeadata só cobre
#' Brasil e Estados. Para vulnerabilidade social usam-se, em substituição, as
#' medidas municipais de pobreza e renda do próprio Atlas (proporção de pobres,
#' extremamente pobres, vulneráveis à pobreza e renda per capita), que são os
#' insumos do próprio IVS. Se o IVS municipal for indispensável, baixe a
#' planilha do Ipea manualmente e salve como dados/brutos/ipea_ivs.csv.
#'
#' @param ano_base ano censitário das séries do Atlas (1991, 2000 ou 2010)
obter_atlas_ipeadata <- function(ufs = PARAMS$ufs, ano_base = 2010,
                                 cache = file.path(PARAMS$dir_processados,
                                                   "atlas_ipeadata.rds")) {

  if (file.exists(cache)) {
    message("Atlas/Ipeadata lido do cache.")
    return(readRDS(cache))
  }
  if (!requireNamespace("ipeadatar", quietly = TRUE)) {
    warning("Pacote 'ipeadatar' ausente — bloco socioeconômico não obtido.")
    return(NULL)
  }

  # code do Ipeadata -> nome da variável no painel
  series <- c(
    ADH_IDHM       = "idhm",
    ADH_IDHM_E     = "idhm_educacao",
    ADH_IDHM_L     = "idhm_longevidade",
    ADH_IDHM_R     = "idhm_renda",
    ADH_GINI       = "gini",
    ADH_PMPOB      = "pct_pobres",
    ADH_PIND       = "pct_extrem_pobres",
    ADH_PPOB       = "pct_vulner_pobreza",
    ADH_RDPC       = "renda_per_capita",
    ADH_T_ANALF15M = "taxa_analfabetismo",
    ADH_T_AGUA     = "pct_agua_encanada_2010"
  )

  cod_uf <- as.character(unname(c(MG = 31, BA = 29)[ufs]))

  puxar <- function(cd, nome) {
    d <- try(ipeadatar::ipeadata(cd, language = "br"), silent = TRUE)
    if (inherits(d, "try-error")) {
      warning("Falha na série ", cd); return(NULL)
    }
    d <- d[nchar(as.character(d$tcode)) == 7, ]          # só municípios
    d$ano <- as.integer(format(d$date, "%Y"))
    anos_disp <- sort(unique(d$ano))
    ano_usar <- if (ano_base %in% anos_disp) ano_base else max(anos_disp)
    d <- d[d$ano == ano_usar &
             substr(as.character(d$tcode), 1, 2) %in% cod_uf, ]
    message(sprintf("  %-16s %s -> %d municípios (anos: %s)",
                    cd, ano_usar, nrow(d), paste(anos_disp, collapse = "/")))
    tibble::tibble(code_muni = pad_cod_muni(d$tcode),
                   !!nome := as.numeric(d$value))
  }

  message("Baixando bloco socioeconômico do Atlas DH (Ipeadata)...")
  partes <- purrr::imap(series, ~ puxar(.y, .x))
  partes <- purrr::compact(partes)
  if (!length(partes)) return(NULL)

  atlas <- purrr::reduce(partes, dplyr::full_join, by = "code_muni")
  attr(atlas, "ano_base") <- ano_base

  saveRDS(atlas, cache)
  message("Atlas: ", nrow(atlas), " municípios, ", ncol(atlas) - 1, " indicadores")
  atlas
}


#' IDHM e componentes (Atlas Brasil — PNUD/Ipea/FJP)
#'
#' Download: http://www.atlasbrasil.org.br/consulta/planilha
#'   Municípios | UF: MG e BA | IDHM, IDHM Renda, IDHM Educação,
#'   IDHM Longevidade, Gini
#'
#' Por padrão usa o ano MAIS RECENTE presente no arquivo. Baixe a série mais
#' nova que o Atlas oferecer no momento da sua extração.
#'
#' LIMITAÇÃO A DECLARAR: mesmo o IDHM mais recente é anterior ao fim do período
#' do estudo (2025) e não varia ano a ano. Ele entra como covariável ESTRUTURAL
#' FIXA — a associação é de contexto socioeconômico, não efeito temporal.
ler_idhm <- function(arquivo = file.path(PARAMS$dir_brutos, "atlas_idhm.csv"),
                     ano_base = "max") {
  if (!file.exists(arquivo)) {
    warning(glue::glue(
      "IDHM NÃO ENCONTRADO ({arquivo}). Prosseguindo sem essa covariável.\n",
      "Baixe em http://www.atlasbrasil.org.br/consulta/planilha e salve como CSV (;)\n",
      "com colunas mínimas: code_muni;ano;idhm;idhm_renda;idhm_educacao;",
      "idhm_longevidade;gini"
    ))
    return(tibble::tibble(code_muni = character()))
  }
  d <- readr::read_csv2(arquivo, show_col_types = FALSE,
                        locale = readr::locale(encoding = "Latin1")) |>
    janitor::clean_names() |>
    dplyr::mutate(code_muni = pad_cod_muni(code_muni))

  d <- selecionar_ano(d, ano_base, rotulo = "IDHM")
  ano_usado <- attr(d, "ano_base")

  out <- d |>
    dplyr::select(code_muni, dplyr::any_of(c("idhm", "idhm_renda", "idhm_educacao",
                                             "idhm_longevidade", "gini")))
  attr(out, "ano_base") <- ano_usado
  out
}

#' IVS — Índice de Vulnerabilidade Social (Ipea)
#' Download: http://ivs.ipea.gov.br/index.php/pt/planilha
ler_ivs <- function(arquivo = file.path(PARAMS$dir_brutos, "ipea_ivs.csv"),
                    ano_base = "max") {
  if (!file.exists(arquivo)) {
    warning(glue::glue("IVS não encontrado ({arquivo}). Prosseguindo sem essa covariável."))
    return(tibble::tibble(code_muni = character(), ivs = numeric()))
  }
  d <- readr::read_csv2(arquivo, show_col_types = FALSE,
                        locale = readr::locale(encoding = "Latin1")) |>
    janitor::clean_names() |>
    dplyr::mutate(code_muni = pad_cod_muni(code_muni))

  d <- selecionar_ano(d, ano_base, rotulo = "IVS")
  ano_usado <- attr(d, "ano_base")

  out <- d |>
    dplyr::select(code_muni, dplyr::any_of(c("ivs", "ivs_infraestrutura_urbana",
                                             "ivs_capital_humano", "ivs_renda_e_trabalho")))
  attr(out, "ano_base") <- ano_usado
  out
}


# =============================================================================
# (C3) SANEAMENTO — Censo 2022 (IBGE) e SNIS
# =============================================================================

#' Saneamento do Censo Demográfico 2022 (corte transversal, base domiciliar)
#'
#' Via SIDRA (tabelas de Características dos Domicílios). Se a tabela mudar,
#' use o CSV manual.
#' Baixa saneamento do Censo 2022 direto da API do SIDRA
#'
#' Tabelas usadas (todas com nível municipal, N6):
#'   6804 · Principal forma de abastecimento de água   (classificação c301)
#'   6805 · Tipo de esgotamento sanitário              (classificação c11558)
#'   6892 · Destino do lixo                            (classificação c67)
#'
#' ATENÇÃO: NÃO usar as tabelas 10341-10345. Elas cobrem apenas setores
#' censitários urbanos com levantamento de entorno — denominador enviesado que
#' exclui a área rural, justamente onde ocorre a transmissão de esquistossomose.
#'
#' Além dos indicadores clássicos de cobertura, extraem-se duas variáveis
#' EPIDEMIOLOGICAMENTE ESPECÍFICAS para esquistossomose, que medem contato e
#' contaminação hídrica diretos — e não cobertura de serviço em geral:
#'   pct_agua_superficial : domicílios abastecidos por rios/açudes/córregos
#'   pct_esgoto_em_corpo_dagua : esgoto lançado em rio, lago, córrego ou mar
obter_saneamento_censo <- function(
    ufs = PARAMS$ufs,
    cache = file.path(PARAMS$dir_processados, "saneamento_censo.rds"),
    arquivo_manual = file.path(PARAMS$dir_brutos, "censo2022_saneamento.csv")) {

  if (file.exists(cache)) {
    message("Saneamento (Censo 2022) lido do cache.")
    return(readRDS(cache))
  }

  # Se a pesquisadora preferir um arquivo próprio, ele tem precedência
  if (file.exists(arquivo_manual)) {
    message("Usando arquivo manual de saneamento: ", arquivo_manual)
    return(readr::read_csv2(arquivo_manual, show_col_types = FALSE,
                            locale = readr::locale(encoding = "Latin1")) |>
             janitor::clean_names() |>
             dplyr::mutate(code_muni = pad_cod_muni(code_muni),
                           dplyr::across(dplyr::starts_with("pct_"), as.numeric)))
  }

  cod_uf <- unname(c(MG = 31, BA = 29)[ufs])

  puxar <- function(tabela, classif, categorias, rotulos) {
    url <- glue::glue(
      "https://apisidra.ibge.gov.br/values/t/{tabela}",
      "/n6/in%20n3%20{paste(cod_uf, collapse=',')}",
      "/v/381/p/2022/c{classif}/{paste(categorias, collapse=',')}"
    )
    message("SIDRA ", tabela, " (Censo 2022) ...")
    txt <- tryCatch(readLines(url, warn = FALSE), error = function(e) NULL)
    if (is.null(txt)) { warning("Falha na tabela ", tabela); return(NULL) }
    js <- jsonlite::fromJSON(paste(txt, collapse = ""))
    cab <- js[1, ]; dados <- js[-1, , drop = FALSE]
    col_cat <- names(cab)[cab == "Município (Código)"]
    col_cls <- names(cab)[grepl("\\(Código\\)$", cab) &
                            !cab %in% c("Município (Código)", "Variável (Código)",
                                        "Ano (Código)", "Nível Territorial (Código)",
                                        "Unidade de Medida (Código)")]
    # Notação do SIDRA: "-" = ZERO; ".." = não se aplica; "..." = indisponível;
    # "X" = omitido por sigilo. Num censo (contagem completa), "-" numa
    # categoria significa nenhum domicílio naquela condição — é zero, não
    # faltante. Tratá-lo como NA descartaria ~100 municípios sem motivo.
    v_bruto <- dados[["V"]]
    valor <- dplyr::case_when(
      v_bruto == "-"                    ~ 0,
      v_bruto %in% c("..", "...", "X")  ~ NA_real_,
      TRUE                              ~ suppressWarnings(as.numeric(v_bruto))
    )

    out <- tibble::tibble(
      code_muni = pad_cod_muni(dados[[col_cat]]),
      cat       = dados[[col_cls[1]]],
      valor     = valor
    ) |>
      dplyr::filter(cat %in% as.character(categorias)) |>
      dplyr::mutate(indicador = rotulos[match(cat, as.character(categorias))]) |>
      dplyr::select(code_muni, indicador, valor) |>
      tidyr::pivot_wider(names_from = indicador, values_from = valor,
                         values_fn = dplyr::first)
    out
  }

  agua <- puxar(6804, 301, c(72053, 31471, 72090),
                c("dom_total_a", "dom_agua_rede", "dom_agua_superficial"))
  esg  <- puxar(6805, 11558, c(46292, 46290, 72114),
                c("dom_total_e", "dom_esgoto_adeq", "dom_esgoto_corpo_dagua"))
  lixo <- puxar(6892, 67, c(10972, 2520),
                c("dom_total_l", "dom_lixo_coletado"))

  if (is.null(agua) || is.null(esg) || is.null(lixo)) {
    stop("Não foi possível obter o saneamento do Censo 2022 pelo SIDRA. ",
         "Baixe manualmente e salve como ", arquivo_manual)
  }

  san <- agua |>
    dplyr::full_join(esg,  by = "code_muni") |>
    dplyr::full_join(lixo, by = "code_muni") |>
    dplyr::mutate(
      pct_agua_rede             = dom_agua_rede          / dom_total_a * 100,
      pct_agua_superficial      = dom_agua_superficial   / dom_total_a * 100,
      pct_esgoto_adeq           = dom_esgoto_adeq        / dom_total_e * 100,
      pct_esgoto_em_corpo_dagua = dom_esgoto_corpo_dagua / dom_total_e * 100,
      pct_lixo_coletado         = dom_lixo_coletado      / dom_total_l * 100,
      dom_total                 = dom_total_a
    ) |>
    dplyr::select(code_muni, dom_total, dplyr::starts_with("pct_"))

  saveRDS(san, cache)
  message("Saneamento Censo 2022: ", nrow(san), " municípios")
  san
}

#' SNIS — série anual de cobertura de água e esgoto
#'
#' Download: https://www.gov.br/cidades/pt-br/acesso-a-informacao/acoes-e-programas/saneamento/snis
#' (Série Histórica > Água e Esgotos > desagregado por município)
#'
#' CUIDADO: o SNIS é declaratório. Municípios pequenos frequentemente não
#' respondem, gerando NA que NÃO deve ser lido como cobertura zero. O script
#' 03 trata isso explicitamente.
ler_snis <- function(arquivo = file.path(PARAMS$dir_brutos, "snis.csv")) {
  if (!file.exists(arquivo)) {
    warning(glue::glue("SNIS não encontrado ({arquivo}). Painel usará só o Censo."))
    return(tibble::tibble(code_muni = character(), ano = integer(),
                          cob_agua_snis = numeric(), cob_esgoto_snis = numeric()))
  }
  readr::read_csv2(arquivo, show_col_types = FALSE,
                   locale = readr::locale(encoding = "Latin1")) |>
    janitor::clean_names() |>
    dplyr::transmute(
      code_muni       = pad_cod_muni(code_muni),
      ano             = as.integer(ano),
      cob_agua_snis   = as.numeric(cob_agua_snis),
      cob_esgoto_snis = as.numeric(cob_esgoto_snis)
    ) |>
    dplyr::filter(ano >= PARAMS$ano_ini, ano <= PARAMS$ano_fim)
}


# =============================================================================
# EXECUÇÃO
# =============================================================================

if (sys.nframe() == 0L) {

  malha <- obter_malha()

  # Compara as duas divisões regionais ANTES de aplicar o recorte, para você
  # ver exatamente o que a escolha custa
  cmp <- try(diagnostico_recortes(malha), silent = TRUE)
  if (inherits(cmp, "try-error")) {
    message("Comparação de recortes indisponível (geobr não retornou uma das ",
            "geografias). Prosseguindo com PARAMS$recorte_mg = '",
            PARAMS$recorte_mg, "'.")
  }

  cods_norte_mg <- municipios_norte_mg(malha)
  populacao <- obter_populacao()

  saveRDS(cods_norte_mg, file.path(PARAMS$dir_processados, "cods_norte_mg.rds"))
  saveRDS(list(recorte = PARAMS$recorte_mg,
               ano_malha = attr(malha, "ano_malha"),
               n_municipios_mg = length(cods_norte_mg)),
          file.path(PARAMS$dir_processados, "metadados_recorte.rds"))

  message("\n--- Bases automáticas prontas ---")
  message("Recorte MG: ", PARAMS$recorte_mg, " | malha: ",
          attr(malha, "ano_malha"), " | ", length(cods_norte_mg), " municípios")
  message("Próximo passo: colocar em dados/brutos/ os arquivos manuais listados ",
          "em dados/brutos/README_dados.md e rodar R/03_preparacao_painel.R")
}
