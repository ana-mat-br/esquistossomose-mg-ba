# =============================================================================
# _run_all.R — Executa o pipeline completo
#
#   Rscript _run_all.R            # tudo
#   Rscript _run_all.R 5 6 7      # apenas os scripts 05, 06 e 07
#
# Os scripts 02 e 03 exigem os arquivos manuais descritos em
# dados/brutos/README_dados.md. Os demais leem apenas dados/processados/.
# =============================================================================

# Execute a partir da RAIZ do projeto (onde está este arquivo).
if (!dir.exists("R")) {
  stop("Rode a partir da raiz do projeto: cd ~/Esquistossomose && Rscript _run_all.R")
}

suppressPackageStartupMessages(library(tibble))

ETAPAS <- tibble::tribble(
  ~n,  ~script,                      ~descricao,                                  ~requer_dados_manuais,
  2,   "R/02_obtencao_dados.R",      "Download de malha e população",             FALSE,
  3,   "R/03_preparacao_painel.R",   "Limpeza e montagem do painel",              TRUE,
  4,   "R/04_descritiva_temporal.R", "Descritiva e tendência (Prais-Winsten)",    FALSE,
  5,   "R/05_bayes_empirico.R",      "Estimadores bayesianos empíricos",          FALSE,
  6,   "R/06_moran_lisa.R",          "Moran global, LISA, Gi*",                   FALSE,
  7,   "R/07_satscan.R",             "Varredura espacial (Kulldorff/SaTScan)",    FALSE,
  8,   "R/08_poisson_robusta.R",     "Poisson com variância robusta",             FALSE,
  9,   "R/09_espacial_sar_gwr.R",    "SAR/SEM/SDM, ESF e GWR",                    FALSE,
  10,  "R/10_bym2_inla.R",           "BYM2 espaço-temporal (requer INLA)",        FALSE,
  11,  "R/11_mapas_tabelas.R",       "Mapas e tabelas finais",                    FALSE,
  12,  "R/12_poder_simulacao.R",     "Poder e precisão (independe dos dados)",    FALSE
)

args <- commandArgs(trailingOnly = TRUE)
alvo <- if (length(args)) as.integer(args) else ETAPAS$n

log_arq <- file.path("saidas", paste0("log_execucao_",
                                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
dir.create("saidas", showWarnings = FALSE)

resultado <- data.frame()

for (i in which(ETAPAS$n %in% alvo)) {
  et <- ETAPAS[i, ]
  cat("\n", strrep("=", 78), "\n",
      sprintf("[%02d] %s — %s\n", et$n, et$script, et$descricao),
      strrep("=", 78), "\n", sep = "")

  t0 <- Sys.time()
  ok <- tryCatch({
    source(et$script, echo = FALSE, local = new.env())
    TRUE
  }, error = function(e) {
    message("\n### ERRO em ", et$script, ":\n", conditionMessage(e), "\n")
    FALSE
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  resultado <- rbind(resultado, data.frame(
    etapa = et$n, script = et$script, sucesso = ok, segundos = round(dt, 1)
  ))

  if (!ok && et$requer_dados_manuais) {
    message("Esta etapa depende de arquivos manuais. ",
            "Consulte dados/brutos/README_dados.md e reexecute.")
    break
  }
}

cat("\n", strrep("=", 78), "\nRESUMO\n", strrep("=", 78), "\n", sep = "")
print(resultado)
utils::write.csv(resultado, log_arq, row.names = FALSE)
cat("\nLog salvo em ", log_arq, "\n", sep = "")

if (all(resultado$sucesso)) {
  cat("\nPipeline concluído com sucesso.\n")
} else {
  cat("\nHouve falhas — ver acima.\n")
  quit(status = 1)
}
