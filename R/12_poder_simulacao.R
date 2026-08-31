# =============================================================================
# 12_poder_simulacao.R — Análise de poder e precisão por simulação
#
# Gera os números da seção 7 do RELATORIO_TECNICO.md. Não depende dos dados do
# estudo: é uma análise a priori, condicional às premissas declaradas abaixo.
#
# Premissas:
#   N = 503 municípios (86 Norte de MG + 417 BA), 8 anos
#   populações log-normais calibradas ao perfil brasileiro (mediana ~15 mil,
#     com uma capital em 2,4 milhões)
#   covariável padronizada (efeito = RR por 1 desvio-padrão)
#   teste bilateral a 5% com variância robusta HC0
#
# Saída: saidas/tabelas/42_poder_*.csv
# =============================================================================

source("R/00_setup.R")
source("R/01_funcoes_auxiliares.R")

set.seed(PARAMS$semente)

N_MUN <- 503L
ANOS  <- PARAMS$ano_fim - PARAMS$ano_ini + 1L
NSIM  <- 400L

gera_pop <- function(n) { p <- rlnorm(n, log(15000), 1.15); p[1] <- 2.4e6; round(p) }

#' Erro-padrão sanduíche HC0 calculado diretamente
#' (equivalente a sandwich::vcovHC(f,"HC0"); verificado com concordância
#'  relativa de 4,5e-06. Usado por ser ~10x mais rápido dentro do laço.)
rob <- function(f) {
  X <- model.matrix(f); w <- f$weights
  br <- solve(crossprod(X * sqrt(w)))
  V  <- br %*% crossprod(X * residuals(f, "response")) %*% br
  b  <- unname(coef(f)["x"]); se <- unname(sqrt(V[2, 2]))
  c(b = b, se = se, p = unname(2 * pnorm(abs(b / se), lower.tail = FALSE)))
}

poder_inc <- function(rr, inc_base, theta, nsim=NSIM){
  sig <- 0; larg <- numeric(nsim)
  for(s in 1:nsim){
    pop <- gera_pop(N_MUN); pa <- pop*ANOS; x <- rnorm(N_MUN)
    mu <- pa*(inc_base/1e5)*exp(log(rr)*x)
    y <- rnbinom(N_MUN, mu=mu, size=theta)
    f <- try(glm(y ~ x + offset(log(pa)), family=poisson()), silent=TRUE)
    if(inherits(f,"try-error")) next
    r <- rob(f)
    if(!is.na(r["p"]) && r["p"] < 0.05) sig <- sig+1
    larg[s] <- exp(r["b"]+1.96*r["se"]) - exp(r["b"]-1.96*r["se"])
  }
  c(poder=sig/nsim, ic=median(larg,na.rm=TRUE))
}

poder_pce <- function(rp, pb, nm, rho, nsim=NSIM){
  sig <- 0; larg <- numeric(nsim)
  for(s in 1:nsim){
    ex <- round(rlnorm(nm, log(800), 1.0))+50; x <- rnorm(nm)
    pi_i <- pmin(pb*exp(log(rp)*x), 0.95)
    if(rho > 0){ k <- (1-rho)/rho; pi_i <- rbeta(nm, pi_i*k, (1-pi_i)*k) }
    y <- rbinom(nm, ex, pi_i)
    f <- try(glm(y ~ x + offset(log(ex)), family=poisson()), silent=TRUE)
    if(inherits(f,"try-error")) next
    r <- rob(f); if(!is.na(r["p"]) && r["p"]<0.05) sig <- sig+1
    larg[s] <- exp(r["b"]+1.96*r["se"]) - exp(r["b"]-1.96*r["se"])
  }
  c(poder=sig/nsim, ic=median(larg,na.rm=TRUE))
}

# =============================================================================
# EXECUÇÃO
# =============================================================================

message("A. Poder — incidência SINAN (superdispersão forte, theta=1)")
gA <- expand.grid(rr = c(1.10, 1.15, 1.20, 1.30, 1.50), inc = c(5, 15, 40))
rA <- t(apply(gA, 1, function(g) poder_inc(g["rr"], g["inc"], theta = 1)))
tabA <- cbind(gA, round(rA, 3)) |>
  setNames(c("RR_verdadeiro", "incidencia_base_100mil", "poder", "largura_IC95")) |>
  dplyr::mutate(cenario = "Superdispersão forte (theta=1)")
print(tabA, row.names = FALSE)

message("\nB. Poder — incidência SINAN (superdispersão moderada, theta=5)")
gB <- expand.grid(rr = c(1.10, 1.15, 1.20, 1.30), inc = 15)
rB <- t(apply(gB, 1, function(g) poder_inc(g["rr"], g["inc"], theta = 5)))
tabB <- cbind(gB, round(rB, 3)) |>
  setNames(c("RR_verdadeiro", "incidencia_base_100mil", "poder", "largura_IC95")) |>
  dplyr::mutate(cenario = "Superdispersão moderada (theta=5)")
print(tabB, row.names = FALSE)

salvar_tabela(dplyr::bind_rows(tabA, tabB), "42_poder_incidencia_sinan")

message("\nC. Poder — positividade PCE (com heterogeneidade extra-binomial)")
gC <- expand.grid(rp = c(1.10, 1.20, 1.30, 1.50), nm = c(150, 300),
                  rho = c(0, 0.02, 0.05))
rC <- t(apply(gC, 1, function(z) poder_pce(z["rp"], 0.06, z["nm"], z["rho"])))
tabC <- cbind(gC, round(rC, 3)) |>
  setNames(c("RP_verdadeiro", "n_municipios_PCE", "rho_intraclasse",
             "poder", "largura_IC95"))
print(tabC[order(tabC$rho, tabC$n_municipios_PCE, tabC$RP_verdadeiro), ],
      row.names = FALSE)
salvar_tabela(tabC, "43_poder_positividade_pce")

message("\n>>> rho=0 (binomial pura) é IRREALISTA — produz poder 1,00 em toda a ",
        "grade.\n    Efeito de desenho com n medio ~1300 exames: ",
        "rho=0,02 -> deff ", round(1 + 1299 * 0.02, 1),
        " | rho=0,05 -> deff ", round(1 + 1299 * 0.05, 1))

message("\nD. Instabilidade da taxa bruta (taxa VERDADEIRA = 15/100 mil em todos)")
set.seed(99)
pop <- gera_pop(N_MUN); pa <- pop * ANOS
y   <- rpois(N_MUN, pa * (15 / 1e5)); tx <- y / pa * 1e5
q   <- cut(pop, quantile(pop, 0:4 / 4), include.lowest = TRUE,
           labels = paste0("Q", 1:4))
tabD <- data.frame(
  quartil_populacao  = levels(q),
  populacao_mediana  = round(tapply(pop, q, median)),
  CV_taxa_bruta      = round(tapply(tx, q, function(z) sd(z) / mean(z)), 2),
  taxa_maxima_observada = round(tapply(tx, q, max), 1)
)
print(tabD, row.names = FALSE)
salvar_tabela(tabD, "44_instabilidade_taxa_bruta")

message("\n>>> Toda a variação da tabela D é ruído: a taxa verdadeira é 15/100 mil\n",
        "    em TODOS os municípios. Justifica o estimador bayesiano empírico.")
