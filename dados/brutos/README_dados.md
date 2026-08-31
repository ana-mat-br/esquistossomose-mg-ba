# Arquivos que você precisa baixar manualmente

Coloque todos nesta pasta (`dados/brutos/`). Separador **`;`**, decimal **`,`**,
codificação **Latin-1 ou UTF-8**. A coluna `code_muni` é sempre o código IBGE de
**7 dígitos**.

O que é baixado automaticamente pelos scripts: malha municipal (`geobr`) e
estimativas populacionais (SIDRA/IBGE). Todo o resto está abaixo.

---

## 1. `ESQUBR18.dbc` … `ESQUBR25.dbc` — SINAN

**Fonte:** FTP do DATASUS
`ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS/FINAIS/` (anos fechados)
`ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS/PRELIM/` (anos recentes)

`ESQUBR` + 2 dígitos do ano. Baixe os 8 arquivos (2018–2025).

**Alternativa:** exportar CSV do TabNet com Linha = *Município de residência*,
Coluna = *Ano*, Conteúdo = *Casos confirmados*, e salvar com prefixo `ESQUBR`.

> Os anos preliminares (provavelmente 2024–2025) sofrem revisão. Declare no
> artigo a data de extração e que esses anos são preliminares.

---

## 2. `pce.csv` — Programa de Controle da Esquistossomose

**Fonte:** TabNet/DATASUS (PCE) e/ou SES-MG e SESAB. Para 2023–2025 a via mais
provável é solicitação direta às secretarias estaduais.

```
code_muni;ano;examinados;positivos
3143302;2018;1240;87
2910800;2018;3105;412
```

Se a exportação vier por faixa etária ou localidade, agregue por
município-ano antes de salvar — o script soma, mas é mais seguro conferir.

---

## 3. `atlas_idhm.csv` — IDHM e componentes

**Fonte:** <http://www.atlasbrasil.org.br/consulta/planilha>
Municípios · UF: MG e BA · Indicadores: IDHM, IDHM Renda, IDHM Educação,
IDHM Longevidade, Gini.

**Baixe TODOS os anos que o Atlas oferecer** e mantenha a coluna `ano`. O script
usa automaticamente o **mais recente** e informa qual foi na tela e no
fluxograma (`saidas/tabelas/00_fluxograma_dados.csv`).

```
code_muni;ano;idhm;idhm_renda;idhm_educacao;idhm_longevidade;gini
3143302;2010;0,770;0,731;0,724;0,861;0,56
3143302;2021;0,812;0,778;0,791;0,880;0,52
```

Para reproduzir a literatura anterior com o IDHM 2010, force o ano:

```r
idhm <- ler_idhm(ano_base = 2010)
```

> **Atenção:** o IDHM é do PNUD/Ipea/FJP, **não do IBGE**. Mesmo a série mais
> recente é anterior ao fim do período do estudo e não varia ano a ano — entra
> como covariável **estrutural fixa**. A associação é de contexto
> socioeconômico, não efeito temporal. Declare isso como limitação.

---

## 4. `ipea_ivs.csv` — Índice de Vulnerabilidade Social *(opcional)*

**Fonte:** <http://ivs.ipea.gov.br/index.php/pt/planilha>

```
code_muni;ano;ivs;ivs_infraestrutura_urbana;ivs_capital_humano;ivs_renda_e_trabalho
```

Mesma regra do IDHM: inclua a coluna `ano` com todos os anos disponíveis — o
script usa o mais recente e reporta qual. Se o arquivo estiver ausente, o
pipeline roda sem essa covariável (com aviso).

---

## 5. `censo2022_saneamento.csv` — Censo Demográfico 2022

**Fonte:** <https://sidra.ibge.gov.br> · Censo 2022 · Características dos
domicílios. Você precisa de quatro extrações e consolida em um arquivo:

| Coluna | Conteúdo |
|---|---|
| `pct_agua_rede` | % de domicílios com abastecimento por rede geral |
| `pct_esgoto_adeq` | % com esgotamento adequado (rede geral/pluvial ou fossa ligada) |
| `pct_lixo_coletado` | % com lixo coletado |
| `pct_rural` | % da população residente em área rural |

```
code_muni;pct_agua_rede;pct_esgoto_adeq;pct_lixo_coletado;pct_rural
3143302;95,2;71,4;98,1;12,7
```

---

## 6. `snis.csv` — Série histórica SNIS *(opcional, painel)*

**Fonte:** <https://www.gov.br/cidades/pt-br/acesso-a-informacao/acoes-e-programas/saneamento/snis>
Série Histórica → Água e Esgotos → desagregado por município.

```
code_muni;ano;cob_agua_snis;cob_esgoto_snis
3143302;2018;92,1;58,3
```

> O SNIS é **declaratório**. Municípios pequenos frequentemente não respondem.
> Deixe as células em branco (NA) — **nunca preencha com 0**. O script 03
> interpola dentro de cada município que tenha ao menos uma resposta e mantém
> NA nos que nunca responderam.

---

## 7. `populacao_manual.csv` — apenas se o SIDRA falhar

O script 02 tenta baixar sozinho. Só crie este arquivo para os anos que ele
reportar como faltantes.

```
code_muni;ano;populacao
3143302;2025;414000
```

---

## 8. `municipios_norte_mg.csv` — apenas se `PARAMS$recorte_mg = "lista"`

Uma coluna `code_muni`, um município por linha. Use se sua definição de "Norte
de Minas" não coincidir com nenhuma das duas divisões oficiais.

**Antes de decidir**, rode o diagnóstico comparativo — ele lista nominalmente
quais municípios entram e saem em cada divisão e salva a tabela em
`saidas/tabelas/00_comparacao_recortes_norte_mg.csv`:

```r
source("R/02_obtencao_dados.R")
malha <- obter_malha()
diagnostico_recortes(malha)
```

O padrão do projeto é `"imed"` (Região Geográfica Intermediária de Montes
Claros, divisão IBGE 2017 — a vigente). A alternativa `"meso"` é a mesorregião
Norte de Minas de 1989, descontinuada pelo IBGE mas ainda dominante na
literatura de esquistossomose.

---

## Checklist rápido

```r
source("R/00_setup.R")
arqs <- c("pce.csv", "atlas_idhm.csv", "censo2022_saneamento.csv",
          "ipea_ivs.csv", "snis.csv")
data.frame(arquivo = arqs,
           existe = file.exists(file.path(PARAMS$dir_brutos, arqs)))
list.files(PARAMS$dir_brutos, pattern = "^ESQUBR")
```
