# =============================================================================
# TCC - SECAO 4.1 - TABELAS, DESVIOS-PADRAO E RECUPERACAO POR SINAPSE
# POS-PROCESSAMENTO DOS CHECKPOINTS EXISTENTES - NAO SIMULA NOVOS DADOS
# R base; compativel com o script v2_info_20261019 fornecido em 20/10/2026.
# =============================================================================

PASTA <- "resultados_4_1_v2_COMPLETA"
ALFA_FIG <- 0.05
LAMBDAS_FIG <- c(0, 0.5, 1)

if (!dir.exists(PASTA)) stop("Pasta de resultados nao encontrada: ", PASTA)
SAIDA <- file.path(PASTA, "complemento_ricardo")
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)

# Mesmo W[i,j] do codigo original: j -> i (linha = alvo; coluna = origem).
matriz_cenario <- function(cenario) {
  if (cenario == 1L) return(matrix(c(
    0,0,1,1,1, 0,0,1,1,1, 1,1,0,1,-4,
    1,1,1,0,-4, 1,1,-4,-4,0), nrow=5L, byrow=TRUE))
  if (cenario == 2L) return(matrix(c(
    0,0,3,3,3, 0,0,3,3,3, 3,3,0,3,-12,
    3,3,3,0,-12, 3,3,-12,-12,0), nrow=5L, byrow=TRUE))
  if (cenario == 3L) return(matrix(c(
    0,0,3,3,3, 0,0,1,1,1, 3,1,0,1,-12,
    3,1,1,0,-4, 3,1,-12,-4,0), nrow=5L, byrow=TRUE))
  stop("Cenario inesperado: ", cenario)
}

salvar_csv <- function(x, nome) {
  write.csv(x, file.path(SAIDA, nome), row.names = FALSE, na = "")
}

chave_par <- function(a, b) paste(a, b, sep = "_")

arquivos <- list.files(PASTA, pattern="^rep_[0-9]+[.]rds$",
                       recursive=TRUE, full.names=TRUE)
if (!length(arquivos)) {
  stop("Nenhum rep_XXX.rds encontrado. Confira PASTA e se a bateria v2 ja rodou.")
}

cat("Localizados", length(arquivos), "checkpoints; iniciando leitura...\n")
lista_resumo <- list()
lista_sinapses <- list()
registro <- list()
n_ok <- 0L
n_ignorados <- 0L

for (arquivo in arquivos) {
  obj <- tryCatch(readRDS(arquivo), error=function(e) NULL)
  if (is.null(obj) || !is.list(obj) || !is.data.frame(obj$resumo) ||
      !is.data.frame(obj$selecionados) || is.null(obj$assinatura)) {
    n_ignorados <- n_ignorados + 1L
    next
  }
  d <- obj$resumo
  obrigatorias <- c("cenario", "T", "replica", "lambda", "alfa",
                    "sensibilidade", "especificidade", "FP", "FN",
                    "proporcao_alvos_exatos", "grafo_exato")
  if (!all(obrigatorias %in% names(d)) || !nrow(d) ||
      length(unique(d$cenario)) != 1L || length(unique(d$T)) != 1L ||
      length(unique(d$replica)) != 1L || anyDuplicated(d[, c("lambda", "alfa")])) {
    n_ignorados <- n_ignorados + 1L
    next
  }
  cenario <- as.integer(d$cenario[1L])
  Tobs <- as.integer(d$T[1L])
  replica <- as.integer(d$replica[1L])
  if (!isTRUE(obj$assinatura$cenario == cenario) ||
      !isTRUE(obj$assinatura$T == Tobs) || !(cenario %in% 1:3)) {
    n_ignorados <- n_ignorados + 1L
    next
  }
  W <- matriz_cenario(cenario)
  n <- nrow(W)
  candidatos <- expand.grid(alvo=seq_len(n), origem=seq_len(n),
                            KEEP.OUT.ATTRS=FALSE)
  candidatos <- candidatos[candidatos$alvo != candidatos$origem, , drop=FALSE]
  candidatos$peso_verdadeiro <- W[cbind(candidatos$alvo, candidatos$origem)]
  candidatos$aresta_real <- as.integer(candidatos$peso_verdadeiro != 0)
  candidatos$id <- chave_par(candidatos$alvo, candidatos$origem)
  
  selecionados <- obj$selecionados
  if (nrow(selecionados) &&
      !all(c("alvo", "candidato", "lambda", "alfa") %in% names(selecionados))) {
    n_ignorados <- n_ignorados + 1L
    next
  }
  # Replica tem sempre uma linha para cada combinacao lambda/alpha em obj$resumo.
  lista_resumo[[length(lista_resumo)+1L]] <- d
  linhas_arquivo <- vector("list", nrow(d))
  for (kk in seq_len(nrow(d))) {
    lam <- d$lambda[kk]
    alfa <- d$alfa[kk]
    s <- if (nrow(selecionados)) {
      selecionados[selecionados$lambda == lam & selecionados$alfa == alfa,
                   , drop=FALSE]
    } else selecionados
    id_selecionadas <- if (nrow(s)) unique(chave_par(s$alvo, s$candidato))
    else character(0)
    tab <- candidatos
    tab$cenario <- cenario
    tab$T <- Tobs
    tab$replica <- replica
    tab$lambda <- lam
    tab$alfa <- alfa
    tab$selecionada <- as.integer(tab$id %in% id_selecionadas)
    tab$correta <- as.integer(tab$selecionada == tab$aresta_real)
    tab$id <- NULL
    linhas_arquivo[[kk]] <- tab
    
    # Checagem: contagens sinapse-a-sinapse devem reproduzir TP, FP, FN, TN.
    tp <- sum(tab$selecionada == 1L & tab$aresta_real == 1L)
    fp <- sum(tab$selecionada == 1L & tab$aresta_real == 0L)
    fn <- sum(tab$selecionada == 0L & tab$aresta_real == 1L)
    tn <- sum(tab$selecionada == 0L & tab$aresta_real == 0L)
    if (all(c("TP", "FP", "FN", "TN") %in% names(d))) {
      if (!identical(as.integer(c(tp, fp, fn, tn)),
                     as.integer(unlist(d[kk, c("TP", "FP", "FN", "TN")],
                                       use.names=FALSE)))) {
        stop("Inconsistencia entre selecionados e resumo no checkpoint: ", arquivo,
             "; lambda=", lam, "; alpha=", alfa)
      }
    }
  }
  lista_sinapses[[length(lista_sinapses)+1L]] <- do.call(rbind, linhas_arquivo)
  n_ok <- n_ok + 1L
  if (n_ok %% 25L == 0L) cat("Checkpoints processados:", n_ok, "\n")
}

if (!length(lista_resumo)) stop("Nenhum checkpoint com estrutura valida.")
resumos <- do.call(rbind, lista_resumo)
sinapses <- do.call(rbind, lista_sinapses)
rownames(resumos) <- NULL
rownames(sinapses) <- NULL

# Evita duplicacao silenciosa de checkpoints/repeticoes.
chave_replicacao <- paste(resumos$cenario, resumos$T, resumos$replica,
                          resumos$lambda, resumos$alfa, sep="|")
if (anyDuplicated(chave_replicacao)) {
  stop("Repeticoes duplicadas encontradas entre checkpoints. Revise a pasta.")
}

# 1) Medias e DPs calculados sobre REPLICAS, por configuracao.
metricas <- intersect(c("TP", "FP", "FN", "TN", "sensibilidade",
                        "especificidade", "proporcao_alvos_exatos", "grafo_exato",
                        "relevancia_media", "redundancia_media",
                        "incremental_media", "n_selecionados"), names(resumos))
agrupamento <- list(cenario=resumos$cenario, T=resumos$T,
                    alfa=resumos$alfa, lambda=resumos$lambda)
medias <- aggregate(resumos[metricas], by=agrupamento,
                    FUN=function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm=TRUE))
dps <- aggregate(resumos[metricas], by=agrupamento,
                 FUN=function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm=TRUE))
names(dps)[-(1:4)] <- paste0(names(dps)[-(1:4)], "_dp")
qt <- aggregate(resumos$replica, by=agrupamento,
                FUN=function(x) length(unique(x)))
names(qt)[5] <- "n_replicas"
tabela <- merge(merge(medias, dps, by=c("cenario","T","alfa","lambda")),
                qt, by=c("cenario","T","alfa","lambda"))
tabela <- tabela[order(tabela$cenario,tabela$T,tabela$alfa,tabela$lambda), ]
salvar_csv(tabela, "metricas_media_dp.csv")

# Apresentacao percentual (DP tambem em pontos percentuais).
para_pct <- c("sensibilidade", "especificidade", "proporcao_alvos_exatos", "grafo_exato")
tabela_pct <- tabela
for (nome in para_pct) {
  if (nome %in% names(tabela_pct)) tabela_pct[[nome]] <- 100*tabela_pct[[nome]]
  nome_dp <- paste0(nome, "_dp")
  if (nome_dp %in% names(tabela_pct)) tabela_pct[[nome_dp]] <- 100*tabela_pct[[nome_dp]]
}
salvar_csv(tabela_pct, "metricas_apresentacao_percentual_media_dp.csv")

# 2) Indicadores por elemento de W. Para W[i,j]!=0, recuperacao = selecao;
#    para W[i,j]==0, identificacao correta = NAO selecao.
chaves_sinapse <- c("cenario", "T", "alfa", "lambda", "alvo", "origem",
                    "peso_verdadeiro", "aresta_real")
agrupar_sinapse <- sinapses[chaves_sinapse]
medias_aresta <- aggregate(sinapses[c("selecionada", "correta")],
                           by=agrupar_sinapse, FUN=mean)
dps_aresta <- aggregate(sinapses[c("selecionada", "correta")],
                        by=agrupar_sinapse,
                        FUN=function(x) if (length(x)<2L) NA_real_ else sd(x))
names(dps_aresta)[-(seq_along(chaves_sinapse))] <- c("selecionada_dp", "correta_dp")
n_aresta <- aggregate(sinapses$replica, by=agrupar_sinapse,
                      FUN=function(x) length(unique(x)))
names(n_aresta)[length(names(n_aresta))] <- "n_replicas"
por_sinapse <- merge(merge(medias_aresta, dps_aresta, by=chaves_sinapse),
                     n_aresta, by=chaves_sinapse)
for (cname in c("selecionada", "correta", "selecionada_dp", "correta_dp"))
  por_sinapse[[paste0(cname,"_pct")]] <- 100*por_sinapse[[cname]]
por_sinapse <- por_sinapse[order(por_sinapse$cenario,por_sinapse$T,
                                 por_sinapse$alfa,por_sinapse$lambda,
                                 por_sinapse$alvo,por_sinapse$origem), ]
salvar_csv(por_sinapse, "proporcoes_por_sinapse.csv")

# Matrizes no formato aproximado das Tabelas 5-7 de Ferreira et al. (2025):
# cada celula informa a proporcao de IDENTIFICACAO CORRETA (0 a 100%).
# Linhas: alvo i; colunas: origem j; diagonal excluida (sem autolacos).
combos <- unique(por_sinapse[, c("cenario","T","alfa","lambda")])
linhas_matriz <- list()
for (kk in seq_len(nrow(combos))) {
  ci <- combos[kk, ]
  ss <- por_sinapse[por_sinapse$cenario==ci$cenario & por_sinapse$T==ci$T &
                      por_sinapse$alfa==ci$alfa & por_sinapse$lambda==ci$lambda, ]
  M <- matrix(NA_real_, 5L, 5L)
  M[cbind(ss$alvo,ss$origem)] <- ss$correta_pct
  df <- data.frame(cenario=ci$cenario, T=ci$T, alfa=ci$alfa,
                   lambda=ci$lambda, alvo=seq_len(5L))
  for (j in seq_len(5L)) df[[paste0("origem_",j)]] <- M[,j]
  linhas_matriz[[kk]] <- df
}
salvar_csv(do.call(rbind,linhas_matriz), "matrizes_identificacao_correta_pct.csv")

# 3) Figuras de medias +/- 1 desvio-padrao Monte Carlo, por cenario.
#    O DP nao e um intervalo de confianca da media.
fig_metricas <- c(sensibilidade="Sensibilidade (%)",
                  especificidade="Especificidade (%)",
                  proporcao_alvos_exatos="Alvos exatos (%)",
                  FP="Falsos positivos (contagem)")
cores <- c("0"="#1B4477", "0.5"="#CC6A28", "1"="#267E78")
for (cen in sort(unique(tabela$cenario))) {
  sub <- tabela[tabela$cenario==cen & abs(tabela$alfa-ALFA_FIG)<1e-10 &
                  tabela$lambda %in% LAMBDAS_FIG, ]
  if (!nrow(sub)) next
  for (metrica in names(fig_metricas)) {
    dp_nome <- paste0(metrica, "_dp")
    if (!all(c(metrica, dp_nome) %in% names(sub))) next
    fator <- if (metrica == "FP") 1 else 100
    y <- fator*sub[[metrica]]
    d <- fator*sub[[dp_nome]]
    y0 <- y-ifelse(is.na(d),0,d)
    y1 <- y+ifelse(is.na(d),0,d)
    ylim <- range(c(y0,y1), finite=TRUE)
    if (diff(ylim) < 1e-8) ylim <- ylim+c(-1,1)
    if (metrica != "FP") ylim <- c(max(0,ylim[1L]-4), min(100,ylim[2L]+4))
    else ylim <- c(max(0, ylim[1L]-0.5),ylim[2L]+0.5)
    ts <- sort(unique(sub$T))
    caminho <- file.path(SAIDA, paste0("fig_cenario_",cen,"_",metrica,".png"))
    png(caminho, width=1450, height=950, res=150)
    par(mar=c(5.4,5.5,3.8,1.5), las=1)
    plot(range(ts), ylim, type="n", xaxt="n", xlab="Tamanho da amostra (T)",
         ylab=fig_metricas[[metrica]],
         main=paste0("Cenario ",cen," | alfa = ",ALFA_FIG," | media +/- 1 DP"))
    axis(1,at=ts,labels=format(ts,big.mark=".",scientific=FALSE))
    abline(h=pretty(ylim), col="#EBEBEB", lty=3)
    for (lam in LAMBDAS_FIG) {
      u <- sub[abs(sub$lambda-lam)<1e-10, ]
      u <- u[order(u$T), ]
      if (!nrow(u)) next
      col <- unname(cores[as.character(lam)])
      yy <- fator*u[[metrica]]
      dd <- fator*u[[dp_nome]]
      lines(u$T, yy, type="b", pch=19, lwd=2, col=col)
      ok <- is.finite(dd)
      if (any(ok)) arrows(u$T[ok], pmax(ylim[1L],yy[ok]-dd[ok]),
                          u$T[ok], pmin(ylim[2L],yy[ok]+dd[ok]),
                          angle=90,code=3,length=.055,col=col)
      # O n por ponto aparece acima da serie, quando n<100.
      poucos <- which(u$n_replicas < 100L)
      if (length(poucos)) text(u$T[poucos], yy[poucos],
                               labels=paste0("n=",u$n_replicas[poucos]),
                               pos=3,cex=.75,col=col)
    }
    legend("bottomright",legend=paste0("lambda = ",LAMBDAS_FIG),
           col=unname(cores[as.character(LAMBDAS_FIG)]),
           lty=1,pch=19,bty="n")
    mtext("Barras: DP entre replicas (nao IC). Pontos com n<100 identificados.",
          side=1,line=4.1,cex=.75)
    dev.off()
  }
}

# 4) Mapas de calor de identificacao correta. Usa o MAIOR T que tenha
#    pelo menos 10 replicas completas (evita mapear um unico ensaio como taxa).
for (cen in sort(unique(por_sinapse$cenario))) {
  ds <- por_sinapse[por_sinapse$cenario==cen &
                      abs(por_sinapse$alfa-ALFA_FIG)<1e-10 &
                      por_sinapse$n_replicas>=10L,]
  if (!nrow(ds)) next
  Tref <- max(ds$T)
  for (lam in c(0,1)) {
    s <- ds[ds$T==Tref & abs(ds$lambda-lam)<1e-10,]
    if (nrow(s) != 20L) next
    M <- matrix(NA_real_,5L,5L)
    M[cbind(s$alvo,s$origem)] <- s$correta_pct
    arq <- file.path(SAIDA,paste0("mapa_sinapses_cenario_",cen,
                                  "_T_",Tref,"_lambda_",lam,".png"))
    png(arq,width=1100,height=980,res=150)
    par(mar=c(5,5,4,1))
    plot(NA, xlim=c(.5,5.5),ylim=c(5.5,.5),asp=1,axes=FALSE,
         xlab="Neuronio de origem (j)", ylab="Neuronio-alvo (i)",
         main=paste0("Cenario ",cen," | T=",Tref," | lambda=",lam,
                     " | alfa=",ALFA_FIG," | n=",s$n_replicas[1L]))
    axis(1,at=1:5,labels=1:5)
    axis(2,at=1:5,labels=1:5,las=1)
    pal <- colorRampPalette(c("#f7fbff","#6baed6","#084594"))(101)
    W <- matriz_cenario(cen)
    for (i in 1:5) for (j in 1:5) {
      pct <- M[i,j]
      cor <- if (i==j || !is.finite(pct)) "#DEDEDE" else
        pal[1L+min(100L,max(0L,round(pct)))]
      rect(j-.48,i-.48,j+.48,i+.48,col=cor,border="white",lwd=2)
      if (i!=j && is.finite(pct)) {
        tc <- if (pct>70) "white" else "#0E1B2D"
        text(j,i,sprintf("%.0f%%",pct),col=tc,cex=1.05,font=2)
        text(j,i+.30,if (W[i,j]==0) "ausente" else "presente",
             cex=.65,col=tc)
      } else if (i==j) text(j,i,"-",cex=1.4)
    }
    mtext("W[i,j] : conexao j -> i. Percentual de identificacao CORRETA.",
          side=1,line=3.5,cex=.78)
    dev.off()
  }
}

cat("\nConcluido! Checkpoints lidos:", n_ok,
    "| ignorados:", n_ignorados,"\n")
cat("Resultados em:", normalizePath(SAIDA),"\n")
cat("Arquivos principais:\n",
    "- metricas_apresentacao_percentual_media_dp.csv\n",
    "- proporcoes_por_sinapse.csv\n",
    "- matrizes_identificacao_correta_pct.csv\n",
    "- figuras PNG (medias +/- DP e mapas por sinapse)\n")
