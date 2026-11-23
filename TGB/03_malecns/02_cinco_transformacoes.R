# =============================================================================
# TCC | Drosophila Male CNS v1.0 | comparacao de CINCO transformacoes de pesos
# Adapta o experimento 01_experimento_pesos.R; R base, sem pacotes adicionais.
#
# Executar:
#   Rscript comparar_5_transformacoes_drosophila.R piloto 12,25,50
#   Rscript comparar_5_transformacoes_drosophila.R completa 12
#   Rscript comparar_5_transformacoes_drosophila.R completa 25,50
# RStudio: ajuste BASE e FASE, execute source("comparar_5_transformacoes_drosophila.R")
#
# W[i,j] = quantidade de sinapses da origem j para o alvo i.
# W transformada e usada SOMENTE para simular; a verdade e W_original > 0.
# Selecao: mRMR com lambda=1 e parada max-T com B surrogates circulares.
# =============================================================================

VERSAO <- "drosophila_pesos_5_metodos_20261104_v1"
BASE <- "/home/murilo/tcc_malecns/resultados"
FASE <- "piloto"                # piloto: 2 alvos por rede; completa: todos
REDES <- c(12L, 25L, 50L)
T_AMOSTRA <- 10000L
BURN_IN <- 500L
B <- 199L
LAMBDA <- 1
ALFA <- 0.05
SEMENTE_CADEIA <- 2026L
SEMENTE_SURROGATES <- 3026L
METODOS <- c("original", "log1p", "zscore_positivos",
             "divisao_maximo", "divisao_soma")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) FASE <- args[1L]
if (length(args) >= 2L) REDES <- as.integer(strsplit(args[2L], ",", fixed=TRUE)[[1L]])
stopifnot(FASE %in% c("piloto", "completa"), length(REDES)>0L,
          all(REDES %in% c(12L,25L,50L)), !anyDuplicated(REDES),
          B >= 1L, T_AMOSTRA >= 3L, BURN_IN >= 0L)

# Conferir/editar os caminhos, caso os arquivos estejam em outro diretorio.
ARQUIVOS_W <- c(
  "12" = file.path(BASE, "W_anatomica_original.rds"),
  "25" = file.path(BASE, "expansao_anatomica", "W_subrede_N25.rds"),
  "50" = file.path(BASE, "expansao_anatomica", "W_subrede_N50.rds")
)
PASTA_RUN <- file.path(BASE, "comparacao_cinco_pesos_T10000_B199_alfa005")
dir.create(PASTA_RUN, recursive=TRUE, showWarnings=FALSE)

ler_W <- function(N) {
  caminho <- ARQUIVOS_W[as.character(N)]
  if (!file.exists(caminho)) stop("Matriz nao localizada: ", caminho,
                                  "\nAjuste ARQUIVOS_W no inicio do script.")
  W <- as.matrix(readRDS(caminho))
  storage.mode(W) <- "double"
  stopifnot(identical(dim(W), c(as.integer(N), as.integer(N))),
            !is.null(rownames(W)), !is.null(colnames(W)),
            identical(rownames(W),colnames(W)),
            all(is.finite(W)), all(W >= 0), all(diag(W) == 0))
  qtd <- c("12"=67L,"25"=261L,"50"=787L)
  somas <- c("12"=4638,"25"=24416,"50"=48920)
  stopifnot(sum(W>0) == qtd[as.character(N)],
            isTRUE(all.equal(sum(W),unname(somas[as.character(N)]))))
  W
}

# Transformacoes: todas preservam zeros estruturais, exceto efeitos
# semanticos da padronizacao em arestas existentes (valores < 0).
transformar_pesos <- function(W, metodo) {
  pos <- W > 0
  W2 <- W
  if (metodo == "original") {
    # sem mudancas
  } else if (metodo == "log1p") {
    W2[pos] <- log1p(W[pos])   # log(1+w), evita log(0); SEM divisao pelo maximo
  } else if (metodo == "zscore_positivos") {
    media <- mean(W[pos]); desvio <- sd(W[pos])
    if (!is.finite(desvio) || desvio <= 0) stop("DP positivo invalido")
    W2[pos] <- (W[pos] - media)/desvio
  } else if (metodo == "divisao_maximo") {
    W2[pos] <- W[pos]/max(W[pos])
  } else if (metodo == "divisao_soma") {
    W2[pos] <- W[pos]/sum(W[pos])   # NORMALIZACAO GLOBAL
  } else stop("Transformacao desconhecida: ", metodo)
  stopifnot(identical(dimnames(W2),dimnames(W)),
            all(W2[!pos]==0), all(is.finite(W2)))
  W2
}

# Galves-Loecherbach: identico a rotina do experimento anterior.
simular_gl <- function(W, tamanho, burn_in) {
  n <- nrow(W)
  total <- tamanho+burn_in
  X <- matrix(0L,total,n)
  V <- matrix(0,tamanho,n)
  P <- matrix(0,tamanho,n)
  memoria <- matrix(0,n,n)
  for (t in seq_len(total)) {
    potencial <- rowSums(W*memoria)
    probabilidade <- plogis(potencial)
    if (t > burn_in) {
      tt <- t-burn_in
      V[tt,] <- potencial
      P[tt,] <- probabilidade
    }
    X[t,] <- rbinom(n, size=1L, prob=probabilidade)
    memoria <- memoria/2 + matrix(X[t,],nrow=n,ncol=n,byrow=TRUE)
    memoria[X[t,]==1L,] <- 0
  }
  X <- X[(burn_in+1L):total,,drop=FALSE]
  colnames(X) <- colnames(W)
  colnames(P) <- colnames(W)
  colnames(V) <- colnames(W)
  list(X=X,V=V,P=P)
}

# Variaveis U: numero de disparos da origem entre o ultimo spike do alvo e t.
criar_variaveis <- function(X, i) {
  Tobs <- nrow(X)
  ultimo <- c(0L,cummax(ifelse(X[seq_len(Tobs-1L),i]==1L,
                               seq_len(Tobs-1L),0L)))
  tempos <- which(ultimo>0L)
  if (!length(tempos)) return(NULL)
  U <- matrix(0L,length(tempos),ncol(X))
  for (j in seq_len(ncol(X))) {
    cumulativa <- cumsum(X[,j])
    U[,j] <- cumulativa[tempos-1L]-cumulativa[ultimo[tempos]]
  }
  list(U=U,y=X[tempos,i],tempos=tempos,ultimo=ultimo[tempos])
}

entropia <- function(dados) {
  if (is.null(dim(dados))) dados <- matrix(dados,ncol=1L)
  padroes <- do.call(paste,c(as.data.frame(dados),sep="|"))
  freq <- as.numeric(table(padroes))
  p <- freq/sum(freq)
  -sum(p*log(p))
}

informacao <- function(u,y,S) {
  if (ncol(S)==0L) {
    val <- entropia(u)+entropia(y)-entropia(cbind(u,y))
  } else {
    val <- entropia(cbind(u,S))+entropia(cbind(y,S))-
      entropia(S)-entropia(cbind(u,y,S))
  }
  max(0,val)
}

# mRMR + maxima estatistica nula dos candidatos em cada etapa.
# Esta funcao replica a regra de escolha/parada do script existente.
selecionar_alvo <- function(X,i,B,alfa,lambda=1) {
  dados <- criar_variaveis(X,i)
  if (is.null(dados)) return(integer(0))
  U <- dados$U; y <- dados$y
  todos <- setdiff(seq_len(ncol(X)),i)
  vazio <- matrix(nrow=length(y),ncol=0L)
  relevancia <- vapply(todos,function(j) informacao(U[,j],y,vazio),numeric(1))
  names(relevancia) <- as.character(todos)
  selecionados <- integer(0)
  repeat {
    candidatos <- setdiff(todos,selecionados)
    if (!length(candidatos)) break
    S <- U[,selecionados,drop=FALSE]
    observado <- vapply(candidatos,function(j) informacao(U[,j],y,S),numeric(1))
    names(observado) <- as.character(candidatos)
    maximos_nulos <- numeric(B)
    for (b in seq_len(B)) {
      valores_nulos <- numeric(length(candidatos))
      for (k in seq_along(candidatos)) {
        j <- candidatos[k]
        serie <- X[,j]
        deslocamento <- sample.int(nrow(X)-1L,size=1L)
        surrogate <- c(tail(serie,deslocamento),
                       head(serie,nrow(X)-deslocamento))
        acumulada <- cumsum(surrogate)
        u_surrogate <- acumulada[dados$tempos-1L]-acumulada[dados$ultimo]
        valores_nulos[k] <- informacao(u_surrogate,y,S)
      }
      maximos_nulos[b] <- max(valores_nulos)
    }
    limiar <- as.numeric(quantile(maximos_nulos,probs=1-alfa,type=1))
    aprovados <- candidatos[observado > limiar]
    if (!length(aprovados)) break
    notas <- (1-lambda)*relevancia[as.character(aprovados)] +
      lambda*observado[as.character(aprovados)]
    selecionados <- c(selecionados,aprovados[which.max(notas)])
  }
  selecionados
}

# No piloto escolhe os mesmos alvos nos cinco tratamentos da MESMA rede.
# Para W12, tenta manter os dois alvos definidos no experimento anterior.
# Para outras redes, usa um alvo de maior forca de entrada e um de
# menor forca positiva; a selecao e feita SOMENTE sobre W_original.
escolher_alvos <- function(W,FASE,N) {
  if (FASE == "completa") return(seq_len(nrow(W)))
  anterior <- c("MBON11_L","PPL103_R")
  if (N == 12L && all(anterior %in% rownames(W)))
    return(match(anterior,rownames(W)))
  entrada <- rowSums(W)
  candidato_maior <- which.max(entrada)
  cands <- which(entrada>0 & seq_along(entrada)!=candidato_maior)
  if (!length(cands)) cands <- setdiff(seq_len(nrow(W)),candidato_maior)
  candidato_menor <- cands[which.min(entrada[cands])]
  c(candidato_maior,candidato_menor)
}

salvar_seguro <- function(obj,caminho) {
  tmp <- paste0(caminho,".tmp")
  saveRDS(obj,tmp)
  if (.Platform$OS.type == "windows" && file.exists(caminho)) file.remove(caminho)
  if (!file.rename(tmp,caminho)) stop("Nao foi possivel gravar ",caminho)
}

# Verificacoes basicas da implementacao.
z <- c(0,1,0,1)
stopifnot(abs(informacao(z,z,matrix(nrow=4,ncol=0L))-log(2)) < 1e-10)
X_teste <- cbind(c(1,0,0,1,0,0),c(0,1,1,0,1,0))
stopifnot(identical(as.integer(criar_variaveis(X_teste,1)$U[,2]),
                    c(0L,1L,2L,0L,1L)))

TODAS_METRICAS <- list(); TODAS_DINAMICAS <- list(); TODAS_PESOS <- list()
TODAS_ARESTAS <- list(); TODAS_ALVOS <- list()
for (N in REDES) {
  W_original <- ler_W(N)
  n <- nrow(W_original)
  nomes <- rownames(W_original)
  indices_alvos <- escolher_alvos(W_original,FASE,N)
  cat("\n=== REDE N=",N," | alvos: ",paste(nomes[indices_alvos],collapse=", ")," ===\n",sep="")
  pasta_rede <- file.path(PASTA_RUN,paste0("N",N))
  dir.create(pasta_rede,recursive=TRUE,showWarnings=FALSE)
  verdade <- W_original>0
  for (metodo in METODOS) {
    cat("\n--- N=",N," | ",metodo," ---\n",sep="")
    W_sim <- transformar_pesos(W_original,metodo)
    stopifnot(all((W_sim==0)[!verdade]),all(diag(W_sim)==0))
    pasta_metodo <- file.path(pasta_rede,metodo)
    dir.create(pasta_metodo,recursive=TRUE,showWarnings=FALSE)
    config <- list(versao=VERSAO,metodo=metodo,W_original=W_original,
                   W_sim=W_sim,T=T_AMOSTRA,burn_in=BURN_IN,B=B,
                   lambda=LAMBDA,alfa=ALFA,semente_cadeia=SEMENTE_CADEIA,
                   semente_surrogates=SEMENTE_SURROGATES)
    arq_cfg <- file.path(pasta_metodo,"configuracao.rds")
    if (file.exists(arq_cfg)) {
      if (!identical(readRDS(arq_cfg),config))
        stop("Checkpoint incompativel: ",arq_cfg)
    } else salvar_seguro(config,arq_cfg)
    
    arq_cadeia <- file.path(pasta_metodo,"cadeia_dinamica.rds")
    if (file.exists(arq_cadeia)) {
      dados <- readRDS(arq_cadeia)
      cat("Cadeia carregada de checkpoint\n")
    } else {
      set.seed(SEMENTE_CADEIA)
      dados <- simular_gl(W_sim,T_AMOSTRA,BURN_IN)
      salvar_seguro(dados,arq_cadeia)
      cat("Cadeia simulada e salva\n")
    }
    X <- dados$X
    stopifnot(identical(dim(X),c(T_AMOSTRA,n)),
              identical(colnames(X),nomes),all(X %in% 0:1),
              all(is.finite(dados$P)),all(dados$P>=0),all(dados$P<=1))
    pos <- verdade
    pesos <- data.frame(N=N,metodo=metodo,n_arestas=sum(pos),
                        soma_global=sum(W_sim),min_positivo_original=min(W_original[pos]),
                        peso_min_existentes=min(W_sim[pos]),
                        peso_max_existentes=max(W_sim[pos]),
                        peso_medio_existentes=mean(W_sim[pos]),
                        negativos_existentes=sum(W_sim[pos]<0),
                        zerados_existentes=sum(W_sim[pos]==0))
    TODAS_PESOS[[length(TODAS_PESOS)+1L]] <- pesos
    
    indices_gerais <- seq_len(n)
    dinamica <- do.call(rbind,lapply(indices_gerais,function(i)
      data.frame(N=N,metodo=metodo,alvo=nomes[i],
                 analisado=as.integer(i %in% indices_alvos),
                 taxa_disparo=mean(X[,i]),
                 potencial_medio=mean(dados$V[,i]),
                 probabilidade_media=mean(dados$P[,i]),
                 fracao_p_099=mean(dados$P[,i]>0.99),
                 fracao_p_0999=mean(dados$P[,i]>0.999),
                 fracao_p_001=mean(dados$P[,i]<0.01))))
    TODAS_DINAMICAS[[length(TODAS_DINAMICAS)+1L]] <- dinamica
    
    W_est <- matrix(NA_integer_,n,n,dimnames=dimnames(W_original))
    for (i in indices_alvos) {
      arq_alvo <- file.path(pasta_metodo,sprintf("alvo_%03d.rds",i))
      if (file.exists(arq_alvo)) {
        salvo <- readRDS(arq_alvo)
        stopifnot(identical(salvo$config,config),identical(salvo$i,i))
        selecionados <- salvo$selecionados
        cat("Alvo ",nomes[i],": checkpoint (",length(selecionados)," arestas)\n",sep="")
      } else {
        cat("Alvo ",nomes[i],": calculando...\n",sep="")
        set.seed(SEMENTE_SURROGATES+i)
        selecionados <- selecionar_alvo(X,i,B=B,alfa=ALFA,lambda=LAMBDA)
        salvar_seguro(list(config=config,i=i,selecionados=selecionados),arq_alvo)
        cat("Alvo ",nomes[i],": ",length(selecionados)," arestas\n",sep="")
      }
      stopifnot(is.integer(selecionados),!anyDuplicated(selecionados),
                all(selecionados %in% setdiff(seq_len(n),i)))
      W_est[i,] <- 0L
      if(length(selecionados)) W_est[i,selecionados] <- 1L
    }
    idx <- matrix(FALSE,n,n)
    idx[indices_alvos,] <- TRUE
    diag(idx) <- FALSE
    TP <- sum(W_est[idx]==1L & verdade[idx])
    FP <- sum(W_est[idx]==1L & !verdade[idx])
    FN <- sum(W_est[idx]==0L & verdade[idx])
    TN <- sum(W_est[idx]==0L & !verdade[idx])
    stopifnot(TP+FP+FN+TN==length(indices_alvos)*(n-1L))
    metricas <- data.frame(N=N,metodo=metodo,fase=FASE,
                           n_alvos=length(indices_alvos),pares_avaliados=sum(idx),
                           T=T_AMOSTRA,burn_in=BURN_IN,B=B,lambda=LAMBDA,alfa=ALFA,
                           TP=TP,FP=FP,FN=FN,TN=TN,
                           sensibilidade=if(TP+FN) TP/(TP+FN) else NA_real_,
                           especificidade=if(TN+FP) TN/(TN+FP) else NA_real_,
                           precisao=if(TP+FP) TP/(TP+FP) else NA_real_,
                           f1=if(2*TP+FP+FN) 2*TP/(2*TP+FP+FN) else NA_real_,
                           alvos_exatos=sum(vapply(indices_alvos,function(i)
                             identical(which(W_est[i,]==1L),which(verdade[i,])),logical(1))))
    TODAS_METRICAS[[length(TODAS_METRICAS)+1L]] <- metricas
    pares <- which(idx,arr.ind=TRUE)
    tab <- data.frame(N=N,metodo=metodo,
                      alvo=nomes[pares[,1L]],origem=nomes[pares[,2L]],
                      peso_original=W_original[pares],peso_transformado=W_sim[pares],
                      verdadeira=as.integer(verdade[pares]),
                      selecionada=as.integer(W_est[pares]))
    tab$classe <- ifelse(tab$verdadeira==1L & tab$selecionada==1L,"TP",
                         ifelse(tab$verdadeira==0L & tab$selecionada==1L,"FP",
                                ifelse(tab$verdadeira==1L,"FN","TN")))
    TODAS_ARESTAS[[length(TODAS_ARESTAS)+1L]] <- tab
    por_alvo <- do.call(rbind,lapply(indices_alvos,function(i){
      te <- tab[tab$alvo==nomes[i],,drop=FALSE]
      d <- dinamica[dinamica$alvo==nomes[i],,drop=FALSE]
      ta <- sum(te$classe=="TP"); fa <- sum(te$classe=="FP")
      na <- sum(te$classe=="FN"); va <- sum(te$classe=="TN")
      data.frame(N=N,metodo=metodo,alvo=nomes[i],
                 TP=ta,FP=fa,FN=na,TN=va,
                 sensibilidade=if(ta+na)ta/(ta+na) else NA_real_,
                 especificidade=if(fa+va)va/(fa+va) else NA_real_,
                 taxa_disparo=d$taxa_disparo,fracao_p_099=d$fracao_p_099,
                 recuperacao_exata=as.integer(fa==0L && na==0L))
    }))
    TODAS_ALVOS[[length(TODAS_ALVOS)+1L]] <- por_alvo
    write.csv(tab,file.path(pasta_metodo,paste0("arestas_",FASE,".csv")),row.names=FALSE)
    cat("N=",N," | ",metodo," | TP=",TP," FP=",FP," FN=",FN," TN=",TN,"\n",sep="")
  }
}

combinar <- function(L) {res<-do.call(rbind,L);rownames(res)<-NULL;res}
metricas <- combinar(TODAS_METRICAS)
dinamica <- combinar(TODAS_DINAMICAS)
pesos <- combinar(TODAS_PESOS)
arestas <- combinar(TODAS_ARESTAS)
por_alvo <- combinar(TODAS_ALVOS)
for (N in REDES) {
  base <- metricas[metricas$N==N & metricas$metodo=="original",,drop=FALSE]
  stopifnot(nrow(base)==1L)
  indices <- which(metricas$N==N)
  metricas$delta_TP_vs_original[indices] <- metricas$TP[indices]-base$TP
  metricas$delta_FP_vs_original[indices] <- metricas$FP[indices]-base$FP
  metricas$delta_FN_vs_original[indices] <- metricas$FN[indices]-base$FN
}
originais <- arestas[arestas$metodo=="original",c("N","alvo","origem","selecionada")]
names(originais)[4] <- "selecionada_original"
mudancas <- merge(arestas[arestas$metodo!="original",],originais,
                  by=c("N","alvo","origem"),all.x=TRUE,sort=FALSE)
mudancas$mudou <- as.integer(mudancas$selecionada!=mudancas$selecionada_original)
mudancas$ganhou_TP <- as.integer(mudancas$verdadeira==1L &
                                   mudancas$selecionada_original==0L & mudancas$selecionada==1L)
mudancas$perdeu_TP <- as.integer(mudancas$verdadeira==1L &
                                   mudancas$selecionada_original==1L & mudancas$selecionada==0L)

saida <- file.path(PASTA_RUN,"comparacoes",FASE)
dir.create(saida,recursive=TRUE,showWarnings=FALSE)
write.csv(metricas,file.path(saida,"comparacao_metricas.csv"),row.names=FALSE)
write.csv(dinamica,file.path(saida,"comparacao_saturacao.csv"),row.names=FALSE)
write.csv(pesos,file.path(saida,"comparacao_pesos.csv"),row.names=FALSE)
write.csv(arestas,file.path(saida,"comparacao_arestas.csv"),row.names=FALSE)
write.csv(por_alvo,file.path(saida,"comparacao_por_alvo.csv"),row.names=FALSE)
write.csv(mudancas,file.path(saida,"mudancas_vs_original.csv"),row.names=FALSE)

pdf(file.path(saida,"graficos_comparativos.pdf"),width=12,height=8)
for(N in REDES) {
  m <- metricas[metricas$N==N,,drop=FALSE]
  d <- dinamica[dinamica$N==N,,drop=FALSE]
  sat <- aggregate(fracao_p_099~metodo,data=d,FUN=mean)
  sat <- sat[match(m$metodo,sat$metodo),,drop=FALSE]
  par(mfrow=c(2,1),mar=c(6,4,3,1))
  barplot(rbind(m$TP,m$FP,m$FN),beside=TRUE,names.arg=m$metodo,las=2,
          col=c("steelblue","firebrick","gray65"),
          main=paste0("N=",N," | acertos e erros (fase: ",FASE,")"),ylab="Arestas")
  legend("topright",legend=c("TP","FP","FN"),
         fill=c("steelblue","firebrick","gray65"),bty="n")
  barplot(sat$fracao_p_099,names.arg=sat$metodo,las=2,ylim=c(0,1),
          col="steelblue",ylab="Fracao p > 0.99 (media entre todos os alvos)",
          main=paste0("N=",N," | diagnostico de saturacao"))
}
dev.off()
cat("\n=== METRICAS DA COMPARACAO ===\n")
print(metricas,row.names=FALSE)
cat("\nSaidas em: ",saida,"\n",sep="")
cat("OBS: uma realizacao por matriz/metodo. As diferencas nao sao ICs nem testes estatisticos.\n")
