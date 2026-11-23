# Reconstituicao do script colado em 23/11/2026; requer fast_cmi.cpp da VM.
ROI <- "LOP_L_col_09_28"
ALFA <- .05
T_OBS <- 40000L
BURN_IN <- 500L
B <- 199L
N_REP <- 100L
N_CORES <- 2L
SEED_GL <- 2026L
SEED_SUR <- 3026L
ROOT <- normalizePath(getwd(),mustWork=TRUE)
DIR_OUT <- file.path(ROOT,"resultados_MC_rapido")
ARQ_W <- file.path(ROOT,"levantamento_regioes","matrizes",ROI,"W_regiao.rds")
ARQ_CPP <- file.path(ROOT,"fast_cmi.cpp")
stopifnot(file.exists(ARQ_W),file.exists(ARQ_CPP),requireNamespace("Rcpp",quietly=TRUE))
dir.create(DIR_OUT,recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(ARQ_CPP,cacheDir=tempdir(),rebuild=FALSE,showOutput=FALSE)
W <- as.matrix(readRDS(ARQ_W)); storage.mode(W) <- "double"
stopifnot(identical(dim(W),c(64L,64L)),sum(W>0)==268L,sum(W)==437,
 all(diag(W)==0),all(W>=0),identical(rownames(W),colnames(W)))
N <- nrow(W)
Wgl <- log1p(W)/log1p(max(W))
simular_gl <- function(W,tamanho,burn_in) {
  n <- nrow(W); total <- tamanho+burn_in
  X <- matrix(0L,total,n); memoria <- matrix(0,n,n)
  for(t in seq_len(total)) {
    X[t,] <- rbinom(n,1L,plogis(rowSums(W*memoria)))
    memoria <- memoria/2+matrix(X[t,],n,n,byrow=TRUE)
    memoria[X[t,]==1L,] <- 0
  }
  X <- X[(burn_in+1L):total,,drop=FALSE]; colnames(X) <- colnames(W); X
}
entropia_rotulos <- function(ids) {
  cont <- as.double(tabulate(match(ids,unique(ids))))
  p <- cont/sum(cont); -sum(p*log(p))
}
codigos_condicionantes <- function(U,S) {
  cod <- rep.int(1L,nrow(U)); k <- 1L
  for(j in S) {
    uj <- match(U[,j],unique(U[,j]))
    chave <- as.double(cod)+(as.double(uj)-1)*as.double(k)
    cod <- match(chave,unique(chave)); k <- max(cod)
  }
  cod
}
preparar_estatistica <- function(U,y,S) {
  id_s <- codigos_condicionantes(U,S); ycod <- match(y,unique(y))
  chave <- as.double(id_s)+(as.double(ycod)-1)*max(id_s)
  id_ys <- match(chave,unique(chave))
  list(id_s=id_s,id_ys=id_ys,base=entropia_rotulos(id_ys)-entropia_rotulos(id_s))
}
criar_variaveis <- function(X,i,prefixos) {
  tt <- nrow(X)
  ultimo <- c(0L,cummax(ifelse(X[seq_len(tt-1L),i]==1L,seq_len(tt-1L),0L)))
  tempos <- which(ultimo>0L)
  if(!length(tempos)) return(NULL)
  last <- ultimo[tempos]
  U <- matrix(0L,length(tempos),ncol(X))
  for(j in seq_len(ncol(X))) U[,j] <- prefixos[[j]][tempos-1L]-prefixos[[j]][last]
  list(U=U,y=X[tempos,i],tempos=as.integer(tempos),ultimo=as.integer(last))
}
selecionar_alvo_rapido <- function(X,i,prefixos,B,alfa) {
  dados <- criar_variaveis(X,i,prefixos)
  if(is.null(dados)) return(integer(0))
  U <- dados$U; y <- dados$y
  todos <- setdiff(seq_len(ncol(X)),i); selecionados <- integer(0)
  tolerados <- as.integer(B-ceiling(B*(1-alfa)))
  repeat {
    candidatos <- setdiff(todos,selecionados)
    if(!length(candidatos)) break
    prep <- preparar_estatistica(U,y,selecionados)
    obs <- vapply(candidatos,function(j)
      cmi_observada_cpp(as.integer(U[,j]),prep$id_s,prep$id_ys,prep$base),numeric(1))
    melhor <- which.max(obs); maior_obs <- obs[melhor]
    if(!is.finite(maior_obs)||maior_obs<=0) break
    excederam <- 0L
    for(b in seq_len(B)) {
      for(j in candidatos) {
        shift <- sample.int(nrow(X)-1L,1L)
        v <- cmi_surrogate_cpp(prefixos[[j]],dados$tempos,dados$ultimo,
          prep$id_s,prep$id_ys,prep$base,shift)
        if(v>=maior_obs) {excederam <- excederam+1L;break}
      }
      if(excederam>tolerados) break
    }
    if(excederam>tolerados) break
    selecionados <- c(selecionados,candidatos[melhor])
  }
  selecionados
}
salvar_seguro <- function(obj,arq) {
  tmp <- paste0(arq,".tmp_",Sys.getpid())
  saveRDS(obj,tmp)
  if(!file.rename(tmp,arq)) stop("Falha ao salvar ",arq)
}
pt <- preparar_estatistica(matrix(c(0L,1L,0L,1L),ncol=1),
 c(0L,1L,0L,1L),integer(0))
stopifnot(abs(cmi_observada_cpp(c(0L,1L,0L,1L),pt$id_s,pt$id_ys,pt$base)-log(2))<1e-8)
Xmini <- c(1L,0L,1L,1L,0L,0L,1L)
pmini <- as.integer(cumsum(Xmini))
for(sh in 1:6) {
  ser <- c(tail(Xmini,sh),head(Xmini,7L-sh))
  tempi <- as.integer(c(3,4,6,7)); ulti <- as.integer(c(1,1,2,2))
  us <- as.integer(cumsum(ser)[tempi-1L]-cumsum(ser)[ulti])
  ptmp <- preparar_estatistica(matrix(us,ncol=1),c(0L,1L,1L,0L),integer(0))
  a <- cmi_observada_cpp(us,ptmp$id_s,ptmp$id_ys,ptmp$base)
  b <- cmi_surrogate_cpp(pmini,tempi,ulti,ptmp$id_s,ptmp$id_ys,ptmp$base,sh)
  stopifnot(isTRUE(all.equal(a,b,tolerance=1e-10)))
}
cl <- parallel::makeCluster(min(N_CORES,N),type="PSOCK")
tryCatch({
  parallel::clusterExport(cl,"ARQ_CPP")
  parallel::clusterEvalQ(cl,{
    Rcpp::sourceCpp(ARQ_CPP,cacheDir=tempdir(),rebuild=TRUE,showOutput=FALSE);NULL
  })
  parallel::clusterExport(cl,c("N","ALFA","B","criar_variaveis","entropia_rotulos",
    "codigos_condicionantes","preparar_estatistica","selecionar_alvo_rapido","salvar_seguro"))
  for(r in seq_len(N_REP)) {
    pasta_r <- file.path(DIR_OUT,sprintf("rep_%03d",r))
    dir.create(pasta_r,recursive=TRUE,showWarnings=FALSE)
    paths <- file.path(pasta_r,sprintf("alvo_%03d.rds",seq_len(N)))
    pendentes <- which(!file.exists(paths))
    if(!length(pendentes)) next
    arq_X <- file.path(pasta_r,"gl_simulado.rds")
    if(file.exists(arq_X)) X <- readRDS(arq_X) else {
      set.seed(as.integer(SEED_GL+10000L*(r-1L)))
      X <- simular_gl(Wgl,T_OBS,BURN_IN)
      salvar_seguro(X,arq_X)
    }
    stopifnot(identical(dim(X),c(T_OBS,N)),identical(colnames(X),colnames(W)))
    prefixos <- lapply(seq_len(N),function(j) as.integer(cumsum(X[,j])))
    parallel::clusterExport(cl,c("X","prefixos","pasta_r","r","SEED_SUR"))
    cat("Replica",r,"| pendentes",length(pendentes),"\n")
    resultados <- parallel::parLapplyLB(cl,pendentes,function(i) {
      set.seed(as.integer(SEED_SUR+10000L*(r-1L)+i))
      inicio <- proc.time()[[3L]]
      pais <- selecionar_alvo_rapido(X,i,prefixos,B,ALFA)
      segundos <- proc.time()[[3L]]-inicio
      salvar_seguro(list(pais=as.integer(pais),segundos=segundos,alvo=as.integer(i),
        replica=as.integer(r)),file.path(pasta_r,sprintf("alvo_%03d.rds",i)))
      c(i,segundos,length(pais))
    })
    rm(X,prefixos);gc(verbose=FALSE)
  }
},finally=parallel::stopCluster(cl))
