# Controlador reconstituido do codigo colado em 23/11/2026.
# Nao altera o script original; executa uma copia parametrizada por T/replica.
ETAPA <- "piloto"
TAMANHOS <- c(10000L,40000L,80000L)
REPLICAS <- 1:5
ROI <- "LOP_L_col_09_28"
if(ETAPA=="piloto") { TAMANHOS <- 10000L; REPLICAS <- 1L } else if(ETAPA!="completo") stop("Modo invalido")
projeto <- normalizePath(getwd(),mustWork=TRUE)
original <- Sys.getenv("TCC_SCRIPT_ROI","TGB/04_roi/03_mrmr_regiao_escolhida.R")
original <- normalizePath(original,mustWork=TRUE)
W <- readRDS(file.path(projeto,"levantamento_regioes","matrizes",ROI,"W_regiao.rds"))
stopifnot(identical(dim(W),c(64L,64L)),sum(W>0)==268L,sum(W)==437,all(diag(W)==0))
raiz <- file.path(projeto,"validacao_profunda_LOP_L_col_09_28")
dir.create(raiz,recursive=TRUE,showWarnings=FALSE)
substituir <- function(linhas,nome,expr) {
  idx <- grep(paste0("^[[:space:]]*",nome,"[[:space:]]*<-[[:space:]]*"),linhas)
  stopifnot(length(idx)==1L)
  linhas[idx] <- paste0(nome," <- ",expr)
  linhas
}
for(tt in TAMANHOS) for(rep in REPLICAS) {
  pasta <- file.path(raiz,paste0("T",tt),sprintf("rep_%03d",rep))
  dir.create(pasta,recursive=TRUE,showWarnings=FALSE)
  link <- file.path(pasta,"levantamento_regioes")
  if(!file.exists(link)) stopifnot(file.symlink(file.path(projeto,"levantamento_regioes"),link))
  valores <- list(REGIOES=sprintf('c("%s")',ROI),METODOS='c("log_maximo")',
    ALFAS="c(0.05,0.10)",LAMBDA="1",T_OBS=paste0(tt,"L"),
    BURN_IN="500L",B="199L",MODO=sprintf('"%s"',ETAPA),
    SEED_GL=paste0(2026L+rep*1000L,"L"),SEED_SUR=paste0(3026L+rep*1000L,"L"))
  linhas <- readLines(original,warn=FALSE)
  for(nm in names(valores)) linhas <- substituir(linhas,nm,valores[[nm]])
  copia <- file.path(pasta,"03_mrmr_regiao_escolhida_configurado.R")
  writeLines(linhas,copia,useBytes=TRUE)
  setwd(pasta)
  erro <- tryCatch({sys.source(copia,envir=new.env(parent=globalenv()));NULL},
    error=function(e) conditionMessage(e),finally=setwd(projeto))
  if(!is.null(erro)) stop(erro)
}
