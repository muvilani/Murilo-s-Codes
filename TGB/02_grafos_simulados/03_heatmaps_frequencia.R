# Reconstituicao do script de figuras utilizado em 23/11/2026.
# Nao executa GL, mRMR ou surrogates; le checkpoints da bateria v2.
library(ggplot2)
PASTA <- Sys.getenv("TCC_RESULTADOS_GRAFOS", "resultados_4_1_v2_COMPLETA")
SALVAR_PNG <- TRUE
TAMANHOS <- c(1200,6000,12000,60000,120000)
ALFAS <- c(.01,.05,.10)
resultados <- list()
# Os tres cenarios possuem o mesmo suporte, mas pesos distintos.
W <- matrix(c(0,0,1,1,1,0,0,1,1,1,1,1,0,1,-4,
              1,1,1,0,-4,1,1,-4,-4,0),5,byrow=TRUE)
for (cen in 1:3) for (tt in TAMANHOS) {
  arquivos <- list.files(file.path(PASTA,paste0("cenario_",cen),
    paste0("T_",tt)),pattern="^rep_[0-9]+\\.rds$",full.names=TRUE)
  stopifnot(length(arquivos)==100L)
  cont <- array(0,c(5,5,3))
  replicas <- integer()
  for (arq in arquivos) {
    obj <- readRDS(arq)
    stopifnot(is.data.frame(obj$resumo),is.data.frame(obj$selecionados),
      all(obj$resumo$cenario==cen),all(obj$resumo$T==tt),
      length(unique(obj$resumo$replica))==1L)
    replicas <- c(replicas,unique(obj$resumo$replica))
    for (a in seq_along(ALFAS)) {
      stopifnot(sum(obj$resumo$lambda==1 &
        abs(obj$resumo$alfa-ALFAS[a])<1e-8)==1L)
      s <- obj$selecionados
      if (!nrow(s)) next
      s <- unique(s[s$lambda==1 & abs(s$alfa-ALFAS[a])<1e-8,
        c("alvo","candidato"),drop=FALSE])
      stopifnot(all(s$alvo %in% 1:5),all(s$candidato %in% 1:5),
        all(s$alvo!=s$candidato))
      for(k in seq_len(nrow(s))) {
        i <- s$alvo[k]; j <- s$candidato[k]
        cont[i,j,a] <- cont[i,j,a]+1
      }
    }
  }
  stopifnot(setequal(replicas,1:100))
  for(a in seq_along(ALFAS)) {
    df <- expand.grid(alvo=1:5,origem=1:5)
    df$cenario <- cen; df$T <- tt; df$alfa <- ALFAS[a]
    df$frequencia <- as.vector(100*cont[,,a]/length(arquivos))
    df$ausente <- W[cbind(df$alvo,df$origem)]==0 & df$alvo!=df$origem
    df$frequencia[df$alvo==df$origem] <- NA_real_
    resultados[[length(resultados)+1L]] <- df
  }
}
dados <- do.call(rbind,resultados)
dados$T <- factor(dados$T,levels=TAMANHOS,
  labels=paste0("T = ",format(TAMANHOS,big.mark=".",trim=TRUE,scientific=FALSE)))
dados$alfa <- factor(dados$alfa,levels=ALFAS,
  labels=c("alfa = 0.01","alfa = 0.05","alfa = 0.10"))
dados$alvo <- factor(dados$alvo,levels=5:1,labels=paste0("N",5:1))
dados$origem <- factor(dados$origem,levels=1:5,labels=paste0("N",1:5))
graficos <- vector("list",3)
for(cen in 1:3) {
  df <- subset(dados,cenario==cen)
  p <- ggplot(df,aes(origem,alvo))+
    geom_tile(aes(fill=frequencia),color="grey80")+
    geom_tile(data=subset(df,ausente),fill=NA,color="red",linewidth=.8)+
    geom_text(data=subset(df,!is.na(frequencia)),
      aes(label=sprintf("%.0f%%",frequencia),color=frequencia>70),size=2.8)+
    scale_color_manual(values=c("FALSE"="black","TRUE"="white"),guide="none")+
    facet_grid(alfa~T)+
    scale_fill_gradient(low="white",high="#2171B5",limits=c(0,100),
      na.value="grey90",name="Frequencia de\nselecao (%)")+
    coord_fixed()+labs(title=paste("Cenario",cen),
      subtitle="lambda = 1 | 100 replicas por configuracao | Borda vermelha: conexao inexistente",
      x="Neuronio de origem (pre)",y="Neuronio-alvo (pos)")+
    theme_minimal()+theme(panel.grid=element_blank(),
      strip.text=element_text(face="bold"),plot.title=element_text(face="bold",hjust=.5),
      plot.subtitle=element_text(hjust=.5))
  graficos[[cen]] <- p
  if(SALVAR_PNG) ggsave(paste0("cenario_",cen,"_heatmaps.png"),
      plot=p,width=18,height=10,dpi=300,device="png")
  else print(p)
}
