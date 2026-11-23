# Inventario MaleCNS v1.0, reconstituido do codigo colado em 23/11/2026.
# Usa NEUPRINT_TOKEN do ambiente. Nao extrai as matrizes W.
library(httr)
library(jsonlite)
servidor <- "https://neuprint.janelia.org"
dataset <- "male-cns:v1.0"
somente_foco <- FALSE
B <- 199L
dir.create("levantamento_regioes",recursive=TRUE,showWarnings=FALSE)
token <- Sys.getenv("NEUPRINT_TOKEN")
consultar <- function(cypher) {
  cabecalhos <- if(nzchar(token)) add_headers(Authorization=paste("Bearer",token)) else add_headers()
  resp <- POST(paste0(servidor,"/api/custom/custom"),cabecalhos,
    content_type_json(),body=list(cypher=cypher,dataset=dataset),
    encode="json",timeout(120))
  stop_for_status(resp)
  fromJSON(content(resp,as="text",encoding="UTF-8"),simplifyVector=FALSE)
}
meta <- consultar("MATCH (m:Meta) RETURN m.roiInfo AS roiInfo, m.primaryRois AS primaryRois LIMIT 1")
roi_info <- meta$data[[1]][[1]]
if(is.character(roi_info)) roi_info <- fromJSON(roi_info,simplifyVector=FALSE)
rois <- sort(unique(c(names(roi_info),unlist(meta$data[[1]][[2]],use.names=FALSE))))
if(somente_foco) rois <- rois[grepl("^(BU|CAN|AME)(\\(|_|$)",rois,ignore.case=TRUE)]
stopifnot(length(rois)>0)
resultados <- data.frame()
for(roi in rois) {
  acento <- intToUtf8(96)
  roi_segura <- gsub(acento,paste0(acento,acento),roi,fixed=TRUE)
  query <- paste0("MATCH (n:Neuron) WHERE n.",acento,roi_segura,acento," = true RETURN count(n) AS N")
  resposta <- tryCatch(consultar(query),error=function(e){message(conditionMessage(e));NULL})
  N <- if(is.null(resposta)||!length(resposta$data)) NA_real_ else as.numeric(resposta$data[[1]][[1]])
  resultados <- rbind(resultados,data.frame(roi=roi,N=N,pares_dirigidos=N*(N-1),
    avaliacoes_mrmr=(N^2*(N-1)/2)*(B+1)))
  write.csv(resultados,"levantamento_regioes/inventario_parcial.csv",row.names=FALSE)
}
resultados <- resultados[order(resultados$N,na.last=TRUE),]
write.csv(resultados,"levantamento_regioes/inventario_regioes.csv",row.names=FALSE)
print(head(resultados,30),row.names=FALSE)
