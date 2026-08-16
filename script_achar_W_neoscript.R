## script obter base 
#install.packages("remotes")
#remotes::install_github("natverse/neuprintr")


library(neuprintr)

conn <- neuprint_login(
  server  = "https://neuprint.janelia.org/",
  dataset = "male-cns:v1.0",
  token   = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6Im11dmlsYW5pMkBnbWFpbC5jb20iLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hL0FDZzhvY0lvVUxYb0trVTh5Ti1lX19EbkVMcjdjanRQc1JlQzhqSDU4VWllR1h0WnJER0p4OEZ6PXM5Ni1jP3N6PTUwP3N6PTUwIiwiZXhwIjoxOTYxOTc1OTI0fQ.F4lexdbe2VkyBSHen-44KgKtnnV2rlY0ZFy5G1DJuA4"
)


#neuprint_dataset_info(conn = conn)
mbon_info <- neuprint_search("MBON", field = "type", fixed = TRUE, conn = conn)
print(dim(mbon_info))
print(head(mbon_info))

algum <- neuprint_search(".*", field = "type", fixed = FALSE, conn = conn)

print(sort(unique(mbon_info$type)))


din_info <- neuprint_search("PAM", field = "type", fixed = TRUE, conn = conn)
print(sort(unique(din_info$type)))

ppl_info <- neuprint_search("PPL", field = "type", fixed = TRUE, conn = conn)
print(sort(unique(ppl_info$type)))

library(dplyr)

# um representante por tipo de MBON, PAM e PPL
mbon_rep <- mbon_info %>%
  group_by(type) %>%
  slice_max(order_by = pre + post, n = 1, with_ties = FALSE) %>%
  ungroup()

pam_rep <- din_info %>%
  group_by(type) %>%
  slice_max(order_by = pre + post, n = 1, with_ties = FALSE) %>%
  ungroup()

ppl_rep <- ppl_info %>%
  group_by(type) %>%
  slice_max(order_by = pre + post, n = 1, with_ties = FALSE) %>%
  ungroup()

entrada_rep <- bind_rows(pam_rep, ppl_rep)

print(nrow(mbon_rep))
print(nrow(entrada_rep))

todos_ids <- c(entrada_rep$bodyid, mbon_rep$bodyid)

conn_mat <- neuprint_get_adjacency_matrix(
  inputids  = todos_ids,
  outputids = todos_ids,
  conn = conn
)

dim(conn_mat)

conn_mat <- as.matrix(conn_mat)

n_entrada <- nrow(entrada_rep)  # 27
n_mbon    <- nrow(mbon_rep)     # 37

ids_entrada <- as.character(entrada_rep$bodyid)
ids_mbon    <- as.character(mbon_rep$bodyid)

bloco_entrada_para_mbon <- conn_mat[ids_entrada, ids_mbon]

rownames(bloco_entrada_para_mbon) <- entrada_rep$type
colnames(bloco_entrada_para_mbon) <- mbon_rep$type

dim(bloco_entrada_para_mbon)

# transforma em formato "longo" para ranquear os pares mais fortes
pares <- expand.grid(
  entrada = rownames(bloco_entrada_para_mbon),
  mbon    = colnames(bloco_entrada_para_mbon),
  stringsAsFactors = FALSE
)
pares$peso <- as.vector(bloco_entrada_para_mbon)

pares_top <- pares[order(-pares$peso), ]
head(pares_top, 20)


tipos_entrada <- c("PPL101", "PPL103", "PPL105")
tipos_saida   <- c("MBON11", "MBON32", "MBON13")

circuito <- bind_rows(
  entrada_rep %>% filter(type %in% tipos_entrada),
  mbon_rep    %>% filter(type %in% tipos_saida)
)

print(circuito[, c("bodyid", "type", "name", "pre", "post")])

ids6 <- circuito$bodyid

M_pre_post <- neuprint_get_adjacency_matrix(
  inputids  = ids6,
  outputids = ids6,
  conn = conn
)

M_pre_post <- as.matrix(M_pre_post)
print(M_pre_post)

ids_ordenados <- c(circuito$bodyid[circuito$type == "PPL101"],
                   circuito$bodyid[circuito$type == "PPL103"],
                   circuito$bodyid[circuito$type == "PPL105"],
                   circuito$bodyid[circuito$type == "MBON11"],
                   circuito$bodyid[circuito$type == "MBON32"],
                   circuito$bodyid[circuito$type == "MBON13"])

M_ordenada <- M_pre_post[as.character(ids_ordenados), as.character(ids_ordenados)]

labels <- c("PPL101", "PPL103", "PPL105", "MBON11", "MBON32", "MBON13")
rownames(M_ordenada) <- colnames(M_ordenada) <- labels

print(M_ordenada)

# transpor a matriz porque no site a linha é o pré e a coluna é o pós 
M_invertida <- t(M_ordenada)

print(M_invertida)

write.csv(M_invertida, "C:/Users/muvil/Downloads/matriz_male_cns_6neurons.csv", row.names = TRUE)

