library(readr)

# ============================================================
# ETAPA 1: CARREGAR DADOS CORRETAMENTE
# ============================================================

labels <- scan(
  "C:/Users/muvil/Downloads/left_cell_labels_original.csv",
  what = character()
)

adj <- as.matrix(
  read.table(
    "C:/Users/muvil/Downloads/left_adjacency_original.csv",
    header = FALSE
  )
)

dim(adj)

ids0 <- c(101,104,114,131,139,149)

ids <- ids0 + 1

labels[ids]


W <- adj[ids, ids]

W


ids0 <- c(107,104,114,131,139,149)
ids  <- ids0 + 1

W2 <- adj[ids, ids]

W2
#107,104,114,131,139,149 --> N1 = 131 N2 = 139 N3 = 149 N4 = 107 N5 = 104 N6 = 114