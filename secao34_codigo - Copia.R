# ============================================================
# TCC - Murilo Cassiavilani
# Secao 3.4: Aplicacao em dados reais
# Base: Wisconsin Breast Cancer (pacote mlbench)
#
# Metodo:
#   - Selecao mRMR via pacote praznik
#   - Primeiro estagio: selecao do conjunto candidato S_{n*}
#   - Segundo estagio: refinamento backward e forward (wrapper)
#
# Metrica:
#   - 1 - F1
#   (base levemente desbalanceada: 444 x 239)
#
# Classificador:
#   - Regressao Logistica
# ============================================================

# ============================================================
# INSTALACAO DE PACOTES (se necessario)
# ============================================================

# install.packages(c(
#   "mlbench",
#   "caret",
#   "e1071",
#   "praznik"
# ))

suppressMessages({
  library(mlbench)
  library(caret)
  library(e1071)
  library(praznik)
})

set.seed(42)

# ============================================================
# ETAPA 1: PREPARACAO DOS DADOS
# ============================================================

# Carregar base Wisconsin Breast Cancer
data(BreastCancer)

# Remover identificador
dados <- BreastCancer[, -1]

# Remover observacoes com NA
dados <- na.omit(dados)

# Variavel resposta:
# benigno = 0
# maligno = 1
dados$Class <- ifelse(dados$Class == "malignant", 1, 0)

dados$Class <- factor(
  dados$Class,
  levels = c(0, 1),
  labels = c("benigno", "maligno")
)

# Nomes das preditoras
preds <- setdiff(names(dados), "Class")

# Converter preditoras para inteiro
for (v in preds) {
  dados[[v]] <- as.integer(as.character(dados[[v]]))
}

cat("====================================\n")
cat("DADOS PREPARADOS\n")
cat("====================================\n")
cat("Observacoes :", nrow(dados), "\n")
cat("Preditoras  :", length(preds), "\n")
cat(
  "Resposta    : benigno =",
  sum(dados$Class == "benigno"),
  "| maligno =",
  sum(dados$Class == "maligno"),
  "\n\n"
)

# ============================================================
# ETAPA 2: FUNCOES AUXILIARES
# ============================================================

# ------------------------------------------------------------
# Informacao mutua empirica
# ------------------------------------------------------------

info_mutua <- function(x, y) {
  
  tbl <- table(x, y)
  
  n <- sum(tbl)
  
  px <- rowSums(tbl) / n
  py <- colSums(tbl) / n
  
  mi <- 0
  
  for (i in seq_along(px)) {
    
    for (j in seq_along(py)) {
      
      pxy <- tbl[i, j] / n
      
      if (pxy > 0) {
        
        mi <- mi + pxy * log(
          pxy / (px[i] * py[j])
        )
      }
    }
  }
  
  max(mi, 0)
}

# ------------------------------------------------------------
# Erro de classificacao:
# retorna 1 - F1
# ------------------------------------------------------------

erro_cv_f1 <- function(
    variaveis,
    dados,
    resposta = "Class",
    k = 5
) {
  
  if (length(variaveis) == 0)
    return(1)
  
  set.seed(42)
  
  folds <- createFolds(
    dados[[resposta]],
    k = k,
    list = TRUE
  )
  
  f1s <- numeric(k)
  
  for (i in seq_len(k)) {
    
    treino <- dados[-folds[[i]], ]
    teste  <- dados[ folds[[i]], ]
    
    fm <- as.formula(
      paste(
        resposta,
        "~",
        paste(variaveis, collapse = "+")
      )
    )
    
    tryCatch({
      
      mod <- glm(
        fm,
        data = treino,
        family = binomial()
      )
      
      prob <- predict(
        mod,
        newdata = teste,
        type = "response"
      )
      
      pred <- factor(
        ifelse(prob >= 0.5,
               "maligno",
               "benigno"),
        levels = levels(dados[[resposta]])
      )
      
      cm <- confusionMatrix(
        pred,
        teste[[resposta]],
        positive = "maligno"
      )
      
      f1s[i] <- cm$byClass["F1"]
      
    }, error = function(e) {
      
      f1s[i] <<- 0
      
    })
  }
  
  1 - mean(f1s, na.rm = TRUE)
}

# ============================================================
# ETAPA 3: PRIMEIRO ESTAGIO -- mRMR (praznik)
# ============================================================

cat("====================================\n")
cat("ESTAGIO 1: mRMR VIA PRAZNIK\n")
cat("====================================\n")

M <- length(preds)

# Matriz de preditoras
X_df <- as.data.frame(dados[, preds])

# Resposta numerica binaria
Y_num <- as.integer(dados$Class == "maligno")

# Aplicar mRMR
res_mrmr <- MRMR(
  X = X_df,
  Y = Y_num,
  k = M
)

# Ordem selecionada
ordem_mrmr <- preds[res_mrmr$selection]

# Escores
escores <- res_mrmr$score

cat("Ordem de selecao mRMR:\n\n")

for (i in seq_along(ordem_mrmr)) {
  
  cat(sprintf(
    "%2d. %-20s escore = %.5f\n",
    i,
    ordem_mrmr[i],
    escores[i]
  ))
}

cat("\n")

# ============================================================
# ETAPA 4: ENCONTRAR n* VIA VALIDACAO CRUZADA
# ============================================================

cat("====================================\n")
cat("VALIDACAO CRUZADA DOS SUBCONJUNTOS\n")
cat("====================================\n")

erros_mrmr <- numeric(M)

for (k in seq_len(M)) {
  
  erros_mrmr[k] <- erro_cv_f1(
    ordem_mrmr[1:k],
    dados
  )
  
  cat(sprintf(
    "k = %d | 1-F1 = %.4f\n",
    k,
    erros_mrmr[k]
  ))
}

# Melhor subconjunto S_{n*}
n_star <- which.min(erros_mrmr)

S_cand <- ordem_mrmr[1:n_star]

e_cand <- erros_mrmr[n_star]

cat("\n")
cat("====================================\n")
cat("MELHOR SUBCONJUNTO mRMR\n")
cat("====================================\n")

cat("n* =", n_star, "\n")

cat(
  "S_{n*} = {",
  paste(S_cand, collapse = ", "),
  "}\n"
)

cat(
  sprintf(
    "Erro (1-F1) = %.4f\n\n",
    e_cand
  )
)

# ============================================================
# ETAPA 5: BASELINE -- MAXIMA RELEVANCIA
# ============================================================

cat("====================================\n")
cat("BASELINE: MAXIMA RELEVANCIA\n")
cat("====================================\n")

# I(X_j ; Y)
mi_alvo <- sapply(
  preds,
  function(v)
    info_mutua(dados[[v]], Y_num)
)

# Ordenacao por relevancia
ordem_mr <- names(
  sort(mi_alvo, decreasing = TRUE)
)

# Avaliar subconjuntos
erros_mr <- sapply(
  seq_len(M),
  function(k)
    erro_cv_f1(ordem_mr[1:k], dados)
)

# Melhor subconjunto
n_mr <- which.min(erros_mr)

S_mr <- ordem_mr[1:n_mr]

e_mr <- erros_mr[n_mr]

cat(
  sprintf(
    "Melhor k = %d | 1-F1 = %.4f\n\n",
    n_mr,
    e_mr
  )
)

# ============================================================
# ETAPA 6: WRAPPER BACKWARD
# ============================================================

cat("====================================\n")
cat("WRAPPER BACKWARD\n")
cat("====================================\n")

S_bw  <- S_cand
e_bw  <- e_cand

mudou <- TRUE

while (mudou && length(S_bw) > 1) {
  
  mudou <- FALSE
  
  for (v in S_bw) {
    
    cand <- setdiff(S_bw, v)
    
    ec <- erro_cv_f1(cand, dados)
    
    if (ec <= e_bw) {
      
      cat(sprintf(
        "Remove '%s': %.4f -> %.4f\n",
        v,
        e_bw,
        ec
      ))
      
      S_bw  <- cand
      e_bw  <- ec
      mudou <- TRUE
      
      break
    }
  }
}

cat("\n")

cat(
  "Subconjunto final = {",
  paste(S_bw, collapse = ", "),
  "}\n"
)

cat(
  sprintf(
    "Erro final (1-F1) = %.4f\n\n",
    e_bw
  )
)

# ============================================================
# ETAPA 7: WRAPPER FORWARD
# ============================================================

cat("====================================\n")
cat("WRAPPER FORWARD\n")
cat("====================================\n")

S_fw <- character(0)

e_fw <- 1

for (v in S_cand) {
  
  cand <- c(S_fw, v)
  
  ec <- erro_cv_f1(cand, dados)
  
  if (ec <= e_fw) {
    
    cat(sprintf(
      "Adiciona '%s': %.4f -> %.4f\n",
      v,
      e_fw,
      ec
    ))
    
    S_fw <- cand
    e_fw <- ec
    
  } else {
    
    cat(sprintf(
      "'%s' nao melhora (%.4f > %.4f)\n",
      v,
      ec,
      e_fw
    ))
    
    next
  }
}

cat("\n")

cat(
  "Subconjunto final = {",
  paste(S_fw, collapse = ", "),
  "}\n"
)

cat(
  sprintf(
    "Erro final (1-F1) = %.4f\n\n",
    e_fw
  )
)

# ============================================================
# ETAPA 8: BENCHMARK
# ============================================================

cat("====================================\n")
cat("BENCHMARK: TODAS AS VARIAVEIS\n")
cat("====================================\n")

e_todas <- erro_cv_f1(preds, dados)

cat(
  sprintf(
    "Erro benchmark (1-F1) = %.4f\n\n",
    e_todas
  )
)

# ============================================================
# ETAPA 9: TABELA RESUMO
# ============================================================

cat("====================================\n")
cat("TABELA RESUMO\n")
cat("====================================\n")

tab <- data.frame(
  
  Metodo = c(
    paste0("Todas as variaveis (k = ", M, ")"),
    paste0("Maxima Relevancia  (k = ", n_mr, ")"),
    paste0("mRMR 1o estagio    (k = ", n_star, ")"),
    paste0("mRMR + Backward    (k = ", length(S_bw), ")"),
    paste0("mRMR + Forward     (k = ", length(S_fw), ")")
  ),
  
  k = c(
    M,
    n_mr,
    n_star,
    length(S_bw),
    length(S_fw)
  ),
  
  Erro_1minusF1 = round(
    c(
      e_todas,
      e_mr,
      e_cand,
      e_bw,
      e_fw
    ),
    4
  )
)

print(tab, row.names = FALSE)

# ============================================================
# VARIAVEIS SELECIONADAS
# ============================================================

cat("\n")
cat("====================================\n")
cat("VARIAVEIS SELECIONADAS\n")
cat("====================================\n")

cat(
  "Maxima Relevancia :",
  paste(S_mr, collapse = ", "),
  "\n"
)

cat(
  "mRMR 1o estagio   :",
  paste(S_cand, collapse = ", "),
  "\n"
)

cat(
  "mRMR + Backward   :",
  paste(S_bw, collapse = ", "),
  "\n"
)

cat(
  "mRMR + Forward    :",
  paste(S_fw, collapse = ", "),
  "\n"
)

