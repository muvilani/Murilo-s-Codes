# ============================================================
# TCC - Murilo Cassiavilani
# Secao 3.4: Aplicacao em dados reais
#
# Base Colon (Alon et al., 1999)
#
# Esta base contem perfis de expressao genica de pacientes
# com cancer colorretal.
#
# Cada variável preditora representa o nível de expressão de um gene,
# medido em unidades de intensidade de fluorescência,
# sendo valores mais altos indicativos de maior atividade
# transcricional daquele gene na amostra.
#
# Estrutura:
#   - 62 pacientes
#   - 2000 genes (preditoras)
#   - resposta binaria:
#         normal  vs  tumor
#
# Trata-se de um problema de alta dimensionalidade:
#
#         p >> n
#
# pois existem muito mais variaveis (2000 genes) do que
# observacoes (62 pacientes).
#
# Nessa situacao, a regressao logistica tradicional torna-se
# instavel ou ate impossivel de ajustar adequadamente.
#
# Por esse motivo, utilizamos Floresta Aleatoria (Random Forest)
# como classificador base para todos os metodos avaliados.
# Trata-se de um metodo nao-parametrico supervisionado, em
# contraste com a regressao logistica (parametrica) utilizada
# na base Wisconsin.
#
# Dessa forma garantimos comparacao justa entre:
#
#   - Maxima Relevancia
#   - mRMR
#   - mRMR + Backward
#   - mRMR + Forward
#   - LASSO (selecao) + Random Forest (classificacao)
#   - Stepwise AIC (selecao) + Random Forest (classificacao)
#
# Todos os metodos serao avaliados pela mesma metrica:
#
#        Erro = 1 - F1
#
# usando validacao cruzada estratificada de 5 dobras.
# ============================================================

suppressMessages({
  library(RaSEn)          # base Colon de Alon et al. (1999)
  library(caret)
  library(praznik)
  library(glmnet)         # LASSO (alpha=1)
  library(MASS)           # stepAIC
  library(randomForest)   # classificador base nao-parametrico
})

set.seed(42)

# ============================================================
# ETAPA 1: PREPARACAO DOS DADOS
# ============================================================

data(colon)

X_orig <- as.data.frame(colon$x)          # 62 x 2000
y_orig <- colon$y                          # 0 = normal, 1 = tumor

Y_fator <- factor(
  ifelse(y_orig == 1, "tumor", "normal"),
  levels = c("normal", "tumor")            # classe positiva: tumor
)

# Discretizacao para o mRMR (10 faixas de igual frequencia)
X_disc <- as.data.frame(
  lapply(X_orig, function(v) {
    as.integer(cut(v,
                   breaks = quantile(v,
                                     probs = seq(0, 1, length.out = 11),
                                     na.rm = TRUE),
                   include.lowest = TRUE,
                   labels = FALSE))
  })
)

preds <- colnames(X_orig)
M     <- length(preds)   # 2000

dados <- cbind(X_orig, Class = Y_fator)

cat("====================================\n")
cat("DADOS PREPARADOS (Colon — Peng 2005)\n")
cat("====================================\n")
cat("Observacoes :", nrow(dados), "\n")
cat("Preditoras  :", M, "\n")
cat("Resposta    : normal =", sum(Y_fator == "normal"),
    "| tumor =", sum(Y_fator == "tumor"), "\n\n")

# ============================================================
# ETAPA 2: FUNCOES AUXILIARES
# ============================================================

info_mutua <- function(x, y) {
  tbl <- table(x, y)
  n   <- sum(tbl)
  px  <- rowSums(tbl) / n
  py  <- colSums(tbl) / n
  mi  <- 0
  for (i in seq_along(px))
    for (j in seq_along(py)) {
      pxy <- tbl[i, j] / n
      if (pxy > 0)
        mi <- mi + pxy * log(pxy / (px[i] * py[j]))
    }
  max(mi, 0)
}

# ------------------------------------------------------------
# Classificador: Random Forest (nao-parametrico)
#
# Substituimos o Ridge Logistico por Random Forest.
# Razoes:
#   1) Nao requer regularizacao manual para p >> n;
#   2) Nao faz suposicoes distribucionais (nao-parametrico);
#   3) Robusto a multicolinearidade entre genes.
#
# Hiperparametros:
#   ntree = 500  (numero de arvores)
#   mtry  = floor(sqrt(k))  (padrao para classificacao)
#
# Retorna lista: erro = 1-F1 medio, dp = DP entre folds.
# ------------------------------------------------------------

erro_cv_rf <- function(variaveis,
                       dados,
                       resposta = "Class",
                       n_folds  = 5) {
  
  if (length(variaveis) == 0)
    return(list(erro = 1, dp = 0))
  
  set.seed(42)
  folds <- createFolds(dados[[resposta]], k = n_folds, list = TRUE)
  f1s   <- numeric(n_folds)
  
  for (i in seq_len(n_folds)) {
    treino <- dados[-folds[[i]], ]
    teste  <- dados[ folds[[i]], ]
    
    X_tr <- treino[, variaveis, drop = FALSE]
    y_tr <- treino[[resposta]]
    X_te <- teste[,  variaveis, drop = FALSE]
    
    tryCatch({
      rf <- randomForest(
        x     = X_tr,
        y     = y_tr,
        ntree = 500,
        mtry  = max(1, floor(sqrt(length(variaveis))))
      )
      pred <- predict(rf, newdata = X_te)
      
      cm     <- confusionMatrix(pred, teste[[resposta]], positive = "tumor")
      f1s[i] <- cm$byClass["F1"]
    }, error = function(e) f1s[i] <<- NA)
  }
  
  validos <- f1s[!is.na(f1s) & !is.nan(f1s)]
  if (length(validos) == 0) return(list(erro = 1, dp = 0))
  list(erro = 1 - mean(validos), dp = sd(validos))
}

# ============================================================
# ETAPA 3: PRIMEIRO ESTAGIO — mRMR
# ============================================================

K_MAX <- 50

Y_num    <- as.integer(Y_fator == "tumor")
res_mrmr <- MRMR(X = X_disc, Y = Y_num, k = K_MAX)

ordem_mrmr <- preds[res_mrmr$selection]
escores    <- res_mrmr$score

cat("Ordem de selecao mRMR (primeiras 20):\n")
for (i in 1:20)
  cat(sprintf("%2d. %-12s escore = %.5f\n",
              i, ordem_mrmr[i], escores[i]))
cat("...\n\n")

# ============================================================
# ETAPA 4: ESCOLHA DE n* VIA VALIDACAO CRUZADA
# ============================================================

cat("Calculando erros para k = 1 a", K_MAX, "...\n")

erros_mrmr <- numeric(K_MAX)
dps_mrmr   <- numeric(K_MAX)

for (j in seq_len(K_MAX)) {
  res            <- erro_cv_rf(ordem_mrmr[1:j], dados)
  erros_mrmr[j] <- res$erro
  dps_mrmr[j]   <- res$dp
  cat(sprintf("k = %2d | 1-F1 = %.4f | DP = %.4f\n",
              j, erros_mrmr[j], dps_mrmr[j]))
}

n_star  <- which.min(ifelse(is.nan(erros_mrmr), Inf, erros_mrmr))
S_cand  <- ordem_mrmr[1:n_star]
e_cand  <- erros_mrmr[n_star]
dp_cand <- dps_mrmr[n_star]

cat("\n====================================\n")
cat("MELHOR SUBCONJUNTO — mRMR 1o estagio\n")
cat("====================================\n")
cat("n* =", n_star, "\n")
cat("S = {", paste(S_cand, collapse = ", "), "}\n")
cat(sprintf("1-F1 = %.4f | DP = %.4f\n\n", e_cand, dp_cand))

# ============================================================
# ETAPA 5: BASELINE — MAXIMA RELEVANCIA
# ============================================================

cat("====================================\n")
cat("BASELINE: MAXIMA RELEVANCIA\n")
cat("====================================\n")

mi_alvo  <- sapply(preds, function(v) info_mutua(X_disc[[v]], Y_num))
ordem_mr <- names(sort(mi_alvo, decreasing = TRUE))

erros_mr <- numeric(K_MAX)
dps_mr   <- numeric(K_MAX)

for (j in seq_len(K_MAX)) {
  res          <- erro_cv_rf(ordem_mr[1:j], dados)
  erros_mr[j] <- res$erro
  dps_mr[j]   <- res$dp
}

n_mr  <- which.min(ifelse(is.nan(erros_mr), Inf, erros_mr))
S_mr  <- ordem_mr[1:n_mr]
e_mr  <- erros_mr[n_mr]
dp_mr <- dps_mr[n_mr]

cat(sprintf("Melhor k = %d | 1-F1 = %.4f | DP = %.4f\n\n",
            n_mr, e_mr, dp_mr))

# ============================================================
# ETAPA 6: REFINAMENTO — WRAPPER BACKWARD
# ============================================================

cat("====================================\n")
cat("REFINAMENTO: WRAPPER BACKWARD\n")
cat("====================================\n")

S_bw   <- S_cand
res_bw <- erro_cv_rf(S_bw, dados)
e_bw   <- res_bw$erro
dp_bw  <- res_bw$dp

mudou <- TRUE
while (mudou && length(S_bw) > 1) {
  mudou <- FALSE
  for (v in S_bw) {
    cand  <- setdiff(S_bw, v)
    res_c <- erro_cv_rf(cand, dados)
    if (res_c$erro <= e_bw) {
      cat(sprintf("Remove '%s': 1-F1 %.4f -> %.4f\n", v, e_bw, res_c$erro))
      S_bw  <- cand
      e_bw  <- res_c$erro
      dp_bw <- res_c$dp
      mudou <- TRUE
      break
    }
  }
}

cat("\nSubconjunto final Backward = {", paste(S_bw, collapse = ", "), "}\n")
cat(sprintf("1-F1 = %.4f | DP = %.4f\n\n", e_bw, dp_bw))

# ============================================================
# ETAPA 7: REFINAMENTO — WRAPPER FORWARD
# ============================================================

cat("====================================\n")
cat("REFINAMENTO: WRAPPER FORWARD\n")
cat("====================================\n")

S_fw  <- character(0)
e_fw  <- 1
dp_fw <- 0

for (v in S_cand) {
  cand  <- c(S_fw, v)
  res_c <- erro_cv_rf(cand, dados)
  ec    <- res_c$erro
  if (is.na(ec) || is.nan(ec)) {
    cat(sprintf("'%s' retornou NA — ignorado\n", v))
    next
  }
  if (ec <= e_fw) {
    cat(sprintf("Adiciona '%s': 1-F1 %.4f -> %.4f\n", v, e_fw, ec))
    S_fw  <- cand
    e_fw  <- ec
    dp_fw <- res_c$dp
  } else {
    cat(sprintf("'%s' nao melhora (%.4f >= %.4f)\n", v, ec, e_fw))
  }
}

cat("\nSubconjunto final Forward = {", paste(S_fw, collapse = ", "), "}\n")
cat(sprintf("1-F1 = %.4f | DP = %.4f\n\n", e_fw, dp_fw))

# ============================================================
# ETAPA 8: LASSO (selecao) + RANDOM FOREST (classificacao)
# ============================================================
#
# O LASSO e mantido como mecanismo de selecao de genes
# via penalizacao l1. A classificacao final, porem, e
# realizada com Random Forest sobre os genes selecionados.
#
# Isso isola o efeito da selecao do efeito do classificador,
# garantindo comparacao justa com os demais metodos.
# ============================================================

cat("====================================\n")
cat("LASSO (selecao) + RANDOM FOREST (classificacao)\n")
cat("====================================\n")

set.seed(42)
folds_ext <- createFolds(Y_fator, k = 5, list = TRUE)
f1s_lasso <- numeric(5)

for (i in seq_len(5)) {
  idx_tr <- setdiff(seq_len(nrow(dados)), folds_ext[[i]])
  
  X_tr   <- as.matrix(X_orig[idx_tr, ])
  y_tr   <- as.numeric(Y_fator[idx_tr] == "tumor")
  X_te   <- X_orig[folds_ext[[i]], ]
  y_te   <- Y_fator[folds_ext[[i]]]
  
  tryCatch({
    # LASSO para selecao de genes
    cv_las     <- cv.glmnet(X_tr, y_tr, family = "binomial",
                            alpha = 1, nfolds = 3,
                            type.measure = "class")
    coef_las   <- coef(cv_las, s = "lambda.min")[-1, 1]
    vars_fold  <- names(coef_las)[coef_las != 0]
    
    # Garantia: ao menos 1 gene selecionado
    if (length(vars_fold) == 0)
      vars_fold <- names(sort(abs(coef_las), decreasing = TRUE))[1]
    
    # Random Forest sobre genes selecionados pelo LASSO
    rf <- randomForest(
      x     = X_orig[idx_tr, vars_fold, drop = FALSE],
      y     = Y_fator[idx_tr],
      ntree = 500,
      mtry  = max(1, floor(sqrt(length(vars_fold))))
    )
    pred         <- predict(rf, newdata = X_te[, vars_fold, drop = FALSE])
    cm           <- confusionMatrix(pred, y_te, positive = "tumor")
    f1s_lasso[i] <- cm$byClass["F1"]
    
  }, error = function(e) f1s_lasso[i] <<- NA)
}

validos_lasso <- f1s_lasso[!is.na(f1s_lasso)]
e_lasso  <- 1 - mean(validos_lasso)
dp_lasso <- sd(validos_lasso)

# Genes selecionados no modelo completo (para relatorio)
cv_full    <- cv.glmnet(as.matrix(X_orig),
                        as.numeric(Y_fator == "tumor"),
                        family = "binomial", alpha = 1, nfolds = 5)
coef_full  <- coef(cv_full, s = "lambda.min")[-1, 1]
vars_lasso <- names(coef_full)[coef_full != 0]

cat(sprintf("1-F1 = %.4f | DP = %.4f\n", e_lasso, dp_lasso))
cat("Genes selecionados pelo LASSO (k =", length(vars_lasso), "):",
    paste(vars_lasso, collapse = ", "), "\n\n")

# ============================================================
# ETAPA 9: STEPWISE (stepAIC, selecao) + RANDOM FOREST
# ============================================================
#
# O Stepwise AIC e aplicado sobre os top-50 genes por MI
# apenas para selecionar variaveis. A classificacao final
# e feita com Random Forest, mantendo consistencia com os
# demais metodos.
# ============================================================

cat("====================================\n")
cat("STEPWISE (selecao, top-50 MI) + RANDOM FOREST\n")
cat("====================================\n")

top50    <- ordem_mr[1:50]
f1s_step <- numeric(5)

for (i in seq_len(5)) {
  idx_tr <- setdiff(seq_len(nrow(dados)), folds_ext[[i]])
  treino <- dados[idx_tr,         c(top50, "Class")]
  teste  <- dados[folds_ext[[i]], c(top50, "Class")]
  
  tryCatch({
    # Stepwise AIC para selecao de genes
    fm_nulo  <- glm(Class ~ 1,  data = treino, family = binomial())
    fm_full  <- glm(Class ~ .,  data = treino, family = binomial())
    mod_step <- suppressWarnings(
      stepAIC(fm_nulo,
              scope     = list(lower = fm_nulo, upper = fm_full),
              direction = "both",
              trace     = FALSE)
    )
    vars_fold <- names(coef(mod_step))[-1]
    
    if (length(vars_fold) == 0) {
      f1s_step[i] <- NA
      next
    }
    
    # Random Forest sobre genes selecionados pelo Stepwise
    rf <- randomForest(
      x     = treino[, vars_fold, drop = FALSE],
      y     = treino$Class,
      ntree = 500,
      mtry  = max(1, floor(sqrt(length(vars_fold))))
    )
    pred        <- predict(rf, newdata = teste[, vars_fold, drop = FALSE])
    cm          <- confusionMatrix(pred, teste$Class, positive = "tumor")
    f1s_step[i] <- cm$byClass["F1"]
    
  }, error = function(e) f1s_step[i] <<- NA)
}

validos_step <- f1s_step[!is.na(f1s_step)]
e_step  <- 1 - mean(validos_step)
dp_step <- sd(validos_step)

# Genes selecionados no modelo completo (para relatorio)
dados_top50   <- dados[, c(top50, "Class")]
mod_step_full <- suppressWarnings(
  stepAIC(
    glm(Class ~ 1, data = dados_top50, family = binomial()),
    scope     = list(
      lower = glm(Class ~ 1, data = dados_top50, family = binomial()),
      upper = glm(Class ~ ., data = dados_top50, family = binomial())
    ),
    direction = "both",
    trace     = FALSE
  )
)
vars_step <- names(coef(mod_step_full))[-1]

cat(sprintf("1-F1 = %.4f | DP = %.4f\n", e_step, dp_step))
cat("Genes selecionados pelo Stepwise (k =", length(vars_step), "):",
    paste(vars_step, collapse = ", "), "\n\n")

# ============================================================
# ETAPA 10: BENCHMARK — top-50 por MI + RANDOM FOREST
# ============================================================

res_top50 <- erro_cv_rf(top50, dados)
e_top50   <- res_top50$erro
dp_top50  <- res_top50$dp

cat(sprintf("Benchmark top-50 MI + RF: 1-F1 = %.4f | DP = %.4f\n\n",
            e_top50, dp_top50))

# ============================================================
# ETAPA 11: TABELA RESUMO FINAL
# ============================================================

cat("\n========================================================\n")
cat("TABELA RESUMO FINAL\n")
cat("========================================================\n")

tab <- data.frame(
  Metodo = c(
    paste0("Top-50 por MI (benchmark)           (k = 50)"),
    paste0("Maxima Relevancia                   (k = ", n_mr,          ")"),
    paste0("mRMR 1o estagio                     (k = ", n_star,        ")"),
    paste0("mRMR + Backward                     (k = ", length(S_bw),  ")"),
    paste0("mRMR + Forward                      (k = ", length(S_fw),  ")"),
    paste0("LASSO (sel.) + RF                   (k = ", length(vars_lasso), ")"),
    paste0("Stepwise AIC (sel.) + RF            (k = ", length(vars_step),  ")")
  ),
  k = c(50, n_mr, n_star, length(S_bw), length(S_fw),
        length(vars_lasso), length(vars_step)),
  Erro_1minusF1 = round(
    c(e_top50, e_mr, e_cand, e_bw, e_fw, e_lasso, e_step), 4),
  DP = round(
    c(dp_top50, dp_mr, dp_cand, dp_bw, dp_fw, dp_lasso, dp_step), 4)
)

print(tab, row.names = FALSE)

# ============================================================
# VARIAVEIS SELECIONADAS
# ============================================================

cat("\n========================================================\n")
cat("VARIAVEIS SELECIONADAS POR METODO\n")
cat("========================================================\n")
cat("Maxima Relevancia       :", paste(S_mr,       collapse = ", "), "\n")
cat("mRMR 1o estagio         :", paste(S_cand,     collapse = ", "), "\n")
cat("mRMR + Backward         :", paste(S_bw,       collapse = ", "), "\n")
cat("mRMR + Forward          :", paste(S_fw,       collapse = ", "), "\n")
cat("LASSO (selecao)         :", paste(vars_lasso, collapse = ", "), "\n")
cat("Stepwise AIC (selecao)  :", paste(vars_step,  collapse = ", "), "\n")