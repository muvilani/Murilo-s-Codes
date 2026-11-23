# GERAÇÃO DA REDE NEURONAL (MODELO GL) COM 6 NEURÔNIOS
# E APLICAÇÃO DO mRMR PARA SELECIONAR VIZINHANÇAS

# -----------------------------------------------------------------------------
# 0. PACOTES
# -----------------------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(infotheo)
library(praznik)

# -----------------------------------------------------------------------------
# 1. PARÂMETROS DO MODELO
# -----------------------------------------------------------------------------

set.seed(42)

n <- 6   # número de neurônios

# Matriz sináptica W (6 x 6), SIMÉTRICA, peso 10 onde há conexão.
# W[i, j] = peso do neurônio j sobre o neurônio i.
# Cada neurônio recebe de exatamente 2 pré-sinápticos:
#   N1 <- N2, N3  |  N2 <- N1, N4  |  N3 <- N1, N5
#   N4 <- N2, N6  |  N5 <- N3, N6  |  N6 <- N4, N5
w <- matrix(
  c(  0, 10, 10,  0,  0,  0,
      10,  0,  0, 10,  0,  0,
      10,  0,  0,  0, 10,  0,
      0, 10,  0,  0,  0, 10,
      0,  0, 10,  0,  0, 10,
      0,  0,  0, 10, 10,  0),
  nrow = n, ncol = n, byrow = TRUE
)

alpha    <- rep(0, n)     # atividade espontânea
t_burnin <- 2500
t_amostra <- 40000
t_total  <- t_burnin + t_amostra

X_init <- matrix(
  c(1, 0, 0, 1,
    0, 1, 1, 0,
    1, 0, 1, 0,
    0, 1, 0, 1,
    1, 1, 0, 0,
    0, 0, 1, 1),
  nrow = n, byrow = TRUE
)

# -----------------------------------------------------------------------------
# 2. FUNÇÕES AUXILIARES DO MODELO GL
# -----------------------------------------------------------------------------

phi <- function(x) exp(x) / (1 + exp(x))

# Retorna a linha 'indice' da matriz m (sempre como matrix)
matriz_x <- function(m, indice) m[indice, , drop = FALSE]

# Retorna todas as linhas de m EXCETO 'indice' (sempre como matrix)
matriz_pesos <- function(m, indice) m[!(1:nrow(m) %in% indice), , drop = FALSE]

# Índice do último disparo até o instante 'tempo_final'
ultimo <- function(x, tempo_inicial, tempo_final) {
  tail(which(x[tempo_inicial:tempo_final] == 1), n = 1)
}

# Potencial de membrana — versão corrigida
# Correções aplicadas em relação ao código original:
#   1. drop = FALSE em todas as submatrizes para evitar colapso para vetor
#   2. No branch 'else', rep() usa (tempo - indice), não (tempo - indice + 1)
V_novo <- function(X, indice_neuronio, tempo) {
  
  x_i    <- matriz_x(X, indice_neuronio)          # linha do neurônio i
  indice <- ultimo(x_i, 1, tempo)                 # último disparo
  
  # Caso 1: disparou no instante atual → potencial = 0
  if (length(indice) == 0 || indice == tempo) {
    return(0)
  }
  
  # Janela de disparos dos outros neurônios: colunas (indice+1) até tempo
  janela <- (indice + 1):tempo
  n_cols <- length(janela)                        # tamanho da janela
  
  # Submatriz dos pré-sinápticos na janela — mantém sempre como matrix
  X_pre  <- matriz_pesos(X, indice_neuronio)[, janela, drop = FALSE]
  
  # Coluna de pesos do neurônio indice_neuronio (coluna fixa do i)
  w_pre  <- matriz_pesos(w, indice_neuronio)[, indice_neuronio, drop = FALSE]
  
  # Fator de decaimento: 1 / 2^(tempo - indice)
  fator  <- 1 / (2^(tempo - indice))
  
  # Potencial: fator * soma dos disparos pré-sinápticos ponderados por w
  # X_pre é (n-1) x n_cols; w_pre é (n-1) x 1
  # rowSums(X_pre) é (n-1); t(w_pre) é 1 x (n-1)
  potencial <- fator * sum(rowSums(X_pre) * as.vector(w_pre))
  
  return(potencial)
}

# Geração da cadeia X
gera_cadeia <- function(X) {
  for (j in 5:ncol(X)) {
    for (i in 1:n) {
      X[i, j] <- ifelse(
        runif(1, 0, 1) <= phi(alpha[i] + V_novo(X, i, j - 1)),
        1, 0
      )
    }
  }
  return(X)
}

# -----------------------------------------------------------------------------
# 3. GERAÇÃO DA CADEIA E BURN-IN
# -----------------------------------------------------------------------------

cat("Gerando cadeia (pode demorar alguns minutos)...\n")

X_full <- cbind(X_init, matrix(0, nrow = n, ncol = t_total - 4))
X_full <- gera_cadeia(X_full)

cat("Cadeia gerada!\n\n")

# Proporção de disparos (cadeia completa)
cat("Proporção de disparos por neurônio (cadeia completa):\n")
for (i in 1:n) cat(sprintf("  N%d: %.3f\n", i, mean(X_full[i, ])))

# Descarta burn-in
X <- X_full[, (t_burnin + 1):t_total]

cat(sprintf("\nAmostra após burn-in: %d neurônios x %d instantes\n",
            nrow(X), ncol(X)))
cat("Proporção de disparos após burn-in:\n")
for (i in 1:n) cat(sprintf("  N%d: %.3f\n", i, mean(X[i, ])))

# -----------------------------------------------------------------------------
# 4. GRÁFICO DE BURN-IN
# -----------------------------------------------------------------------------

prop_acum <- matrix(0, nrow = n, ncol = ncol(X_full))
for (i in 1:n)
  prop_acum[i, ] <- cumsum(X_full[i, ]) / seq_along(X_full[i, ])

df_bi <- as.data.frame(t(prop_acum))
colnames(df_bi) <- paste0("N", 1:n)
df_bi$Tempo <- seq_len(ncol(X_full))

df_long <- pivot_longer(df_bi, cols = -Tempo,
                        names_to = "Neuronio", values_to = "Proporcao")

p_bi <- ggplot(df_long, aes(x = Tempo, y = Proporcao, color = Neuronio)) +
  geom_line(linewidth = 0.4) +
  geom_vline(xintercept = t_burnin, linetype = "dashed", color = "black") +
  annotate("text", x = t_burnin + 300, y = 0.93,
           label = "Burn-in", hjust = 0, size = 3.5) +
  labs(x = "Tempo (bin)", y = "Proporção acumulada de disparos",
       title = "Verificação do burn-in — rede com 6 neurônios") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("burnin_6neuronios.pdf", plot = p_bi, width = 8, height = 4)
cat("\nGráfico salvo em 'burnin_6neuronios.pdf'\n")

# -----------------------------------------------------------------------------
# 5. CONSTRUÇÃO DAS VARIÁVEIS U_t^{(j -> i)}
# -----------------------------------------------------------------------------
# U[[i]][[j]]: vetor de tamanho T_obs com a atividade recente de j
# em relação ao neurônio alvo i.

T_obs <- ncol(X)

cat("\nCalculando variáveis U_t^{(j->i)}...\n")

# Pré-calcula tau_{i,t} = último instante de disparo de i antes de t
calc_tau <- function(X_i) {
  tau        <- integer(length(X_i))
  ultimo_d   <- 0L
  for (t in seq_along(X_i)) {
    tau[t]   <- ultimo_d
    if (X_i[t] == 1L) ultimo_d <- t
  }
  tau
}

U <- vector("list", n)
for (i in 1:n) {
  tau_i  <- calc_tau(X[i, ])
  U[[i]] <- vector("list", n)
  for (j in (1:n)[-i]) {
    U_ij <- integer(T_obs)
    for (t in 2:T_obs) {
      tau_it <- tau_i[t]
      if (tau_it > 0 && tau_it < t - 1) {
        U_ij[t] <- sum(X[j, (tau_it + 1):(t - 1)])
      }
    }
    U[[i]][[j]] <- U_ij
  }
}

cat("Variáveis U calculadas.\n")

# -----------------------------------------------------------------------------
# 6. mRMR: FUNÇÕES DE RELEVÂNCIA E REDUNDÂNCIA
# -----------------------------------------------------------------------------

# Relevância marginal: I(U_{j->i} ; X_i)
calc_relevancia <- function(U, X, i, n) {
  rel <- rep(NA_real_, n)
  Xi  <- as.factor(X[i, ])
  for (j in (1:n)[-i]) {
    Uij    <- as.factor(U[[i]][[j]])
    rel[j] <- mutinformation(Uij, Xi, method = "emp")
  }
  rel
}

# Relevância incremental: I(U_{j->i} ; X_i | U_{S->i})
calc_inc <- function(U, X, i, j, S) {
  Xi  <- as.factor(X[i, ])
  Uij <- as.factor(U[[i]][[j]])
  
  if (length(S) == 0) {
    return(mutinformation(Uij, Xi, method = "emp"))
  }
  
  # Estado conjunto de U_S
  US_mat   <- sapply(S, function(k) as.integer(U[[i]][[k]]))
  US_joint <- apply(US_mat, 1, paste, collapse = "_")
  
  UjS_joint <- paste(as.integer(U[[i]][[j]]), US_joint, sep = "_")
  
  H_Xi_US  <- condentropy(Xi, as.factor(US_joint),  method = "emp")
  H_Xi_UjS <- condentropy(Xi, as.factor(UjS_joint), method = "emp")
  
  max(H_Xi_US - H_Xi_UjS, 0)
}

# -----------------------------------------------------------------------------
# 7. ALGORITMO mRMR PARA CADA NEURÔNIO
# -----------------------------------------------------------------------------

mrmr_neuronal <- function(U, X, i, n,
                          lambda = 1, M_max = 2) {
  candidatos <- (1:n)[-i]
  rel        <- calc_relevancia(U, X, i, n)
 
  S <- integer(0)
  
  for (m in seq_len(M_max)) {
    disponiveis <- setdiff(candidatos, S)
    if (length(disponiveis) == 0) break
    
    if (m == 1) {
      escore <- rel[disponiveis]
    } else {
      escore <- sapply(disponiveis, function(j) {
        inc_j <- calc_inc(U, X, i, j, S)
        red_j <- rel[j] - inc_j
        rel[j] - lambda * red_j
      })
    }
    
    j_star <- disponiveis[which.max(escore)]
    
    
    S <- c(S, j_star)
  }
  
  list(S = S, relevancia = rel)
}

cat("\nAplicando mRMR para cada neurônio...\n")

resultados <- vector("list", n)
for (i in 1:n) {
  cat(sprintf("  Neurônio %d...\n", i))
  resultados[[i]] <- mrmr_neuronal(U, X, i, n,
                                   lambda = 1, M_max = 2)
}

# -----------------------------------------------------------------------------
# 8. RESULTADOS: VIZINHANÇA ESTIMADA vs. VERDADEIRA
# -----------------------------------------------------------------------------

vizinhanca_verd <- lapply(1:n, function(i) sort(which(w[i, ] != 0)))

cat("\n=== RESULTADOS DO mRMR ===\n")
cat(sprintf("%-10s %-20s %-20s %-8s\n",
            "Neurônio", "Viz. Verdadeira", "Viz. Estimada", "Correto?"))
cat(strrep("-", 60), "\n")

for (i in 1:n) {
  vv  <- vizinhanca_verd[[i]]
  ve  <- sort(resultados[[i]]$S)
  ok  <- ifelse(identical(vv, ve), "SIM", "NAO")
  cat(sprintf("N%-9d {%s}%*s {%s}%*s %-8s\n",
              i,
              paste(vv, collapse=","), max(0, 14 - nchar(paste(vv,collapse=","))), "",
              paste(ve, collapse=","), max(0, 14 - nchar(paste(ve,collapse=","))), "",
              ok))
}

# -----------------------------------------------------------------------------
# 9. ESCORES DE RELEVÂNCIA INCREMENTAL W_ij
# -----------------------------------------------------------------------------

cat("\n=== ESCORES W_ij (relevância incremental) ===\n")
for (i in 1:n) {
  S_i <- resultados[[i]]$S
  if (length(S_i) == 0) next
  cat(sprintf("Neurônio %d (pós-sináptico):\n", i))
  escores <- numeric(length(S_i))
  names(escores) <- S_i
  for (k in seq_along(S_i)) {
    j       <- S_i[k]
    S_sem_j <- setdiff(S_i, j)
    w_ij    <- calc_inc(U, X, i, j, S_sem_j)
    escores[k] <- w_ij
    cat(sprintf("  N%d -> N%d : W = %.6f\n", j, i, w_ij))
  }
  # Classificação forte/fraca
  thr <- quantile(escores, 0.75)
  cat(sprintf("  Limiar (Q75) = %.6f\n", thr))
  for (k in seq_along(S_i)) {
    j    <- S_i[k]
    cls  <- ifelse(escores[k] >= thr, "FORTE", "fraca")
    cat(sprintf("  N%d -> N%d : %s\n", j, i, cls))
  }
}

# -----------------------------------------------------------------------------
# 10. MATRIZ DE ADJACÊNCIA E MÉTRICAS
# -----------------------------------------------------------------------------

A_hat  <- matrix(0, nrow = n, ncol = n)
for (i in 1:n)
  for (j in resultados[[i]]$S)
    A_hat[j, i] <- 1

A_true <- (w != 0) * 1

cat("\n=== MATRIZ DE ADJACÊNCIA ESTIMADA ===\n")
print(A_hat)

cat("\n=== MATRIZ DE ADJACÊNCIA VERDADEIRA ===\n")
print(A_true)

# Métricas (exclui diagonal)
diag_idx <- which(diag(n) == 1)
VP <- sum(A_hat[-diag_idx] == 1 & A_true[-diag_idx] == 1)
FP <- sum(A_hat[-diag_idx] == 1 & A_true[-diag_idx] == 0)
FN <- sum(A_hat[-diag_idx] == 0 & A_true[-diag_idx] == 1)

prec <- VP / max(VP + FP, 1)
sens <- VP / max(VP + FN, 1)
f1   <- ifelse(prec + sens > 0,
               2 * prec * sens / (prec + sens), 0)

cat(sprintf("\nVP=%d  FP=%d  FN=%d\n", VP, FP, FN))
cat(sprintf("Precisão     : %.3f\n", prec))
cat(sprintf("Sensibilidade: %.3f\n", sens))
cat(sprintf("F1-score     : %.3f\n", f1))
cat("\n=== CONCLUÍDO ===\n")