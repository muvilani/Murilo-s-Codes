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

# Matriz sináptica W (6 x 6), SIMÉTRICA, peso 10 onde há conexão e 0 c.c..
# W[i, j] = peso do neurônio j sobre o neurônio i. -> influência do neurônio j sobre o neurônio i
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

alpha    <- rep(0, n)     # atividade espontânea (O parâmetro de atividade espontânea foi fixado em zero para todos os neurônios)
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
) # estados iniciais dos neurônios nos quatro primeiros instantes.

# -----------------------------------------------------------------------------
# 2. FUNÇÕES AUXILIARES DO MODELO GL
# -----------------------------------------------------------------------------

phi <- function(x) exp(x) / (1 + exp(x))

# Função auxiliar que devolve a linha correspondente ao neurônio analisado.
matriz_x <- function(m, indice) m[indice, , drop = FALSE]

# Retorna todos os outros neurônios, exceto o neurônio alvo.
matriz_pesos <- function(m, indice) m[!(1:nrow(m) %in% indice), , drop = FALSE]

# Localiza o último instante em que um neurônio disparou.
ultimo <- function(x, tempo_inicial, tempo_final) {
  tail(which(x[tempo_inicial:tempo_final] == 1), n = 1)
}
# exemplo: 0 1 0 0 1 0 (posição 5)
# Potencial de membrana — versão corrigida
# Correções aplicadas em relação ao código original:
#   1. drop = FALSE em todas as submatrizes para evitar colapso para vetor
#   2. No branch 'else', rep() usa (tempo - indice), não (tempo - indice + 1)
V_novo <- function(X, indice_neuronio, tempo) {
  
  x_i    <- matriz_x(X, indice_neuronio)          # pegando o histórico do neurônio i
  indice <- ultimo(x_i, 1, tempo)                 # último disparo
  
  # Caso 1: Se ele acabou de disparar, o potencial é zerado.0
  if (length(indice) == 0 || indice == tempo) {
    return(0)
  }
  
  # Janela de disparos dos outros neurônios: Defino a janela de observação após o último disparo
  janela <- (indice + 1):tempo
  n_cols <- length(janela)                        # tamanho da janela
  
  # Recupero os disparos dos neurônios pré-sinápticos nessa janela.
  X_pre  <- matriz_pesos(X, indice_neuronio)[, janela, drop = FALSE]
  
  # Recupero os pesos das conexões.
  w_pre  <- matriz_pesos(w, indice_neuronio)[, indice_neuronio, drop = FALSE]
  
  # Fator de decaimento: 1 / 2^(tempo - indice) -> Quanto mais antigo o estímulo, menor sua influência.
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

ggsave("C:\\Users\\muvil\\Downloads\\burnin_6neuronios.jpeg", plot = p_bi, width = 8, height = 4)
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
# 6. mRMR VIA PRAZNIK
# -----------------------------------------------------------------------------

cat("\nAplicando mRMR (praznik)...\n")

resultados <- vector("list", n)

for(i in 1:n){
  
  candidatos <- setdiff(1:n, i)
  
  X_df <- as.data.frame(
    lapply(
      candidatos,
      function(j) as.integer(U[[i]][[j]])
    )
  )
  
  colnames(X_df) <- paste0("N", candidatos)
  
  Y <- as.integer(X[i, ])
  
  res <- MRMR(
    X = X_df,
    Y = Y,
    k = length(candidatos)
  )
  
  ranking <- candidatos[res$selection]
  
  resultados[[i]] <- list(
    ranking = ranking,
    selecionados = ranking[1:2],
    score = res$score
  )
  
  cat(
    sprintf(
      "Neurônio %d -> %s\n",
      i,
      paste(ranking, collapse = " ")
    )
  )
}

# -----------------------------------------------------------------------------
# 7. RESULTADOS
# -----------------------------------------------------------------------------

vizinhanca_verd <- lapply(
  1:n,
  function(i) sort(which(w[i, ] != 0))
)

cat("\n=== RESULTADOS DO mRMR ===\n")

for(i in 1:n){
  
  vv <- sort(vizinhanca_verd[[i]])
  
  ve <- sort(resultados[[i]]$selecionados)
  
  ok <- identical(vv, ve)
  
  cat("\n")
  
  cat(sprintf("Neurônio %d\n", i))
  
  cat(
    "Verdadeira:",
    paste(vv, collapse = ","),
    "\n"
  )
  
  cat(
    "Estimada  :",
    paste(ve, collapse = ","),
    "\n"
  )
  
  cat(
    "Correto   :",
    ifelse(ok, "SIM", "NAO"),
    "\n"
  )
}

# -----------------------------------------------------------------------------
# 8. RANKING COMPLETO
# -----------------------------------------------------------------------------

cat("\n=== RANKING COMPLETO ===\n")

for(i in 1:n){
  
  cat(sprintf("\nNeurônio %d\n", i))
  
  print(
    data.frame(
      Neuronio = resultados[[i]]$ranking,
      Score = resultados[[i]]$score
    )
  )
}

# -----------------------------------------------------------------------------
# 9. MATRIZ DE ADJACÊNCIA
# -----------------------------------------------------------------------------

A_hat <- matrix(0, n, n)

for(i in 1:n){
  
  viz <- resultados[[i]]$selecionados
  
  A_hat[viz, i] <- 1
}

A_true <- (w != 0) * 1

diag(A_hat) <- 0
diag(A_true) <- 0

VP <- sum(A_hat == 1 & A_true == 1)
FP <- sum(A_hat == 1 & A_true == 0)
FN <- sum(A_hat == 0 & A_true == 1)

prec <- VP/(VP + FP)
sens <- VP/(VP + FN)

f1 <- ifelse(
  prec + sens > 0,
  2 * prec * sens/(prec + sens),
  0
)

cat("\n")

cat("VP =", VP, "\n")
cat("FP =", FP, "\n")
cat("FN =", FN, "\n\n")

cat("Precisão     =", round(prec, 4), "\n")
cat("Sensibilidade=", round(sens, 4), "\n")
cat("F1           =", round(f1, 4), "\n")

resultados
vizinhanca_verd

# -----------------------------------------------------------------------------
# 10. MATRIZ DE PESOS W_tilde (scores mRMR)
# -----------------------------------------------------------------------------

W_tilde <- matrix(0, nrow = n, ncol = n)
rownames(W_tilde) <- paste0("N", 1:n)
colnames(W_tilde) <- paste0("N", 1:n)

for (i in 1:n) {
  candidatos <- setdiff(1:n, i)
  scores     <- resultados[[i]]$score        # nomeado por "N1", "N2", ...
  for (j in candidatos) {
    nome_j          <- paste0("N", j)
    W_tilde[j, i]  <- scores[nome_j]         # W_tilde[j, i] = score de j -> i
  }
}

cat("\n=== MATRIZ W_tilde (scores mRMR, j -> i) ===\n")
print(round(W_tilde, 4))

cat("\n=== COMPARAÇÃO: W_tilde vs W_verdadeira ===\n")
cat("\nConexões verdadeiras (W != 0):\n")
print((w != 0) * 1)

cat("\nConexões recuperadas (W_tilde > 0):\n")
print((W_tilde > 0) * 1)

