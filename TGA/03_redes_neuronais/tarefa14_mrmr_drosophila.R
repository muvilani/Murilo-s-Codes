# =============================================================================
# APLICAÇÃO DO mRMR EM DADOS REAIS DO CONECTOMA DA DROSOPHILA
#
# Circuito reduzido formado por 6 neurônios extraídos do conectoma público
# do sistema nervoso central da drosófila macho (male-cns:v1.0).
#
# Fonte dos dados:
# neuPrint (Janelia Research Campus)
# Dataset: male-cns:v1.0
#
# Seleção dos neurônios:
#
# 1) Foram identificados neurônios dos tipos PPL, PAM e MBON.
# 2) Para cada tipo neuronal foi selecionado um representante,
#    escolhido como o neurônio com maior grau total de conectividade.
# 3) As conexões entre neurônios de entrada (PPL/PAM) e neurônios MBON
#    foram ranqueadas segundo o peso sináptico.
# 4) Foram selecionados os três pares mais fortemente conectados:
#
#       PPL101 -> MBON11   (610 sinapses)
#       PPL103 -> MBON32   (242 sinapses)
#       PPL105 -> MBON13   (199 sinapses)
#
# O circuito final utilizado nos experimentos é composto pelos neurônios:
#
#       PPL101
#       PPL103
#       PPL105
#       MBON11
#       MBON32
#       MBON13
#
# A matriz W contém os pesos sinápticos reais observados nesse subcircuito.
# O restante do código segue exatamente a mesma estrutura utilizada no
# exemplo simulado do modelo de Galves-Löcherbach.
# =============================================================================

library(readr)
drosophila_mb_connections <- read_csv("C:/Users/muvil/Downloads/drosophila_mb_connections.csv")
drosophila_mb_circuit_6neurons <- read_csv("C:/Users/muvil/Downloads/drosophila_mb_circuit_6neurons.csv")

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

# Matriz sináptica W (6 x 6) obtida do conectoma male-cns:v1.0.
#
# Ordem dos neurônios:
#
# 1 = PPL101
# 2 = PPL103
# 3 = PPL105
# 4 = MBON11
# 5 = MBON32
# 6 = MBON13
#
# Cada elemento W[i,j] representa o número de sinapses do neurônio j
# sobre o neurônio i.
#
# Convenção utilizada:
#
#       coluna j = neurônio pré-sináptico
#       linha  i = neurônio pós-sináptico
#
# Essa convenção coincide com a utilizada no exemplo simulado do modelo GL.
#
# A matriz foi inicialmente extraída via neuPrint e posteriormente
# transposta para adequação à notação adotada neste trabalho.
#
# Os pesos correspondem ao número de contatos sinápticos observados
# entre os neurônios selecionados no conectoma real.


library(readr)
tmp <- read_csv(
  "C:/Users/muvil/Downloads/matriz_male_cns_6neurons.csv",
  show_col_types = FALSE
)

w <- as.matrix(tmp[, -1])
storage.mode(w) <- "numeric"

# Identificação dos neurônios selecionados no conectoma.

neuron_labels <- c(
  "PPL101",
  "PPL103",
  "PPL105",
  "MBON11",
  "MBON32",
  "MBON13"
)
# Classificação utilizada ao longo do trabalho:
# I = neurônio de entrada do circuito
# O = neurônio de saída do circuito

neuron_types  <- c("I", "I", "I", "O", "O", "O")

alpha    <- rep(0, n)     # atividade espontânea (fixada em zero)
t_burnin  <- 2500
t_amostra <- 40000
t_total   <- t_burnin + t_amostra

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

matriz_x     <- function(m, indice) m[indice, , drop = FALSE]
matriz_pesos <- function(m, indice) m[!(1:nrow(m) %in% indice), , drop = FALSE]

ultimo <- function(x, tempo_inicial, tempo_final) {
  tail(which(x[tempo_inicial:tempo_final] == 1), n = 1)
}

V_novo <- function(X, indice_neuronio, tempo) {
  
  x_i    <- matriz_x(X, indice_neuronio)
  indice <- ultimo(x_i, 1, tempo)
  
  if (length(indice) == 0 || indice == tempo) {
    return(0)
  }
  
  janela <- (indice + 1):tempo
  X_pre  <- matriz_pesos(X, indice_neuronio)[, janela, drop = FALSE]
  w_pre  <- matriz_pesos(w, indice_neuronio)[, indice_neuronio, drop = FALSE]
  fator  <- 1 / (2^(tempo - indice))
  potencial <- fator * sum(rowSums(X_pre) * as.vector(w_pre))
  
  return(potencial)
}

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

# Simulação de uma cadeia do modelo GL utilizando os pesos
# sinápticos observados no conectoma real.

X_full <- cbind(X_init, matrix(0, nrow = n, ncol = t_total - 4))
X_full <- gera_cadeia(X_full)

cat("Cadeia gerada!\n\n")

cat("Proporção de disparos por neurônio (cadeia completa):\n")
for (i in 1:n)
  cat(sprintf("  %s: %.3f\n", neuron_labels[i], mean(X_full[i, ])))

X <- X_full[, (t_burnin + 1):t_total]

cat(sprintf("\nAmostra após burn-in: %d neurônios x %d instantes\n",
            nrow(X), ncol(X)))
cat("Proporção de disparos após burn-in:\n")
for (i in 1:n)
  cat(sprintf("  %s: %.3f\n", neuron_labels[i], mean(X[i, ])))

# -----------------------------------------------------------------------------
# 4. GRÁFICO DE BURN-IN
# -----------------------------------------------------------------------------

prop_acum <- matrix(0, nrow = n, ncol = ncol(X_full))
for (i in 1:n)
  prop_acum[i, ] <- cumsum(X_full[i, ]) / seq_along(X_full[i, ])

df_bi <- as.data.frame(t(prop_acum))
colnames(df_bi) <- neuron_labels
df_bi$Tempo <- seq_len(ncol(X_full))

df_long <- pivot_longer(df_bi, cols = -Tempo,
                        names_to = "Neuronio", values_to = "Proporcao")

p_bi <- ggplot(
  df_long,
  aes(
    x = Tempo,
    y = Proporcao,
    color = Neuronio,
    linetype = Neuronio
  )
) +
  geom_line(linewidth = 0.4) +
  geom_vline(
    xintercept = t_burnin,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = t_burnin + 300,
    y = 0.93,
    label = "Burn-in"
  ) +
  labs(
    title = "Verificação do burn-in — circuito real da Drosophila (6 neurônios)",
    x = "Tempo",
    y = "Proporção acumulada de disparos",
    color = "Neurônio",
    linetype = "Neurônio"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "C:/Users/muvil/Downloads/burnin_tcc.jpeg",
  plot = p_bi,
  width = 8,
  height = 4
)


# -----------------------------------------------------------------------------
# 5. CONSTRUÇÃO DAS VARIÁVEIS U_t^{(j -> i)}
# -----------------------------------------------------------------------------

T_obs <- ncol(X)

cat("\nCalculando variáveis U_t^{(j->i)}...\n")

calc_tau <- function(X_i) {
  tau      <- integer(length(X_i))
  ultimo_d <- 0L
  for (t in seq_along(X_i)) {
    tau[t] <- ultimo_d
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

for (i in 1:n) {
  
  candidatos <- setdiff(1:n, i)
  
  X_df <- as.data.frame(
    lapply(candidatos, function(j) as.integer(U[[i]][[j]]))
  )
  colnames(X_df) <- paste0("N", candidatos)
  
  Y <- as.integer(X[i, ])
  
  res <- MRMR(
    X = X_df,
    Y = Y,
    k = length(candidatos)
  )
  
  ranking <- candidatos[res$selection]
  
  # Tamanho da vizinhança verdadeira de i
  viz_verd  <- which(w[i, ] != 0)
  M_i       <- length(viz_verd)
  
  resultados[[i]] <- list(
    ranking      = ranking,
    selecionados = ranking[1:M_i],   # seleciona M_i neurônios (viz. verdadeira)
    M_i          = M_i,
    score        = res$score
  )
  
  cat(sprintf("%s -> ranking: %s  (top-%d selecionados: %s)\n",
              neuron_labels[i],
              paste(ranking, collapse = " "),
              M_i,
              paste(ranking[1:M_i], collapse = " ")))
}

# -----------------------------------------------------------------------------
# 7. RESULTADOS
# -----------------------------------------------------------------------------
# A vizinhança verdadeira de cada neurônio é determinada diretamente
# pela matriz W: um neurônio j pertence à vizinhança de i quando
# W[i,j] > 0, ou seja, quando existe conexão sináptica j -> i.

vizinhanca_verd <- lapply(1:n, function(i) sort(which(w[i, ] != 0)))

cat("\n=== RESULTADOS DO mRMR ===\n")

for (i in 1:n) {
  
  vv <- sort(vizinhanca_verd[[i]])
  ve <- sort(resultados[[i]]$selecionados)
  ok <- identical(vv, ve)
  
  cat("\n")
  cat(sprintf("%s\n", neuron_labels[i]))
  cat("Verdadeira:", paste(vv, collapse = ","), "\n")
  cat("Estimada  :", paste(ve, collapse = ","), "\n")
}

# -----------------------------------------------------------------------------
# 8. RANKING COMPLETO
# -----------------------------------------------------------------------------

cat("\n=== RANKING COMPLETO ===\n")

for (i in 1:n) {
  cat(sprintf("\n%s\n", neuron_labels[i]))
  print(
    data.frame(
      Neuronio = resultados[[i]]$ranking,
      Tipo     = neuron_types[resultados[[i]]$ranking],
      Score    = resultados[[i]]$score
    )
  )
}

# -----------------------------------------------------------------------------
# 9. MATRIZ DE ADJACÊNCIA E MÉTRICAS
# -----------------------------------------------------------------------------

A_hat <- matrix(0, n, n)
for (i in 1:n) {
  viz <- resultados[[i]]$selecionados
  A_hat[i, viz] <- 1   # linha = alvo (i), coluna = pré-sináptico (viz) — mesma convenção de A_true
}

# A_true representa a estrutura real do conectoma reduzido.
# A_hat representa a estrutura inferida pelo mRMR.
A_true <- (w != 0) * 1
diag(A_hat)  <- 0
diag(A_true) <- 0

VP <- sum(A_hat == 1 & A_true == 1)
FP <- sum(A_hat == 1 & A_true == 0)
FN <- sum(A_hat == 0 & A_true == 1)

prec <- VP / (VP + FP)
sens <- VP / (VP + FN)
f1   <- ifelse(prec + sens > 0,
               2 * prec * sens / (prec + sens), 0)

cat("\n")
cat("VP =", VP, "\n")
cat("FP =", FP, "\n")
cat("FN =", FN, "\n\n")
cat("Precisão     =", round(prec, 4), "\n")
cat("Sensibilidade=", round(sens, 4), "\n")
cat("F1           =", round(f1,   4), "\n")

# -----------------------------------------------------------------------------
# 10. COMPARAÇÃO ENTRE MATRIZ VERDADEIRA E ESTIMADA
# -----------------------------------------------------------------------------

cat("\n=== MATRIZES DE ADJACÊNCIA ===\n")
cat("\nVerdadeira (A_true):\n")
rownames(A_true) <- colnames(A_true) <- neuron_labels
print(A_true)

cat("\nEstimada (A_hat):\n")
rownames(A_hat) <- colnames(A_hat) <- neuron_labels
print(A_hat)

# Matriz recuperada com pesos sinápticos REAIS nas posições detectadas
W_recuperada <- matrix(0, n, n)
rownames(W_recuperada) <- colnames(W_recuperada) <- neuron_labels

for (i in 1:n) {
  viz <- resultados[[i]]$selecionados
  W_recuperada[viz, i] <- w[i, viz]   # coloca o peso REAL onde detectou conexão
}

cat("\nMatriz W verdadeira:\n")
print(w)

cat("\nMatriz W recuperada (pesos reais nas posições detectadas):\n")
print(W_recuperada)

W_recuperada <- matrix(0, n, n)
rownames(W_recuperada) <- colnames(W_recuperada) <- neuron_labels

for (i in 1:n) {
  viz <- resultados[[i]]$selecionados
  W_recuperada[i, viz] <- w[i, viz]   # linha = pós-sináptico, coluna = pré-sináptico
}

cat("\nMatriz W verdadeira:\n")
print(w)

cat("\nMatriz W recuperada:\n")
print(W_recuperada)



