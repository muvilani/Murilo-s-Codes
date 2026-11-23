# =============================================================================
# GERAÇÃO DA REDE NEURONAL (MODELO GL) COM 6 NEURÔNIOS
# E APLICAÇÃO DO mRMR PARA SELECIONAR VIZINHANÇAS
#
# Autor: Murilo Cassiavilani
# Orientador: Prof. Dr. Ricardo Felipe Ferreira
#
# Referência do modelo GL: Galves & Löcherbach (2013, 2016)
# Referência do mRMR:      Peng et al. (2005)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACOTES
# -----------------------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(infotheo)   # mutinformation(), discretize()

# -----------------------------------------------------------------------------
# 1. PARÂMETROS DO MODELO
# -----------------------------------------------------------------------------

set.seed(42)

n <- 6   # número de neurônios

# Matriz sináptica W (6 x 6), SIMÉTRICA, peso 10 onde há conexão.
# Leitura: W[i, j] = peso do neurônio j sobre o neurônio i.
# Cada neurônio recebe de exatamente 2 pré-sinápticos:
#   N1 <- N2, N3
#   N2 <- N1, N4
#   N3 <- N1, N5
#   N4 <- N2, N6
#   N5 <- N3, N6
#   N6 <- N4, N5
w <- matrix(
  c(  0, 10, 10,  0,  0,  0,
     10,  0,  0, 10,  0,  0,
     10,  0,  0,  0, 10,  0,
      0, 10,  0,  0,  0, 10,
      0,  0, 10,  0,  0, 10,
      0,  0,  0, 10, 10,  0),
  nrow = n, ncol = n, byrow = TRUE
)

# Atividade espontânea de cada neurônio (alpha_i = 0 para todos)
alpha <- rep(0, n)

# Comprimento do burn-in e da amostra final
t_burnin <- 2500
t_amostra <- 40000
t_total   <- t_burnin + t_amostra

# Condições iniciais (4 colunas de inicialização, como nos códigos do João)
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

# Função taxa de disparo (logística padrão)
phi <- function(x) exp(x) / (1 + exp(x))

# Retorna a linha 'indice' da matriz m
matriz_x <- function(m, indice) m[indice, , drop = FALSE]

# Retorna todas as linhas de m EXCETO 'indice'
matriz_pesos <- function(m, indice) m[!(1:nrow(m) %in% indice), , drop = FALSE]

# Índice do último disparo do neurônio 'indice_neuronio' até o instante 'tempo'
ultimo <- function(x, tempo_inicial, tempo_final) {
  u <- tail(which(x[tempo_inicial:tempo_final] == 1), n = 1)
  return(u)
}

# Potencial de membrana do neurônio 'indice_neuronio' no instante 'tempo - 1'
V_novo <- function(X, indice_neuronio, tempo) {
  x      <- matriz_x(X, indice_neuronio)
  indice <- ultimo(x, 1, tempo)

  if (indice == tempo) {
    potencial <- 0
  } else if (indice != (tempo - 1)) {
    adjusted_w <- matriz_pesos(w, indice_neuronio)[,
                    rep(indice_neuronio, tempo - indice)]
    potencial <- sum(
      (1 / (2^(tempo - indice))) *
      matriz_pesos(X, indice_neuronio)[, (indice + 1):tempo] %*% adjusted_w
    )
  } else {
    adjusted_w <- matriz_pesos(w, indice_neuronio)[,
                    rep(indice_neuronio, tempo - indice + 1)]
    potencial <- sum(
      (1 / (2^(tempo - indice))) *
      matriz_pesos(X, indice_neuronio)[, (indice + 1):tempo] %*% adjusted_w
    )
  }
  return(potencial)
}

# Geração da cadeia X (modelo GL com atividade espontânea alpha)
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

# Monta a matriz inicial: 4 colunas fixas + (t_total - 4) colunas a preencher
X_full <- cbind(X_init, matrix(0, nrow = n, ncol = t_total - 4))

# Gera a cadeia completa
X_full <- gera_cadeia(X_full)

# Proporção de disparos por neurônio (para verificar burn-in)
prop_disparos <- rowMeans(X_full)
cat("\nProporção de disparos por neurônio (cadeia completa):\n")
for (i in 1:n) cat(sprintf("  N%d: %.3f\n", i, prop_disparos[i]))

# Descarta burn-in: mantém colunas (t_burnin + 1) até t_total
X <- X_full[, (t_burnin + 1):t_total]

cat(sprintf("\nDimensão da amostra após burn-in: %d neurônios x %d instantes\n",
            nrow(X), ncol(X)))

# Proporção de disparos após burn-in
cat("\nProporção de disparos após burn-in:\n")
for (i in 1:n) cat(sprintf("  N%d: %.3f\n", i, mean(X[i, ])))

# -----------------------------------------------------------------------------
# 4. GRÁFICO DE BURN-IN (proporção acumulada ao longo do tempo)
# -----------------------------------------------------------------------------

# Calcula proporção acumulada em toda a cadeia (para visualização)
proporcao_acumulada <- matrix(0, nrow = n, ncol = ncol(X_full))
for (i in 1:n) {
  proporcao_acumulada[i, ] <- cumsum(X_full[i, ]) / seq_along(X_full[i, ])
}

df_burnin <- as.data.frame(t(proporcao_acumulada))
colnames(df_burnin) <- paste0("N", 1:n)
df_burnin$Tempo <- 1:ncol(X_full)

df_longo <- pivot_longer(df_burnin, cols = -Tempo,
                         names_to = "Neuronio", values_to = "Proporcao")

p_burnin <- ggplot(df_longo, aes(x = Tempo, y = Proporcao, color = Neuronio)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = t_burnin, linetype = "dashed", color = "black") +
  annotate("text", x = t_burnin + 500, y = 0.95,
           label = "Burn-in", hjust = 0, size = 3.5) +
  labs(x = "Tempo (bin)", y = "Proporção acumulada de disparos",
       title = "Verificação do burn-in — rede com 6 neurônios") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("burnin_6neuronios.pdf", plot = p_burnin,
       width = 8, height = 4, device = "pdf")
cat("\nGráfico de burn-in salvo em 'burnin_6neuronios.pdf'\n")

# -----------------------------------------------------------------------------
# 5. CONSTRUÇÃO DAS VARIÁVEIS U_t^{(j -> i)} PARA O mRMR
# -----------------------------------------------------------------------------
# Para cada neurônio pós-sináptico i e cada instante t:
#   - tau_{i,t} = último instante em que i disparou antes de t
#   - L_{i,t}   = t - tau_{i,t} - 1 (comprimento da janela)
#   - U_t^{(j->i)} = soma dos disparos de j no intervalo (tau_{i,t}, t)
#
# Retorna uma lista: U[[i]][[j]] é o vetor de comprimento T com
# as atividades recentes de j em relação a i.

T_obs <- ncol(X)   # número de instantes observados

cat("\nCalculando variáveis U_t^{(j->i)}...\n")

# Pré-calcula, para cada neurônio i, o vetor tau_{i,t} para t = 1,...,T_obs
# (índice do último disparo antes de t; usa 0 se nunca disparou antes)
calcular_tau <- function(X_i, T_obs) {
  tau <- integer(T_obs)
  ultimo_disp <- 0L
  for (t in 1:T_obs) {
    tau[t]     <- ultimo_disp
    if (X_i[t] == 1L) ultimo_disp <- t
  }
  return(tau)
}

# Calcula U_t^{(j->i)} para todos os pares (i, j)
calcular_U <- function(X, n, T_obs) {
  tau_list <- lapply(1:n, function(i) calcular_tau(X[i, ], T_obs))

  U <- vector("list", n)
  for (i in 1:n) {
    U[[i]] <- vector("list", n)
    tau_i  <- tau_list[[i]]
    for (j in (1:n)[-i]) {
      U_ij <- integer(T_obs)
      for (t in 1:T_obs) {
        tau_it <- tau_i[t]
        if (t > 1 && tau_it < t - 1) {
          # Soma disparos de j no intervalo (tau_it, t-1)
          U_ij[t] <- sum(X[j, (tau_it + 1):(t - 1)])
        } else {
          U_ij[t] <- 0L
        }
      }
      U[[i]][[j]] <- U_ij
    }
  }
  return(U)
}

U <- calcular_U(X, n, T_obs)
cat("Variáveis U calculadas.\n")

# -----------------------------------------------------------------------------
# 6. mRMR: RELEVÂNCIA E REDUNDÂNCIA VIA INFORMAÇÃO MÚTUA EMPÍRICA
# -----------------------------------------------------------------------------
# Usa o pacote infotheo para calcular I(U_j; X_i) e I(U_j; X_i | U_S)
# As variáveis são discretizadas (U_ij já é discreta — contagem inteira;
# X_i é binária).

# Função: calcula I(U_j; X_i) para todos os candidatos j != i
calcular_relevancia <- function(U, X, i, n, T_obs) {
  rel <- rep(NA_real_, n)
  X_i <- as.integer(X[i, ])
  for (j in (1:n)[-i]) {
    U_ij <- U[[i]][[j]]
    # discretize: já são inteiros, usamos equalwidth com bins = max+1
    df_ij  <- data.frame(U = as.factor(U_ij), Xi = as.factor(X_i))
    rel[j] <- mutinformation(df_ij$U, df_ij$Xi, method = "emp")
  }
  return(rel)
}

# Função: calcula I(U_j; X_i | U_S) via regra da cadeia:
#   I(U_j; X_i | U_S) = H(X_i | U_S) - H(X_i | U_j, U_S)
# Usa mutinformation do pacote infotheo com condicionamento manual
calcular_relevancia_incremental <- function(U, X, i, j, S, T_obs) {
  X_i  <- as.integer(X[i, ])
  U_j  <- as.integer(U[[i]][[j]])

  if (length(S) == 0) {
    # Sem conjunto selecionado: I(U_j ; X_i)
    return(mutinformation(as.factor(U_j), as.factor(X_i), method = "emp"))
  }

  # Constrói variável conjunta U_S (concatenação dos valores)
  U_S_mat <- do.call(cbind, lapply(S, function(k) as.integer(U[[i]][[k]])))
  # Índice de estado conjunto para U_S
  U_S_joint <- apply(U_S_mat, 1, function(row) paste(row, collapse = "_"))

  # I(U_j ; X_i | U_S) = H(X_i | U_S) - H(X_i | U_j, U_S)
  U_jS_joint <- paste(U_j, U_S_joint, sep = "_")

  H_Xi_dado_US  <- condentropy(as.factor(X_i),   as.factor(U_S_joint),  method = "emp")
  H_Xi_dado_UjS <- condentropy(as.factor(X_i),   as.factor(U_jS_joint), method = "emp")

  inc <- H_Xi_dado_US - H_Xi_dado_UjS
  return(max(inc, 0))   # informação mútua >= 0
}

# Algoritmo mRMR adaptado para redes neuronais
# Parâmetros:
#   lambda  : peso da penalização de redundância
#   M_max   : número máximo de neurônios pré-sinápticos a selecionar
#   tau_q   : quantil para regra de parada por relevância (ex. 0.90)
mrmr_neuronal <- function(U, X, i, n, T_obs,
                           lambda = 1, M_max = n - 1, tau_q = 0.90) {
  candidatos <- (1:n)[-i]

  # Relevâncias marginais
  rel <- calcular_relevancia(U, X, i, n, T_obs)

  # Regra de parada: limiar pelo quantil das relevâncias
  tau_i <- quantile(rel[candidatos], probs = tau_q, na.rm = TRUE)

  S     <- integer(0)   # conjunto selecionado
  ordem <- integer(0)   # ordem de entrada

  for (m in 1:M_max) {
    disponiveis <- setdiff(candidatos, S)
    if (length(disponiveis) == 0) break

    if (m == 1) {
      # Primeiro passo: máxima relevância marginal
      escore  <- rel[disponiveis]
      j_star  <- disponiveis[which.max(escore)]
    } else {
      # Passos seguintes: relevância - lambda * redundância média
      escore <- sapply(disponiveis, function(j) {
        inc_j <- calcular_relevancia_incremental(U, X, i, j, S, T_obs)
        red_j <- rel[j] - inc_j
        rel[j] - lambda * red_j
      })
      j_star <- disponiveis[which.max(escore)]
    }

    # Regra de parada: relevância do próximo candidato abaixo do limiar
    if (rel[j_star] < tau_i && m > 1) break

    S     <- c(S, j_star)
    ordem <- c(ordem, j_star)
  }

  return(list(S = S, ordem = ordem, relevancia = rel))
}

# -----------------------------------------------------------------------------
# 7. APLICAÇÃO DO mRMR PARA TODOS OS 6 NEURÔNIOS
# -----------------------------------------------------------------------------

cat("\nAplicando mRMR para cada neurônio...\n")

resultados <- vector("list", n)
for (i in 1:n) {
  cat(sprintf("  Neurônio %d...\n", i))
  resultados[[i]] <- mrmr_neuronal(
    U, X, i, n, T_obs,
    lambda  = 1,
    M_max   = 2,    # verdadeira vizinhança tem tamanho 2
    tau_q   = 0.90
  )
}

# -----------------------------------------------------------------------------
# 8. RESULTADOS: VIZINHANÇA ESTIMADA vs. VERDADEIRA
# -----------------------------------------------------------------------------

# Vizinhança verdadeira (a partir de W)
vizinhanca_verdadeira <- lapply(1:n, function(i) which(w[i, ] != 0))

cat("\n=== RESULTADOS DO mRMR ===\n")
cat(sprintf("%-10s %-20s %-20s %-10s\n",
            "Neurônio", "Viz. Verdadeira", "Viz. Estimada", "Correto?"))
cat(strrep("-", 62), "\n")

for (i in 1:n) {
  vv  <- sort(vizinhanca_verdadeira[[i]])
  ve  <- sort(resultados[[i]]$S)
  ok  <- ifelse(identical(vv, ve), "SIM", "NÃO")
  cat(sprintf("N%-9d %-20s %-20s %-10s\n",
              i,
              paste(vv, collapse = ", "),
              paste(ve, collapse = ", "),
              ok))
}

# -----------------------------------------------------------------------------
# 9. ESCORES DE RELEVÂNCIA INCREMENTAL (força das arestas)
# -----------------------------------------------------------------------------

cat("\n=== ESCORES DE RELEVÂNCIA INCREMENTAL (W_ij) ===\n")
for (i in 1:n) {
  S_i <- resultados[[i]]$S
  if (length(S_i) == 0) next
  cat(sprintf("\nNeurônio %d (pós-sináptico):\n", i))
  for (j in S_i) {
    S_sem_j <- setdiff(S_i, j)
    w_ij    <- calcular_relevancia_incremental(U, X, i, j, S_sem_j, T_obs)
    cat(sprintf("  N%d -> N%d : W = %.6f\n", j, i, w_ij))
  }
}

# -----------------------------------------------------------------------------
# 10. MATRIZ DE ADJACÊNCIA ESTIMADA
# -----------------------------------------------------------------------------

A_hat <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  for (j in resultados[[i]]$S) {
    A_hat[j, i] <- 1   # aresta j -> i
  }
}

cat("\n=== MATRIZ DE ADJACÊNCIA ESTIMADA (A_hat) ===\n")
cat("Linhas = pré-sináptico (j), Colunas = pós-sináptico (i)\n")
print(A_hat)

cat("\n=== MATRIZ DE ADJACÊNCIA VERDADEIRA (A) ===\n")
A_true <- (w != 0) * 1
print(A_true)

# Métricas de desempenho
VP <- sum(A_hat == 1 & A_true == 1)   # verdadeiros positivos
FP <- sum(A_hat == 1 & A_true == 0)   # falsos positivos
FN <- sum(A_hat == 0 & A_true == 1)   # falsos negativos
VN <- sum(A_hat == 0 & A_true == 0)   # verdadeiros negativos (excl. diagonal)

precisao    <- VP / (VP + FP)
sensib      <- VP / (VP + FN)
f1          <- 2 * precisao * sensib / (precisao + sensib)

cat(sprintf("\nVP=%d  FP=%d  FN=%d  VN=%d\n", VP, FP, FN, VN))
cat(sprintf("Precisão    : %.3f\n", precisao))
cat(sprintf("Sensibilidade: %.3f\n", sensib))
cat(sprintf("F1-score    : %.3f\n", f1))

cat("\n=== CONCLUÍDO ===\n")
