# =============================================================================
# TCC - DROSOPHILA MALE CNS v1.0
# Experimento controlado: efeito da transformacao dos pesos anatomicos
# sobre a simulacao GL e a selecao mRMR com surrogates.
#
# W[i,j] = peso da conexao j -> i. Os pesos NAO entram no mRMR:
# afetam a cadeia X pela simulacao GL. O algoritmo de selecao e mantido.
# Somente R base. Nenhum arquivo original e sobrescrito.
# =============================================================================

# 0. CONFIGURACAO --------------------------------------------------------------
VERSAO_CODIGO <- "comparacao_W_2026-11-04_v1"
BASE <- "/home/murilo/tcc_malecns/resultados"
ARQUIVO_W <- file.path(BASE, "W_anatomica_original.rds")
PASTA_EXPERIMENTO <- file.path(BASE, "experimento_padronizacao_pesos")

# TRUE: verificacao operacional curta. FALSE: experimento principal.
TESTE_RAPIDO <- TRUE
T_AMOSTRA <- if (TESTE_RAPIDO) 1200L else 10000L
BURN_IN <- 500L
N_SURROGATES <- if (TESTE_RAPIDO) 19L else 199L
ALFAS <- 0.05        # Nivel sugerido pelo orientador
LAMBDA <- 1          # Escore = informacao condicional
SEMENTE_CADEIA <- 2026L
SEMENTE_SURROGATES <- 3026L

# Dois alvos mais investigados. Use ALVOS <- NULL para avaliar todos os 12.
ALVOS <- c("MBON11_L", "PPL103_R")

# O z-score e exploratorio: transforma pesos positivos pequenos em negativos,
# o que NAO demonstra que tais sinapses sejam biologicamente inibitorias.
METODOS <- c("original", "divisao_maximo", "log_maximo", "zscore_positivos")

# 1. MATRIZ ANATOMICA E TRANSFORMACOES -----------------------------------------
if (!file.exists(ARQUIVO_W)) {
  stop("Arquivo nao encontrado: ", ARQUIVO_W,
       "\nCopie este script para o servidor onde esta a matriz .rds.")
}
W_original <- as.matrix(readRDS(ARQUIVO_W))
storage.mode(W_original) <- "double"
stopifnot(identical(dim(W_original), c(12L, 12L)),
          identical(rownames(W_original), colnames(W_original)),
          !anyNA(W_original), all(is.finite(W_original)),
          all(W_original >= 0), all(diag(W_original) == 0),
          sum(W_original > 0) == 67L, sum(W_original) == 4638)

n <- nrow(W_original)
nomes <- rownames(W_original)
if (is.null(ALVOS)) ALVOS <- nomes
stopifnot(all(ALVOS %in% nomes), !anyDuplicated(ALVOS))
indices_alvos <- match(ALVOS, nomes)

transformar_pesos <- function(W, metodo) {
  W2 <- W
  positivo <- W > 0
  if (metodo == "original") {
    return(W2)
  } else if (metodo == "divisao_maximo") {
    W2[positivo] <- W[positivo] / max(W[positivo])
  } else if (metodo == "log_maximo") {
    W2[positivo] <- log1p(W[positivo]) / log1p(max(W[positivo]))
  } else if (metodo == "zscore_positivos") {
    # Centralizacao APENAS das 67 sinapses existentes: zeros seguem zeros.
    # Pesos negativos sao resultados matematicos, nao rotulos biologicos.
    dp <- sd(W[positivo])
    if (!is.finite(dp) || dp == 0) stop("DP dos pesos positivos invalido.")
    W2[positivo] <- (W[positivo] - mean(W[positivo])) / dp
  } else {
    stop("Metodo desconhecido: ", metodo)
  }
  stopifnot(identical(dimnames(W2), dimnames(W)),
            all(W2[!positivo] == 0), all(is.finite(W2)))
  W2
}

# 2. SIMULACAO GL (DINAMICA ORIGINAL, SEM ALTERACAO) ---------------------------
simular_gl <- function(W, tamanho, burn_in) {
  n <- nrow(W)
  total <- tamanho + burn_in
  X <- matrix(0L, total, n)
  V <- matrix(0, tamanho, n)
  P <- matrix(0, tamanho, n)
  memoria <- matrix(0, n, n)
  
  for (t in seq_len(total)) {
    potencial <- rowSums(W * memoria)
    probabilidade <- plogis(potencial)
    if (t > burn_in) {
      tt <- t - burn_in
      V[tt, ] <- potencial
      P[tt, ] <- probabilidade
    }
    X[t, ] <- rbinom(n, size = 1L, prob = probabilidade)
    memoria <- memoria / 2 + matrix(X[t, ], nrow = n, ncol = n,
                                    byrow = TRUE)
    memoria[X[t, ] == 1L, ] <- 0
  }
  X <- X[(burn_in + 1L):total, , drop = FALSE]
  colnames(X) <- colnames(W)
  colnames(V) <- colnames(W)
  colnames(P) <- colnames(W)
  list(X = X, V = V, P = P)
}

# 3. VARIAVEIS U (DEFINICAO ORIGINAL DO PILOTO) -------------------------------
criar_variaveis <- function(X, i) {
  Tobs <- nrow(X)
  ultimo <- c(0L, cummax(ifelse(X[seq_len(Tobs - 1L), i] == 1L,
                                seq_len(Tobs - 1L), 0L)))
  tempos <- which(ultimo > 0L)
  if (!length(tempos)) return(NULL)
  U <- matrix(0L, nrow = length(tempos), ncol = ncol(X))
  for (j in seq_len(ncol(X))) {
    acumulada <- cumsum(X[, j])
    U[, j] <- acumulada[tempos - 1L] - acumulada[ultimo[tempos]]
  }
  list(U = U, y = X[tempos, i], tempos = tempos,
       ultimo = ultimo[tempos])
}

# 4. INFORMACAO MUTUA EMPIRICA (MESMA FORMULA) ---------------------------------
entropia <- function(dados) {
  if (is.null(dim(dados))) dados <- matrix(dados, ncol = 1L)
  padroes <- do.call(paste, c(as.data.frame(dados), sep = "|"))
  frequencias <- as.numeric(table(padroes))
  p <- frequencias / sum(frequencias)
  -sum(p * log(p))
}

informacao <- function(u, y, S) {
  if (ncol(S) == 0L) {
    valor <- entropia(u) + entropia(y) - entropia(cbind(u, y))
  } else {
    valor <- entropia(cbind(u, S)) + entropia(cbind(y, S)) -
      entropia(S) - entropia(cbind(u, y, S))
  }
  max(0, valor)
}

# 5. mRMR COM SURROGATES (REGRA ORIGINAL) --------------------------------------
selecionar_alvo <- function(X, i, B, alfa, lambda = 1) {
  dados <- criar_variaveis(X, i)
  if (is.null(dados)) return(integer(0))
  U <- dados$U
  y <- dados$y
  todos <- setdiff(seq_len(ncol(X)), i)
  vazio <- matrix(nrow = length(y), ncol = 0L)
  
  # Relevancia marginal so participa da nota se lambda < 1.
  relevancia <- vapply(todos, function(j) informacao(U[, j], y, vazio),
                       numeric(1))
  names(relevancia) <- as.character(todos)
  selecionados <- integer(0)
  
  repeat {
    candidatos <- setdiff(todos, selecionados)
    if (!length(candidatos)) break
    S <- U[, selecionados, drop = FALSE]
    observado <- vapply(candidatos,
                        function(j) informacao(U[, j], y, S), numeric(1))
    names(observado) <- as.character(candidatos)
    
    maximos_nulos <- numeric(B)
    for (b in seq_len(B)) {
      valores_nulos <- numeric(length(candidatos))
      for (k in seq_along(candidatos)) {
        j <- candidatos[k]
        serie <- X[, j]
        deslocamento <- sample.int(nrow(X) - 1L, size = 1L)
        surrogate <- c(tail(serie, deslocamento),
                       head(serie, nrow(X) - deslocamento))
        acumulada <- cumsum(surrogate)
        u_surrogate <- acumulada[dados$tempos - 1L] -
          acumulada[dados$ultimo]
        valores_nulos[k] <- informacao(u_surrogate, y, S)
      }
      maximos_nulos[b] <- max(valores_nulos)
    }
    limiar <- as.numeric(quantile(maximos_nulos, probs = 1 - alfa,
                                  type = 1))
    aprovados <- candidatos[observado > limiar]
    if (!length(aprovados)) break
    notas <- (1 - lambda) * relevancia[as.character(aprovados)] +
      lambda * observado[as.character(aprovados)]
    selecionados <- c(selecionados, aprovados[which.max(notas)])
  }
  selecionados
}

# 6. TESTES RAPIDOS DAS FUNCOES ------------------------------------------------
X_teste <- cbind(c(1, 0, 0, 1, 0, 0), c(0, 1, 1, 0, 1, 0))
stopifnot(identical(as.integer(criar_variaveis(X_teste, 1)$U[, 2]),
                    c(0L, 1L, 2L, 0L, 1L)))
z <- c(0, 1, 0, 1)
stopifnot(abs(informacao(z, z, matrix(nrow = 4, ncol = 0L)) - log(2)) < 1e-10)
stopifnot(all(transformar_pesos(W_original, "divisao_maximo")[W_original == 0] == 0))

# 7. EXECUCAO REPRODUTIVEL ----------------------------------------------------
modo <- sprintf("%s_T%d_B%d", if (TESTE_RAPIDO) "teste" else "principal",
                T_AMOSTRA, N_SURROGATES)
PASTA_RUN <- file.path(PASTA_EXPERIMENTO, modo)
dir.create(PASTA_RUN, recursive = TRUE, showWarnings = FALSE)

resumos_dinamica <- list()
resumos_metricas <- list()
resumos_arestas <- list()
resumos_pesos <- list()

for (metodo in METODOS) {
  cat("\n==================== ", metodo, " ====================\n", sep = "")
  W_sim <- transformar_pesos(W_original, metodo)
  pasta_metodo <- file.path(PASTA_RUN, metodo)
  dir.create(pasta_metodo, recursive = TRUE, showWarnings = FALSE)
  
  cfg <- list(versao_codigo = VERSAO_CODIGO,
              W_original = W_original, W_sim = W_sim,
              metodo = metodo, T = T_AMOSTRA, burn_in = BURN_IN,
              B = N_SURROGATES, alfas = ALFAS, lambda = LAMBDA,
              alvos = ALVOS, semente_cadeia = SEMENTE_CADEIA,
              semente_surrogates = SEMENTE_SURROGATES)
  arquivo_cfg <- file.path(pasta_metodo, "configuracao.rds")
  if (file.exists(arquivo_cfg) && !identical(readRDS(arquivo_cfg), cfg)) {
    stop("Configuracao alterada na pasta ", pasta_metodo,
         ". Use outra pasta ou renomeie o experimento antigo.")
  }
  if (!file.exists(arquivo_cfg)) saveRDS(cfg, arquivo_cfg)
  
  arq_cadeia <- file.path(pasta_metodo, "cadeia_dinamica.rds")
  if (file.exists(arq_cadeia)) {
    dados <- readRDS(arq_cadeia)
    cat("Cadeia recuperada.\n")
  } else {
    # Mesma semente e mesma sequencia de numeros aleatorios nos metodos.
    set.seed(SEMENTE_CADEIA)
    dados <- simular_gl(W_sim, T_AMOSTRA, BURN_IN)
    saveRDS(dados, arq_cadeia)
  }
  X <- dados$X
  stopifnot(identical(dim(X), c(T_AMOSTRA, n)),
            identical(colnames(X), nomes), all(X %in% 0:1))
  
  # Se existir o piloto original, exigir reproducao IDENTICA da cadeia base.
  arq_piloto <- file.path(BASE, "drosophila_piloto_final", "dados",
                          "cadeia_principal.rds")
  if (metodo == "original" && !TESTE_RAPIDO && file.exists(arq_piloto)) {
    X_piloto <- as.matrix(readRDS(arq_piloto))
    if (!identical(unname(X), unname(X_piloto))) {
      stop("Cadeia original nao coincide com a do piloto; revise o codigo.")
    }
    cat("Cadeia original reproduzida exatamente.\n")
  }
  
  # Guardar pesos e diagnostico de saturacao de CADA alvo avaliado.
  resumos_pesos[[metodo]] <- data.frame(
    metodo = metodo, n_sinapses = sum(W_original > 0),
    peso_minimo = min(W_sim[W_original > 0]),
    peso_maximo = max(W_sim[W_original > 0]),
    peso_medio = mean(W_sim[W_original > 0]),
    n_pesos_negativos = sum(W_sim[W_original > 0] < 0),
    n_pesos_zerados = sum(W_sim[W_original > 0] == 0)
  )
  resumos_dinamica[[metodo]] <- do.call(rbind, lapply(indices_alvos, function(i) {
    data.frame(metodo = metodo, alvo = nomes[i],
               taxa_disparo = mean(X[, i]),
               potencial_medio = mean(dados$V[, i]),
               probabilidade_media = mean(dados$P[, i]),
               proporcao_p_099 = mean(dados$P[, i] > 0.99),
               proporcao_p_0999 = mean(dados$P[, i] > 0.999),
               proporcao_p_001 = mean(dados$P[, i] < 0.01))
  }))
  
  W_est <- matrix(NA_integer_, n, n, dimnames = dimnames(W_original))
  estimados <- vector("list", length(indices_alvos))
  names(estimados) <- ALVOS
  for (k in seq_along(indices_alvos)) {
    i <- indices_alvos[k]
    arq_alvo <- file.path(pasta_metodo, sprintf("alvo_%02d.rds", i))
    if (file.exists(arq_alvo)) {
      selecionados <- readRDS(arq_alvo)
      cat(sprintf("%s: selecao recuperada (%d pais).\n",
                  nomes[i], length(selecionados)))
    } else {
      cat(sprintf("%s: executando mRMR com %d surrogates...\n",
                  nomes[i], N_SURROGATES))
      set.seed(SEMENTE_SURROGATES + i)
      # alpha fixo em 0.05 e lambda=1. Nenhum peso entra aqui.
      selecionados <- selecionar_alvo(X, i, B = N_SURROGATES,
                                      alfa = ALFAS, lambda = LAMBDA)
      saveRDS(selecionados, arq_alvo)
    }
    estimados[[k]] <- selecionados
    W_est[i, ] <- 0L
    if (length(selecionados)) W_est[i, selecionados] <- 1L
  }
  write.csv(W_est, file.path(pasta_metodo, "W_estimada_alvos.csv"))
  saveRDS(estimados, file.path(pasta_metodo, "pais_estimados.rds"))
  
  verdade <- W_original > 0
  idx <- matrix(FALSE, n, n)
  idx[indices_alvos, ] <- TRUE
  diag(idx) <- FALSE
  TP <- sum(W_est[idx] == 1L & verdade[idx])
  FP <- sum(W_est[idx] == 1L & !verdade[idx])
  FN <- sum(W_est[idx] == 0L & verdade[idx])
  TN <- sum(W_est[idx] == 0L & !verdade[idx])
  stopifnot(TP + FP + FN + TN == length(indices_alvos) * (n - 1L))
  
  resumos_metricas[[metodo]] <- data.frame(
    metodo = metodo, alfa = ALFAS, lambda = LAMBDA,
    n_alvos = length(ALVOS), n_pares_avaliados = sum(idx),
    T = T_AMOSTRA, B = N_SURROGATES, TP = TP, FP = FP, FN = FN, TN = TN,
    sensibilidade = if ((TP + FN) > 0) TP / (TP + FN) else NA_real_,
    especificidade = if ((TN + FP) > 0) TN / (TN + FP) else NA_real_,
    precisao = if ((TP + FP) > 0) TP / (TP + FP) else NA_real_,
    n_selecionadas = TP + FP,
    alvos_exatos = sum(vapply(indices_alvos, function(i)
      identical(which(W_est[i, ] == 1L), which(verdade[i, ])), logical(1)))
  )
  
  pares <- which(idx, arr.ind = TRUE)
  tab <- data.frame(metodo = metodo, alfa = ALFAS,
                    alvo = nomes[pares[, 1]], origem = nomes[pares[, 2]],
                    peso_original = W_original[pares],
                    peso_transformado = W_sim[pares],
                    verdadeira = as.integer(verdade[pares]),
                    selecionada = as.integer(W_est[pares]))
  tab$classe <- ifelse(tab$verdadeira == 1L & tab$selecionada == 1L, "TP",
                       ifelse(tab$verdadeira == 0L & tab$selecionada == 1L, "FP",
                              ifelse(tab$verdadeira == 1L, "FN", "TN")))
  resumos_arestas[[metodo]] <- tab
  write.csv(tab, file.path(pasta_metodo, "arestas.csv"), row.names = FALSE)
  cat(sprintf("TP=%d | FP=%d | FN=%d | TN=%d\n", TP, FP, FN, TN))
}

# 8. CONSOLIDACAO E GRAFICOS --------------------------------------------------
combinar <- function(lista) { tab <- do.call(rbind, lista); rownames(tab) <- NULL; tab }
pesos <- combinar(resumos_pesos)
dinamica <- combinar(resumos_dinamica)
metricas <- combinar(resumos_metricas)
arestas <- combinar(resumos_arestas)

baseline <- metricas[metricas$metodo == "original", , drop = FALSE]
stopifnot(nrow(baseline) == 1L)
metricas$delta_TP_vs_original <- metricas$TP - baseline$TP
metricas$delta_FP_vs_original <- metricas$FP - baseline$FP
metricas$delta_FN_vs_original <- metricas$FN - baseline$FN

# Indicadores por alvo evitam que melhorias pontuais se percam no agregado.
resumo_por_alvo <- do.call(rbind, lapply(METODOS, function(m) {
  do.call(rbind, lapply(ALVOS, function(a) {
    tab <- arestas[arestas$metodo == m & arestas$alvo == a, ]
    din <- dinamica[dinamica$metodo == m & dinamica$alvo == a, ]
    TP_a <- sum(tab$classe == "TP")
    FP_a <- sum(tab$classe == "FP")
    FN_a <- sum(tab$classe == "FN")
    TN_a <- sum(tab$classe == "TN")
    data.frame(metodo = m, alvo = a,
               TP = TP_a, FP = FP_a, FN = FN_a, TN = TN_a,
               sensibilidade = if (TP_a + FN_a > 0) TP_a / (TP_a + FN_a)
               else NA_real_,
               especificidade = if (TN_a + FP_a > 0) TN_a / (TN_a + FP_a)
               else NA_real_,
               recuperacao_exata = as.integer(FP_a == 0L && FN_a == 0L),
               taxa_disparo = din$taxa_disparo,
               proporcao_p_099 = din$proporcao_p_099)
  }))
}))
rownames(resumo_por_alvo) <- NULL

write.csv(resumo_por_alvo, file.path(PASTA_RUN, "comparacao_por_alvo.csv"),
          row.names = FALSE)
write.csv(pesos, file.path(PASTA_RUN, "comparacao_pesos.csv"), row.names = FALSE)
write.csv(dinamica, file.path(PASTA_RUN, "comparacao_saturacao.csv"), row.names = FALSE)
write.csv(metricas, file.path(PASTA_RUN, "comparacao_metricas.csv"), row.names = FALSE)
write.csv(arestas, file.path(PASTA_RUN, "comparacao_arestas.csv"), row.names = FALSE)

# Mudancas de decisao versus o metodo original, considerando as mesmas arestas.
base <- arestas[arestas$metodo == "original",
                c("alvo", "origem", "selecionada")]
names(base)[3] <- "selecionada_original"
mais <- merge(arestas[arestas$metodo != "original", ], base,
              by = c("alvo", "origem"), all.x = TRUE, sort = FALSE)
mais$mudou <- as.integer(mais$selecionada != mais$selecionada_original)
mais$ganho_verdadeiro <- as.integer(mais$verdadeira == 1L &
                                      mais$selecionada_original == 0L &
                                      mais$selecionada == 1L)
mais$perda_verdadeiro <- as.integer(mais$verdadeira == 1L &
                                      mais$selecionada_original == 1L &
                                      mais$selecionada == 0L)
write.csv(mais, file.path(PASTA_RUN, "mudancas_vs_original.csv"),
          row.names = FALSE)

# PDF apenas como saida grafica do script R no computador do usuario.
pdf(file.path(PASTA_RUN, "comparacao_graficos.pdf"), width = 11, height = 8)
par(mfrow = c(2, 1), mar = c(5, 4, 3, 1))
barplot(rbind(metricas$TP, metricas$FP, metricas$FN), beside = TRUE,
        names.arg = metricas$metodo, las = 2,
        col = c("steelblue", "firebrick", "gray65"),
        ylab = "Numero de arestas", main = "Acertos e erros por metodo")
legend("topright", legend = c("TP", "FP", "FN"),
       fill = c("steelblue", "firebrick", "gray65"), bty = "n")
barplot(dinamica$proporcao_p_099, names.arg = paste(dinamica$metodo,
                                                    dinamica$alvo, sep = "\n"), las = 2, cex.names = 0.7,
        col = "steelblue", ylim = c(0, 1),
        ylab = "Fracao com P > 0.99", main = "Saturacao logistica por alvo")
dev.off()

cat("\n=== COMPARACAO FINAL ===\n")
print(metricas, row.names = FALSE)
cat("\n=== SATURACAO ===\n")
print(dinamica, row.names = FALSE)
cat("\nResultados em: ", PASTA_RUN, "\n", sep = "")
