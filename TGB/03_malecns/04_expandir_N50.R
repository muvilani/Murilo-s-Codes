# DROSOPHILA N=50 - T x ALFA x TRANSFORMACAO DOS PESOS
# No RStudio: Ctrl+A, Ctrl+Enter. Sem argumentos e sem outro script.
# Continua automaticamente dos checkpoints do experimento anterior.

Ts <- c(40000L, 80000L)
ALFAS <- c(0.01, 0.05, 0.10)
METODOS <- c("original", "divisao_maximo", "log_maximo", "zscore_positivos")
LAMBDA <- 1
BURN_IN <- 500L
N_SURROGATES <- 199L
SEMENTE_CADEIA <- 2026L
SEMENTE_SURROGATES <- 170000L
VERSAO <- "N50_T_alfa_v2"  # preserva compatibilidade com checkpoints anteriores

# A matriz anatomica original N50 ja esta no servidor.
arquivo_W <- "/home/murilo/tcc_malecns/resultados/expansao_anatomica/W_subrede_N50.rds"
if (!file.exists(arquivo_W)) stop("Matriz N50 ausente do local usado no experimento anterior.")
W_original <- as.matrix(readRDS(arquivo_W))
storage.mode(W_original) <- "double"
stopifnot(identical(dim(W_original), c(50L, 50L)),
          !is.null(rownames(W_original)),
          identical(rownames(W_original), colnames(W_original)),
          all(is.finite(W_original)), all(W_original >= 0),
          all(diag(W_original) == 0), sum(W_original > 0) == 787L)
nomes <- rownames(W_original)
verdade <- W_original > 0
SAIDA <- "resultados_N50_T"  # mesma pasta da execucao anterior
# Se o diretorio de trabalho mudar, abra o projeto no RStudio antes de executar.
dir.create(SAIDA, showWarnings = FALSE, recursive = TRUE)

informar <- function(...) {
  cat(format(Sys.time(), "%F %T"), "|", sprintf(...), "\n")
  flush.console()
}
salvar_seguro <- function(obj, destino) {
  dir.create(dirname(destino), showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(".temp_", tmpdir = dirname(destino))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(obj, tmp)
  if (file.exists(destino)) stop("Arquivo ja existe, nao sobrescrevendo: ", destino)
  if (!file.rename(tmp, destino)) stop("Falha ao salvar: ", destino)
}

# Atualiza tabelas de progresso com os CSVs finalizados, mesmo apos reinicio.
resumir <- function() {
  arquivos <- list.files(SAIDA, pattern = "^metricas\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(arquivos)) return(invisible(NULL))
  tabs <- lapply(arquivos, function(f) tryCatch(read.csv(f), error = function(e) NULL))
  tabs <- Filter(Negate(is.null), tabs)
  if (!length(tabs)) return(invisible(NULL))
  campos <- c("T", "alfa", "lambda", "metodo", "TP", "FP", "FN", "TN",
              "sensibilidade", "especificidade", "precisao", "F1")
  tabs <- lapply(tabs, function(x) {
    faltam <- setdiff(campos, names(x))
    for (nome in faltam) x[[nome]] <- NA
    x[, campos]
  })
  tabela <- do.call(rbind, tabs)
  tabela <- tabela[order(tabela$T, match(tabela$metodo, METODOS), tabela$alfa), ]
  rownames(tabela) <- NULL
  write.csv(tabela, file.path(SAIDA, "resumo_geral.csv"), row.names = FALSE)
  invisible(tabela)
}

transformar_W <- function(W, metodo) {
  positivas <- W > 0
  vals <- W[positivas]
  if (metodo == "original") return(W)
  novo <- W
  if (metodo == "divisao_maximo") {
    novo[positivas] <- vals / max(vals)
  } else if (metodo == "log_maximo") {
    novo[positivas] <- log1p(vals) / log1p(max(vals))
  } else if (metodo == "zscore_positivos") {
    d <- sd(vals)
    if (!is.finite(d) || d <= 0) stop("Z-score impossivel: DP nulo.")
    novo[positivas] <- (vals - mean(vals)) / d
  } else stop("Metodo desconhecido: ", metodo)
  stopifnot(identical(dimnames(novo), dimnames(W)),
            all(is.finite(novo)), all(novo[!positivas] == 0))
  novo
}

# Modelo GL original: W[i,j] representa j -> i; memoria decai por 1/2 e
# reinicia a linha quando o alvo i dispara. Nao ha termo de intercepto.
simular_gl <- function(W, T, burn_in) {
  n <- nrow(W)
  total <- T + burn_in
  X <- matrix(0L, total, n)
  V <- matrix(0, T, n)
  P <- matrix(0, T, n)
  mem <- matrix(0, n, n)
  for (t in seq_len(total)) {
    v <- rowSums(W * mem)
    p <- plogis(v)
    if (t > burn_in) {
      k <- t - burn_in
      V[k, ] <- v
      P[k, ] <- p
    }
    X[t, ] <- rbinom(n, size = 1L, prob = p)
    mem <- mem / 2 + matrix(X[t, ], nrow = n, ncol = n, byrow = TRUE)
    mem[X[t, ] == 1L, ] <- 0
  }
  X <- X[(burn_in + 1L):total, , drop = FALSE]
  colnames(X) <- colnames(W)
  colnames(V) <- colnames(W)
  colnames(P) <- colnames(W)
  list(X = X, V = V, P = P)
}

criar_variaveis <- function(X, i) {
  T <- nrow(X)
  ultimo <- c(0L, cummax(ifelse(X[seq_len(T - 1L), i] == 1L,
                                seq_len(T - 1L), 0L)))
  tempos <- which(ultimo > 0L)
  if (!length(tempos)) return(NULL)
  U <- matrix(0L, length(tempos), ncol(X))
  for (j in seq_len(ncol(X))) {
    cs <- cumsum(X[, j])
    U[, j] <- cs[tempos - 1L] - cs[ultimo[tempos]]
  }
  list(U = U, y = X[tempos, i], tempos = tempos, ultimo = ultimo[tempos])
}

# Entropias plug-in. O codigo rapido usa identificadores inteiros para cada
# configuracao discreta; os agrupamentos sao os mesmos de paste/table originais.
entropia_rotulos <- function(ids) {
  f <- tabulate(match(ids, unique(ids)))
  p <- f / sum(f)
  -sum(p * log(p))
}

codigos_condicionantes <- function(U, S) {
  cod <- rep.int(1L, nrow(U))
  k <- 1L
  for (j in S) {
    uj <- match(U[, j], unique(U[, j]))
    chave <- cod + (uj - 1L) * k
    cod <- match(chave, unique(chave))
    k <- max(cod)
  }
  cod
}

preparar_estatistica <- function(U, y, S) {
  id_s <- codigos_condicionantes(U, S)
  ycod <- match(y, unique(y))
  chave <- id_s + (ycod - 1L) * max(id_s)
  id_ys <- match(chave, unique(chave))
  list(id_s = id_s, id_ys = id_ys, ks = max(id_s), kys = max(id_ys),
       Hs = entropia_rotulos(id_s), Hys = entropia_rotulos(id_ys))
}

informacao_rapida <- function(u, prep) {
  uc <- match(u, unique(u))
  id_us <- prep$id_s + (uc - 1L) * prep$ks
  id_uys <- prep$id_ys + (uc - 1L) * prep$kys
  max(0, entropia_rotulos(id_us) + prep$Hys - prep$Hs -
        entropia_rotulos(id_uys))
}

# Implementacao literalmente equivalente do piloto, para testes de regressao.
entropia_legada <- function(dados) {
  if (is.null(dim(dados))) dados <- matrix(dados, ncol = 1L)
  keys <- do.call(paste, c(as.data.frame(dados), sep = "|"))
  freq <- as.numeric(table(keys))
  p <- freq / sum(freq)
  -sum(p * log(p))
}
informacao_legada <- function(u, y, S) {
  if (ncol(S) == 0L) {
    v <- entropia_legada(u) + entropia_legada(y) -
      entropia_legada(cbind(u, y))
  } else {
    v <- entropia_legada(cbind(u, S)) + entropia_legada(cbind(y, S)) -
      entropia_legada(S) - entropia_legada(cbind(u, y, S))
  }
  max(0, v)
}

# Verifica a equivalencia numerica do estimador em varios condicionamentos,
# incluindo zero, um e multiplos condicionantes, antes da simulacao principal.
testar_estatistica <- function() {
  set.seed(3719L)
  Ut <- matrix(sample(0:5, 500L * 5L, replace = TRUE), 500L, 5L)
  yt <- sample(0:1, 500L, replace = TRUE)
  for (S in list(integer(0), 1L, c(1L, 2L), c(1L, 2L, 3L))) {
    base <- preparar_estatistica(Ut, yt, S)
    for (j in setdiff(seq_len(5L), S)) {
      a <- informacao_rapida(Ut[, j], base)
      b <- informacao_legada(Ut[, j], yt, Ut[, S, drop = FALSE])
      if (!isTRUE(all.equal(a, b, tolerance = 1e-10)))
        stop("Teste de equivalencia informacional falhou: S=",
             paste(S, collapse = ","), " j=", j, " rapido=", a,
             " legado=", b)
    }
  }
  Xtoy <- cbind(c(1, 0, 0, 1, 0, 0), c(0, 1, 1, 0, 1, 0))
  stopifnot(identical(as.integer(criar_variaveis(Xtoy, 1L)$U[, 2L]),
                      c(0L, 1L, 2L, 0L, 1L)))
  informar("Autotestes: estatistica plug-in equivalente ao codigo legado.")
}
testar_estatistica()

# Mesma selecao do original: CMI observada, B surrogates com deslocamento
# circular por candidato, limiar max-T e selecao do maior escore Q_lambda.
selecionar_alvo <- function(X, i, B, alfa, lambda) {
  dados <- criar_variaveis(X, i)
  if (is.null(dados)) return(list(pais = integer(0), traco = data.frame()))
  U <- dados$U
  y <- dados$y
  todos <- setdiff(seq_len(ncol(X)), i)
  prep_vazio <- preparar_estatistica(U, y, integer(0))
  relevancia <- vapply(todos, function(j)
    informacao_rapida(U[, j], prep_vazio), numeric(1))
  names(relevancia) <- as.character(todos)
  selecionados <- integer(0)
  historico <- list()
  repita <- 0L
  repeat {
    candidatos <- setdiff(todos, selecionados)
    if (!length(candidatos)) break
    repita <- repita + 1L
    prep <- preparar_estatistica(U, y, selecionados)
    observado <- vapply(candidatos, function(j)
      informacao_rapida(U[, j], prep), numeric(1))
    names(observado) <- as.character(candidatos)
    maximos <- numeric(B)
    for (b in seq_len(B)) {
      res <- numeric(length(candidatos))
      for (k in seq_along(candidatos)) {
        j <- candidatos[k]
        serie <- X[, j]
        d <- sample.int(nrow(X) - 1L, size = 1L)
        serie_s <- c(tail(serie, d), head(serie, nrow(X) - d))
        cs <- cumsum(serie_s)
        U_s <- cs[dados$tempos - 1L] - cs[dados$ultimo]
        res[k] <- informacao_rapida(U_s, prep)
      }
      maximos[b] <- max(res)
    }
    limiar <- as.numeric(quantile(maximos, probs = 1 - alfa, type = 1))
    aprovados <- candidatos[observado > limiar]
    escolhido <- NA_integer_
    if (length(aprovados)) {
      notas <- (1 - lambda) * relevancia[as.character(aprovados)] +
        lambda * observado[as.character(aprovados)]
      escolhido <- aprovados[which.max(notas)]
    }
    historico[[repita]] <- data.frame(
      passo = repita, n_candidatos = length(candidatos),
      maximo_observado = max(observado), limiar = limiar,
      n_aprovados = length(aprovados),
      escolhido = if (is.na(escolhido)) "" else colnames(X)[escolhido],
      stringsAsFactors = FALSE
    )
    if (!length(aprovados)) break
    selecionados <- c(selecionados, escolhido)
  }
  traco <- if (length(historico)) do.call(rbind, historico) else data.frame()
  list(pais = as.integer(selecionados), traco = traco)
}


# Primeiro o metodo original para recuperar os 6 resultados ja produzidos.
for (metodo in METODOS) {
  W_sim <- transformar_W(W_original, metodo)
  informar("METODO: %s | positivos=%d | negativos=%d | zeros=%d", metodo,
           sum(W_sim[verdade] > 0), sum(W_sim[verdade] < 0), sum(W_sim[verdade] == 0))
  for (T_AMOSTRA in Ts) {
    pasta_T <- file.path(SAIDA, paste0("T", T_AMOSTRA), metodo)
    dir.create(pasta_T, showWarnings = FALSE, recursive = TRUE)
    arq_cadeia <- file.path(pasta_T, "cadeia.rds")
    if (file.exists(arq_cadeia)) {
      cadeia <- readRDS(arq_cadeia)
      informar("T=%d %s: cadeia recuperada", T_AMOSTRA, metodo)
    } else {
      informar("T=%d %s: simulando cadeia", T_AMOSTRA, metodo)
      set.seed(SEMENTE_CADEIA)
      cadeia <- simular_gl(W_sim, T_AMOSTRA, BURN_IN)
      salvar_seguro(cadeia, arq_cadeia)
      informar("T=%d %s: cadeia salva", T_AMOSTRA, metodo)
    }
    X <- cadeia$X
    stopifnot(identical(dim(X), c(T_AMOSTRA, 50L)),
              identical(colnames(X), nomes), all(X %in% c(0L, 1L)))
    # Diagnostico de saturacao do processo de disparos para cada transformacao.
    dinamica <- data.frame(metodo = metodo, T = T_AMOSTRA, alvo = nomes,
                           taxa_disparo = colMeans(X),
                           prop_p099 = colMeans(cadeia$P > 0.99),
                           prop_p0999 = colMeans(cadeia$P > 0.999),
                           media_p = colMeans(cadeia$P))
    write.csv(dinamica, file.path(pasta_T, "dinamica.csv"), row.names = FALSE)
    for (alfa in ALFAS) {
      pasta <- file.path(pasta_T, paste0("alfa_", gsub("\\.", "_", sprintf("%.2f", alfa))))
      dir.create(pasta, showWarnings = FALSE, recursive = TRUE)
      cfg <- list(versao = VERSAO, T = T_AMOSTRA, alfa = alfa, lambda = LAMBDA,
                  B = N_SURROGATES, burn_in = BURN_IN, metodo = metodo,
                  W_original = W_original, W_sim = W_sim,
                  semente_cadeia = SEMENTE_CADEIA,
                  semente_surrogates = SEMENTE_SURROGATES)
      arq_cfg <- file.path(pasta, "configuracao.rds")
      if (file.exists(arq_cfg)) {
        if (!identical(readRDS(arq_cfg), cfg))
          stop("Configuracao diferente dos checkpoints existentes: ", pasta)
      } else salvar_seguro(cfg, arq_cfg)
      
      # Caso a configuracao toda ja esteja concluida, nao recalcula nada.
      arquivos_alvos <- file.path(pasta, sprintf("alvo_%02d.rds", 1:50))
      arq_metricas <- file.path(pasta, "metricas.csv")
      if (all(file.exists(arquivos_alvos)) && file.exists(arq_metricas)) {
        informar("JA CONCLUIDO: %s T=%d alfa=%.2f", metodo, T_AMOSTRA, alfa)
        resumir()
        next
      }
      informar("INICIO: %s | T=%d | alfa=%.2f", metodo, T_AMOSTRA, alfa)
      W_est <- matrix(0L, 50L, 50L, dimnames = dimnames(W_original))
      for (i in seq_len(50L)) {
        arq_alvo <- arquivos_alvos[i]
        if (file.exists(arq_alvo)) {
          obj <- readRDS(arq_alvo)
          stopifnot(identical(obj$config, cfg), identical(obj$i, i))
          selecao <- obj$selecao
          informar("RECUPERADO: %s T=%d alfa=%.2f alvo=%d/50",
                   metodo, T_AMOSTRA, alfa, i)
        } else {
          inicio <- Sys.time()
          informar("CALCULANDO: %s T=%d alfa=%.2f alvo=%d/50",
                   metodo, T_AMOSTRA, alfa, i)
          set.seed(SEMENTE_SURROGATES + i)
          selecao <- selecionar_alvo(X, i, N_SURROGATES, alfa, LAMBDA)
          salvar_seguro(list(config = cfg, i = i, selecao = selecao), arq_alvo)
          informar("SALVO: alvo=%d/50 | %.2f min | %d pais", i,
                   as.numeric(difftime(Sys.time(), inicio, units = "mins")),
                   length(selecao$pais))
        }
        if (length(selecao$pais)) W_est[i, selecao$pais] <- 1L
      }
      idx <- row(W_est) != col(W_est)
      v <- verdade[idx]
      e <- W_est[idx] != 0L
      TP <- sum(v & e); FP <- sum(!v & e)
      FN <- sum(v & !e); TN <- sum(!v & !e)
      stopifnot(TP + FP + FN + TN == 2450, TP + FN == 787)
      met <- data.frame(T = T_AMOSTRA, alfa = alfa, lambda = LAMBDA,
                        metodo = metodo, TP = TP, FP = FP, FN = FN, TN = TN,
                        sensibilidade = TP / (TP + FN),
                        especificidade = TN / (TN + FP),
                        precisao = if (TP + FP > 0) TP / (TP + FP) else NA_real_,
                        F1 = if (2 * TP + FP + FN > 0) 2 * TP / (2 * TP + FP + FN) else NA_real_)
      # Primeiro salva a rede, depois as metricas: sem falsos resultados completos.
      saveRDS(W_est, file.path(pasta, "W_estimada.rds"))
      write.csv(met, arq_metricas, row.names = FALSE)
      informar("CONCLUIDO: %s T=%d alfa=%.2f TP=%d FP=%d", metodo, T_AMOSTRA,
               alfa, TP, FP)
      resumir()
    }
    rm(cadeia, X)
    gc(verbose = FALSE)
  }
}
resumo <- resumir()
informar("EXPERIMENTO FINALIZADO - %d configuracoes em resumo_geral.csv", nrow(resumo))
print(resumo, row.names = FALSE)
