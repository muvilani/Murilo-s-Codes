# =============================================================================
# TCC — reconstrução integral de uma ROI do MaleCNS v1.0
# mRMR original: MI condicional plug-in, surrogates circulares e max-T.
# W[alvo, origem] representa origem -> alvo.
# Varia os pesos usados NO GL e alpha usado NA PARADA; lambda fica em 1.
# A execução salva checkpoints por alvo e pode ser retomada com source().
# =============================================================================

# 1. CONFIGURAÇÃO -------------------------------------------------------------
REGIOES <- c("LOP_L_col_09_28")  # Principal: N=64, 268 arestas
# Para comparar depois: c("LOP_L_col_09_28", "LOP_R_col_02_05")
METODOS <- c("zscore_positivos", "log_maximo")
ALFAS <- c(0.01, 0.05, 0.10)
LAMBDA <- 1
T_OBS <- 40000L
BURN_IN <- 500L
B <- 199L
MODO <- "completo"            # "piloto" = só 3 alvos; "completo" = todos
SEED_GL <- 2026L
SEED_SUR <- 3026L
PASTA_W <- file.path(getwd(), "levantamento_regioes", "matrizes")
SAIDA <- file.path(getwd(), "resultados_mrmr_roi")
VERSAO <- "roi_multialfa_lambda1_v1"
stopifnot(LAMBDA == 1, B >= 19, T_OBS > 2,
          all(ALFAS > 0 & ALFAS < 1), !anyDuplicated(ALFAS),
          MODO %in% c("piloto", "completo"))
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)

salvar_seguro <- function(obj, caminho) {
  tmp <- paste0(caminho, ".tmp")
  saveRDS(obj, tmp)
  if (file.exists(caminho)) unlink(caminho)
  if (!file.rename(tmp, caminho)) stop("Não foi possível salvar: ", caminho)
}

# 2. PESOS PARA O GL ----------------------------------------------------------
transformar_pesos <- function(W, metodo) {
  positivo <- W > 0
  if (!any(positivo)) stop("Rede sem arestas.")
  Z <- W
  if (metodo == "log_maximo") {
    Z[positivo] <- log1p(W[positivo]) / log1p(max(W[positivo]))
  } else if (metodo == "zscore_positivos") {
    d <- sd(W[positivo])
    if (!is.finite(d) || d == 0) stop("Desvio padrão de pesos nulo.")
    Z[positivo] <- (W[positivo] - mean(W[positivo])) / d
  } else stop("Transformação desconhecida.")
  stopifnot(all(Z[!positivo] == 0), all(is.finite(Z)))
  Z
}

# 3. SIMULADOR GL UTILIZADO NO EXPERIMENTO ANTERIOR -----------------------------
simular_gl <- function(W, tamanho, burn_in) {
  n <- nrow(W); total <- tamanho + burn_in
  X <- matrix(0L, total, n)
  P <- matrix(0, tamanho, n)
  memoria <- matrix(0, n, n)
  for (t in seq_len(total)) {
    potencial <- rowSums(W * memoria)
    probabilidade <- plogis(potencial)
    if (t > burn_in) P[t - burn_in, ] <- probabilidade
    X[t, ] <- rbinom(n, 1L, probabilidade)
    memoria <- memoria / 2 + matrix(X[t, ], n, n, byrow = TRUE)
    memoria[X[t, ] == 1L, ] <- 0
  }
  X <- X[(burn_in + 1L):total, , drop = FALSE]
  colnames(X) <- colnames(W)
  colnames(P) <- colnames(W)
  list(X = X, P = P)
}

# 4. MESMAS ESTATÍSTICAS PLUG-IN (IMPLEMENTAÇÃO RÁPIDA EQUIVALENTE) ------------
criar_variaveis <- function(X, i) {
  tt <- nrow(X)
  ultimo <- c(0L, cummax(ifelse(X[seq_len(tt - 1L), i] == 1L,
                                seq_len(tt - 1L), 0L)))
  tempos <- which(ultimo > 0L)
  if (!length(tempos)) return(NULL)
  U <- matrix(0L, length(tempos), ncol(X))
  for (j in seq_len(ncol(X))) {
    cs <- cumsum(X[, j])
    U[, j] <- cs[tempos - 1L] - cs[ultimo[tempos]]
  }
  list(U = U, y = X[tempos, i], tempos = tempos,
       ultimo = ultimo[tempos])
}

entropia_rotulos <- function(ids) {
  cont <- as.double(tabulate(match(ids, unique(ids))))
  p <- cont / sum(cont)
  -sum(p * log(p))
}

codigos_condicionantes <- function(U, S) {
  cod <- rep.int(1L, nrow(U)); k <- 1L
  for (j in S) {
    uj <- match(U[, j], unique(U[, j]))
    chave <- as.double(cod) + (as.double(uj) - 1) * as.double(k)
    cod <- match(chave, unique(chave)); k <- max(cod)
  }
  cod
}

preparar_estatistica <- function(U, y, S) {
  id_s <- codigos_condicionantes(U, S)
  ycod <- match(y, unique(y))
  chave <- as.double(id_s) + (as.double(ycod) - 1) * max(id_s)
  id_ys <- match(chave, unique(chave))
  list(id_s = id_s, id_ys = id_ys,
       ks = max(id_s), kys = max(id_ys),
       Hs = entropia_rotulos(id_s), Hys = entropia_rotulos(id_ys))
}

informacao_rapida <- function(u, prep) {
  uc <- match(u, unique(u))
  id_us <- as.double(prep$id_s) + (as.double(uc) - 1) * prep$ks
  id_uys <- as.double(prep$id_ys) + (as.double(uc) - 1) * prep$kys
  max(0, entropia_rotulos(id_us) + prep$Hys - prep$Hs -
        entropia_rotulos(id_uys))
}

# Teste independente da matriz real: estimativa exata para u=y binários.
Ut <- matrix(c(0, 1, 0, 1), ncol = 1L)
pt <- preparar_estatistica(Ut, c(0, 1, 0, 1), integer(0))
stopifnot(abs(informacao_rapida(Ut[, 1], pt) - log(2)) < 1e-10)

# 5. MAX-T MULTIALFA ---------------------------------------------------------
# Para lambda=1, entre candidatos aprovados é escolhido aquele com MAIOR CMI.
# O candidato escolhido é o mesmo para qualquer alpha ainda ativo. Assim,
# os três alpha compartilham os deslocamentos e um único caminho progressivo.
selecionar_multialfa <- function(X, i, B, alfas) {
  dados <- criar_variaveis(X, i)
  vazio <- setNames(replicate(length(alfas), integer(0), simplify = FALSE),
                    as.character(alfas))
  if (is.null(dados)) return(list(pais = vazio, traco = data.frame()))
  U <- dados$U; y <- dados$y
  todos <- setdiff(seq_len(ncol(X)), i)
  escolhidos <- integer(0)
  ativo <- rep(TRUE, length(alfas))
  nomes_a <- as.character(alfas)
  traco <- list(); passo <- 0L
  repeat {
    candidatos <- setdiff(todos, escolhidos)
    if (!length(candidatos) || !any(ativo)) break
    passo <- passo + 1L
    prep <- preparar_estatistica(U, y, escolhidos)
    obs <- vapply(candidatos,
                  function(j) informacao_rapida(U[, j], prep), numeric(1))
    melhor <- which.max(obs); maior_obs <- obs[melhor]
    maximos <- numeric(B)
    for (b in seq_len(B)) {
      valores <- numeric(length(candidatos))
      for (k in seq_along(candidatos)) {
        j <- candidatos[k]
        serie <- X[, j]
        shift <- sample.int(nrow(X) - 1L, 1L)
        serie_s <- c(tail(serie, shift), head(serie, nrow(X) - shift))
        cs <- cumsum(serie_s)
        us <- cs[dados$tempos - 1L] - cs[dados$ultimo]
        valores[k] <- informacao_rapida(us, prep)
      }
      maximos[b] <- max(valores)
    }
    limiares <- vapply(alfas,
                       function(a) as.numeric(quantile(maximos, 1 - a, type = 1)), numeric(1))
    aprova <- ativo & (maior_obs > limiares)
    traco[[passo]] <- data.frame(
      passo = passo, melhor_indice = candidatos[melhor],
      maior_cmi = maior_obs,
      p_maxT = (1 + sum(maximos >= maior_obs)) / (B + 1),
      alpha = alfas, limiar = limiares, aprovado = aprova
    )
    for (h in which(aprova)) vazio[[h]] <- c(vazio[[h]], candidatos[melhor])
    ativo <- aprova
    if (!any(ativo)) break
    escolhidos <- c(escolhidos, candidatos[melhor])
  }
  # Checagem da propriedade de aninhamento para lambda=1.
  ord <- order(alfas)
  for (h in seq_len(length(ord) - 1L)) {
    menor <- vazio[[ord[h]]]; maior <- vazio[[ord[h + 1L]]]
    stopifnot(identical(menor, head(maior, length(menor))))
  }
  list(pais = vazio, traco = do.call(rbind, traco))
}

# 6. EXECUÇÃO COM CHECKPOINT POR ALVO -----------------------------------------
linhas <- list()
for (roi in REGIOES) {
  nome <- gsub("[^A-Za-z0-9_-]", "_", roi)
  arq_W <- file.path(PASTA_W, nome, "W_regiao.rds")
  if (!file.exists(arq_W)) stop("Matriz não encontrada: ", arq_W)
  W <- as.matrix(readRDS(arq_W))
  storage.mode(W) <- "double"
  stopifnot(nrow(W) == ncol(W),
            identical(rownames(W), colnames(W)),
            all(is.finite(W)), all(W >= 0), all(diag(W) == 0))
  n <- nrow(W); verdade <- W > 0
  cat(sprintf("\nROI=%s | N=%d | arestas=%d\n", roi, n, sum(verdade)))
  k <- rowSums(verdade)
  ALVOS <- if (MODO == "completo") seq_len(n) else
    unique(c(which.min(k), which.min(abs(k - median(k))), which.max(k)))
  for (metodo in METODOS) {
    pasta <- file.path(SAIDA, nome, metodo,
                       paste0("T", T_OBS, "_B", B))
    dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
    Wgl <- transformar_pesos(W, metodo)
    config <- list(versao = VERSAO, roi = roi, metodo = metodo, W = W,
                   T = T_OBS, burn_in = BURN_IN, B = B, alfas = ALFAS,
                   lambda = LAMBDA, seed_gl = SEED_GL,
                   seed_sur = SEED_SUR)
    arq_cfg <- file.path(pasta, "configuracao.rds")
    if (file.exists(arq_cfg)) {
      if (!identical(config, readRDS(arq_cfg)))
        stop("Configuração incompatível com checkpoint em ", pasta)
    } else salvar_seguro(config, arq_cfg)
    arq_gl <- file.path(pasta, "gl_simulado.rds")
    if (file.exists(arq_gl)) {
      gl <- readRDS(arq_gl)
    } else {
      set.seed(SEED_GL)
      gl <- simular_gl(Wgl, T_OBS, BURN_IN)
      salvar_seguro(gl, arq_gl)
    }
    X <- gl$X
    stopifnot(identical(dim(X), c(T_OBS, n)),
              identical(colnames(X), colnames(W)))
    frac_saturacao <- mean(gl$P > 0.99 | gl$P < 0.01)
    cat(sprintf("Método=%s | saturação=%.4f | alvos=%d\n",
                metodo, frac_saturacao, length(ALVOS)))
    for (i in ALVOS) {
      ck <- file.path(pasta, sprintf("alvo_%03d.rds", i))
      if (file.exists(ck)) {
        sel <- readRDS(ck)
      } else {
        t0 <- proc.time()[[3L]]
        cat(sprintf("  Calculando alvo %d/%d (pais=%d)...\n", i, n, k[i]))
        set.seed(SEED_SUR + i)
        sel <- selecionar_multialfa(X, i, B, ALFAS)
        sel$segundos <- proc.time()[[3L]] - t0
        salvar_seguro(sel, ck)
      }
      for (h in seq_along(ALFAS)) {
        pais <- sel$pais[[as.character(ALFAS[h])]]
        v <- setdiff(which(verdade[i, ]), i)
        TP <- length(intersect(v, pais))
        FP <- length(setdiff(pais, v))
        FN <- length(setdiff(v, pais))
        TN <- n - 1L - TP - FP - FN
        linhas[[length(linhas) + 1L]] <- data.frame(
          roi = roi, metodo = metodo, alfa = ALFAS[h],
          alvo = rownames(W)[i], N = n, n_alvos_previstos = length(ALVOS),
          n_pais_verdadeiros = length(v), n_selecionados = length(pais),
          TP = TP, FP = FP, FN = FN, TN = TN,
          exato = as.integer(setequal(v, pais)),
          segundos_alvo = sel$segundos, saturacao = frac_saturacao
        )
      }
      parcial <- do.call(rbind, linhas)
      write.csv(parcial, file.path(SAIDA, "metricas_parciais.csv"),
                row.names = FALSE)
      cat(sprintf("  Salvo alvo %d/%d | %.1f min\n",
                  i, n, sel$segundos / 60))
    }
    rm(gl, X); gc(verbose = FALSE)
  }
}

# 7. RESUMO POR REGIÃO / TRANSFORMAÇÃO / ALPHA -------------------------------
df <- do.call(rbind, linhas)
chaves <- unique(df[, c("roi", "metodo", "alfa")])
res <- do.call(rbind, lapply(seq_len(nrow(chaves)), function(z) {
  s <- chaves[z, ]
  x <- df[df$roi == s$roi & df$metodo == s$metodo & df$alfa == s$alfa, ]
  TP <- sum(x$TP); FP <- sum(x$FP); FN <- sum(x$FN); TN <- sum(x$TN)
  data.frame(roi = s$roi, metodo = s$metodo, alfa = s$alfa,
             N = x$N[1], n_alvos = nrow(x),
             TP = TP, FP = FP, FN = FN, TN = TN,
             sensibilidade = if (TP + FN > 0) TP / (TP + FN) else NA_real_,
             precisao = if (TP + FP > 0) TP / (TP + FP) else NA_real_,
             F1 = if (2 * TP + FP + FN > 0) 2 * TP / (2 * TP + FP + FN)
             else NA_real_,
             fracao_exatos = mean(x$exato),
             saturacao = x$saturacao[1],
             horas_alvos = sum(x$segundos_alvo) / 3600)
}))
rownames(res) <- NULL
write.csv(df, file.path(SAIDA, "metricas_por_alvo.csv"), row.names = FALSE)
write.csv(res, file.path(SAIDA, "resumo_final.csv"), row.names = FALSE)
print(res, row.names = FALSE)
cat("\nArquivos salvos em: ", normalizePath(SAIDA), "\n", sep = "")
