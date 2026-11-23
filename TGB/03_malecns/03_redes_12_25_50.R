# ============================================================================
# TCC - DROSOPHILA | COMPARACAO DE 3 REDES x 4 TRANSFORMACOES DOS PESOS
# Cassiavilani / Ferreira - experimento de sensibilidade, novembro de 2026
#
# REDES: N12, N25, N50
# METODOS: original, divisao_maximo, log_maximo, zscore_positivos
# T=10000 apos burn-in=500; surrogates=199 por passo; lambda=1; alfa=0.05
#
# Rodar: Rscript 01_experimento_comparativo.R
# Ou:    Rscript 01_experimento_comparativo.R N12,N25,N50
# Ou:    Rscript 01_experimento_comparativo.R N25
#
# O peso so entra na simulacao GL. O estimador de informacao e o criterio de
# selecao sobre os disparos sao identicos ao piloto em sentido matematico.
# O codigo usa um agrupamento mais eficiente, validado contra a formula
# plug-in original antes do inicio de cada execucao.
# ============================================================================

VERSAO <- "comparacao_N12_N25_N50_4pesos_v1_20261104"
BASE <- "/home/murilo/tcc_malecns/resultados"
SAIDA <- file.path(BASE, "comparacao_tres_matrizes_pesos_T10000_B199")

REDES <- c("N12", "N25", "N50")
METODOS <- c("original", "divisao_maximo", "log_maximo", "zscore_positivos")
T_AMOSTRA <- 10000L
BURN_IN <- 500L
N_SURROGATES <- 199L
LAMBDA <- 1
ALFA <- 0.05
SEMENTE_CADEIA <- 2026L
# Replicamos a convencao das sementes empregada em cada experimento anterior:
SEMENTE_SURROGATES <- c(N12 = 3026L, N25 = 170000L, N50 = 170000L)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) REDES <- trimws(strsplit(args[[1L]], ",", fixed = TRUE)[[1L]])
if (!length(REDES) || anyDuplicated(REDES) ||
    !all(REDES %in% c("N12", "N25", "N50"))) {
  stop("Redes validas: N12,N25,N50 (ou uma delas).")
}

# E possivel ajustar os caminhos abaixo conforme o servidor, sem modificar
# nenhum parametro estatistico.
ARQUIVOS_W <- list(
  N12 = file.path(BASE, "W_anatomica_original.rds"),
  N25 = file.path(BASE, "expansao_anatomica", "W_subrede_N25.rds"),
  N50 = file.path(BASE, "expansao_anatomica", "W_subrede_N50.rds")
)

pasta_script <- function() {
  arq <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arq)) return(dirname(normalizePath(sub("^--file=", "", arq[1L]))))
  frames <- sys.frames()
  for (k in rev(seq_along(frames))) {
    if (!is.null(frames[[k]]$ofile))
      return(dirname(normalizePath(frames[[k]]$ofile)))
  }
  getwd()
}
DIR_SCRIPT <- pasta_script()

informar <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")),
                              sprintf(...), "\n", sep = "")
salvar_rds <- function(obj, caminho) {
  dir.create(dirname(caminho), recursive = TRUE, showWarnings = FALSE)
  temporario <- tempfile(pattern = ".checkpoint_", tmpdir = dirname(caminho))
  on.exit(if (file.exists(temporario)) unlink(temporario), add = TRUE)
  saveRDS(obj, temporario)
  if (.Platform$OS.type == "windows" && file.exists(caminho))
    stop("Destino ja existe no Windows, evitando sobrescrita: ", caminho)
  if (!file.rename(temporario, caminho)) stop("Nao foi possivel salvar: ", caminho)
}
escrever_csv <- function(tab, arquivo) {
  dir.create(dirname(arquivo), recursive = TRUE, showWarnings = FALSE)
  write.csv(tab, arquivo, row.names = FALSE, na = "")
}

carregar_W <- function(rede) {
  caminho <- ARQUIVOS_W[[rede]]
  if (rede == "N12" && !file.exists(caminho)) {
    ref <- file.path(DIR_SCRIPT, "referencias", "configuracao_original_teste.rds")
    if (file.exists(ref)) {
      obj <- readRDS(ref)
      if (!is.list(obj) || is.null(obj$W_original))
        stop("Referencia N12 sem matriz W_original.")
      W <- as.matrix(obj$W_original)
      caminho <- ref
    }
  }
  if (!exists("W", inherits = FALSE)) {
    if (!file.exists(caminho)) {
      stop("Matriz ", rede, " nao encontrada em ", caminho,
           ". Ajuste ARQUIVOS_W no inicio do script.")
    }
    W <- as.matrix(readRDS(caminho))
  }
  storage.mode(W) <- "double"
  n <- as.integer(sub("N", "", rede, fixed = TRUE))
  stopifnot(is.matrix(W), identical(dim(W), c(n, n)),
            !is.null(rownames(W)), !is.null(colnames(W)),
            identical(rownames(W), colnames(W)),
            !anyDuplicated(rownames(W)),
            all(is.finite(W)), all(W >= 0), all(diag(W) == 0),
            sum(W > 0) > 0)
  if (rede == "N12") stopifnot(sum(W > 0) == 67L, sum(W) == 4638)
  if (rede == "N50") stopifnot(sum(W > 0) == 787L, sum(W) == 48920)
  informar("%s: matriz %dx%d carregada (%d conexoes) de %s",
           rede, n, n, sum(W > 0), caminho)
  W
}

# Falhar no inicio se uma das tres matrizes estiver ausente: nao desperdicamos
# horas de calculo antes de descobrir um arquivo que faltou.
matrizes <- setNames(lapply(REDES, carregar_W), REDES)

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

# As matrizes sao carregadas antes; so iniciamos apos confirmar as 3 fontes.
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)
TODAS_METRICAS <- list()
TODAS_ARESTAS <- list()
TODOS_PESOS <- list()
TODA_DINAMICA <- list()
TODOS_ALVOS <- list()
TRACOS <- list()

for (rede in REDES) {
  W_original <- matrizes[[rede]]
  n <- nrow(W_original)
  nomes <- rownames(W_original)
  verdade <- W_original > 0
  n_verdadeiras <- sum(verdade)
  dir_rede <- file.path(SAIDA, rede)
  dir.create(dir_rede, recursive = TRUE, showWarnings = FALSE)
  for (metodo in METODOS) {
    informar("INICIO | rede=%s | transformacao=%s", rede, metodo)
    W_sim <- transformar_W(W_original, metodo)
    pasta <- file.path(dir_rede, metodo)
    dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
    cfg <- list(versao = VERSAO, rede = rede, metodo = metodo,
                W_original = W_original, W_sim = W_sim,
                T = T_AMOSTRA, burn_in = BURN_IN, B = N_SURROGATES,
                lambda = LAMBDA, alfa = ALFA,
                semente_cadeia = SEMENTE_CADEIA,
                semente_surrogates = unname(SEMENTE_SURROGATES[[rede]]))
    arq_cfg <- file.path(pasta, "configuracao.rds")
    if (file.exists(arq_cfg)) {
      if (!identical(readRDS(arq_cfg), cfg))
        stop("Checkpoint pertence a outra configuracao: ", pasta)
    } else salvar_rds(cfg, arq_cfg)
    
    arq_cadeia <- file.path(pasta, "cadeia_dinamica.rds")
    if (file.exists(arq_cadeia)) {
      cadeia <- readRDS(arq_cadeia)
      informar("Cadeia reutilizada: %s / %s", rede, metodo)
    } else {
      set.seed(SEMENTE_CADEIA)
      cadeia <- simular_gl(W_sim, T_AMOSTRA, BURN_IN)
      salvar_rds(cadeia, arq_cadeia)
      informar("Cadeia salva: %s / %s", rede, metodo)
    }
    X <- cadeia$X
    stopifnot(identical(dim(X), c(T_AMOSTRA, n)),
              identical(colnames(X), nomes), all(X %in% 0:1),
              identical(dim(cadeia$P), c(T_AMOSTRA, n)),
              all(cadeia$P >= 0 & cadeia$P <= 1))
    
    if (metodo == "original") {
      anterior <- switch(rede,
                         N12 = file.path(BASE, "drosophila_piloto_final", "dados", "cadeia_principal.rds"),
                         N25 = file.path(BASE, "aplicacao_mrmr_N25", "cadeia.rds"),
                         N50 = file.path(BASE, "aplicacao_mrmr_N50", "cadeia.rds"))
      if (file.exists(anterior)) {
        Xant <- as.matrix(readRDS(anterior))
        if (!identical(unname(X), unname(Xant)))
          warning("A cadeia original para ", rede,
                  " nao coincide com a anterior; confira sementes e parametros.")
        else informar("Cadeia original reproduziu o piloto anterior: %s", rede)
      }
    }
    
    diagnostico_pesos <- data.frame(
      rede = rede, metodo = metodo, n_neuronios = n,
      n_conexoes_anatomicas = n_verdadeiras,
      peso_positivo_original_min = min(W_original[verdade]),
      peso_positivo_original_max = max(W_original[verdade]),
      peso_transformado_min = min(W_sim[verdade]),
      peso_transformado_max = max(W_sim[verdade]),
      peso_transformado_medio = mean(W_sim[verdade]),
      n_pesos_negativos = sum(W_sim[verdade] < 0),
      n_pesos_zerados = sum(W_sim[verdade] == 0))
    TODOS_PESOS[[paste(rede, metodo)]] <- diagnostico_pesos
    taxa <- colMeans(X)
    por_neuronio <- data.frame(
      rede = rede, metodo = metodo, alvo = nomes,
      taxa_disparo = as.numeric(taxa),
      potencial_medio = as.numeric(colMeans(cadeia$V)),
      probabilidade_media = as.numeric(colMeans(cadeia$P)),
      fracao_p099 = as.numeric(colMeans(cadeia$P > 0.99)),
      fracao_p0999 = as.numeric(colMeans(cadeia$P > 0.999)),
      fracao_p001 = as.numeric(colMeans(cadeia$P < 0.01)))
    TODA_DINAMICA[[paste(rede, metodo)]] <- por_neuronio
    
    W_est <- matrix(0L, n, n, dimnames = dimnames(W_original))
    passos <- list()
    for (i in seq_len(n)) {
      arq_alvo <- file.path(pasta, sprintf("alvo_%02d.rds", i))
      if (file.exists(arq_alvo)) {
        resultado <- readRDS(arq_alvo)
        if (!identical(resultado$config, cfg) || !identical(resultado$i, i))
          stop("Checkpoint de alvo incompatível: ", arq_alvo)
        selecao <- resultado$selecao
        informar("REUTILIZADO %s/%s alvo %d/%d", rede, metodo, i, n)
      } else {
        informar("SELECAO %s/%s alvo %d/%d: %s", rede, metodo, i, n, nomes[i])
        inicio <- Sys.time()
        set.seed(unname(SEMENTE_SURROGATES[[rede]]) + i)
        selecao <- selecionar_alvo(X, i, N_SURROGATES, ALFA, LAMBDA)
        salvar_rds(list(config = cfg, i = i, selecao = selecao), arq_alvo)
        informar("CONCLUIDO em %.2f min: %d conexoes", as.numeric(
          difftime(Sys.time(), inicio, units = "mins")), length(selecao$pais))
      }
      if (length(selecao$pais)) W_est[i, selecao$pais] <- 1L
      if (nrow(selecao$traco)) {
        tr <- selecao$traco
        tr$rede <- rede
        tr$metodo <- metodo
        tr$alvo <- nomes[i]
        passos[[length(passos) + 1L]] <- tr
      }
    }
    
    pos <- which(row(W_original) != col(W_original), arr.ind = TRUE)
    v <- as.integer(verdade[pos])
    e <- as.integer(W_est[pos])
    TP <- sum(v == 1L & e == 1L)
    FP <- sum(v == 0L & e == 1L)
    FN <- sum(v == 1L & e == 0L)
    TN <- sum(v == 0L & e == 0L)
    stopifnot(TP + FP + FN + TN == n * (n - 1L), TP + FN == n_verdadeiras)
    F1 <- if ((2 * TP + FP + FN) > 0L) 2 * TP / (2 * TP + FP + FN) else NA_real_
    resultado_metricas <- data.frame(
      rede = rede, metodo = metodo, n_neuronios = n,
      T = T_AMOSTRA, burn_in = BURN_IN, B = N_SURROGATES,
      lambda = LAMBDA, alfa = ALFA,
      TP = TP, FP = FP, FN = FN, TN = TN,
      sensibilidade = if (TP + FN) TP / (TP + FN) else NA_real_,
      especificidade = if (TN + FP) TN / (TN + FP) else NA_real_,
      precisao = if (TP + FP) TP / (TP + FP) else NA_real_,
      f1 = F1,
      n_selecionadas = TP + FP,
      alvos_exatos = sum(vapply(seq_len(n), function(i)
        identical(which(W_est[i, ] != 0L), which(verdade[i, ])), logical(1))),
      prop_saturacao_099 = mean(cadeia$P > 0.99),
      prop_saturacao_0999 = mean(cadeia$P > 0.999),
      prop_p001 = mean(cadeia$P < 0.01),
      taxa_disparo_media = mean(X))
    TODAS_METRICAS[[paste(rede, metodo)]] <- resultado_metricas
    
    resultado_arestas <- data.frame(
      rede = rede, metodo = metodo,
      alvo = nomes[pos[, 1]], origem = nomes[pos[, 2]],
      peso_original = W_original[pos], peso_transformado = W_sim[pos],
      verdadeira = v, selecionada = e,
      classe = ifelse(v == 1L & e == 1L, "TP",
                      ifelse(v == 0L & e == 1L, "FP",
                             ifelse(v == 1L, "FN", "TN"))),
      stringsAsFactors = FALSE)
    TODAS_ARESTAS[[paste(rede, metodo)]] <- resultado_arestas
    
    resultado_alvos <- do.call(rbind, lapply(seq_len(n), function(i) {
      idx <- pos[, 1] == i
      vi <- v[idx]
      ei <- e[idx]
      tpi <- sum(vi == 1L & ei == 1L)
      fpi <- sum(vi == 0L & ei == 1L)
      fni <- sum(vi == 1L & ei == 0L)
      tni <- sum(vi == 0L & ei == 0L)
      data.frame(rede = rede, metodo = metodo, alvo = nomes[i],
                 TP = tpi, FP = fpi, FN = fni, TN = tni,
                 sensibilidade = if (tpi + fni) tpi / (tpi + fni) else NA_real_,
                 especificidade = if (tni + fpi) tni / (tni + fpi) else NA_real_,
                 recuperacao_exata = as.integer(fpi == 0 && fni == 0))
    }))
    TODOS_ALVOS[[paste(rede, metodo)]] <- resultado_alvos
    if (length(passos)) TRACOS[[paste(rede, metodo)]] <- do.call(rbind, passos)
    
    escrever_csv(resultado_metricas, file.path(pasta, "metricas.csv"))
    escrever_csv(resultado_arestas, file.path(pasta, "arestas.csv"))
    escrever_csv(resultado_alvos, file.path(pasta, "por_alvo.csv"))
    escrever_csv(por_neuronio, file.path(pasta, "dinamica.csv"))
    escrever_csv(diagnostico_pesos, file.path(pasta, "pesos.csv"))
    if (length(passos)) escrever_csv(do.call(rbind, passos),
                                     file.path(pasta, "passos_mrmr.csv"))
    salvar_rds(W_est, file.path(pasta, "W_estimada.rds"))
    informar("RESULTADO %s %s: TP=%d, FP=%d, FN=%d, TN=%d",
             rede, metodo, TP, FP, FN, TN)
  }
}

# ============================================================================
# CONSOLIDACAO: diferencas emparelhadas por matriz e por par dirigido.
# ============================================================================
concatenar <- function(x) {
  a <- do.call(rbind, x)
  rownames(a) <- NULL
  a
}
metricas <- concatenar(TODAS_METRICAS)
arestas <- concatenar(TODAS_ARESTAS)
pesos <- concatenar(TODOS_PESOS)
dinamica <- concatenar(TODA_DINAMICA)
por_alvo <- concatenar(TODOS_ALVOS)

base <- metricas[metricas$metodo == "original",
                 c("rede", "TP", "FP", "FN", "TN", "sensibilidade",
                   "especificidade", "f1", "prop_saturacao_099")]
names(base)[-1L] <- paste0(names(base)[-1L], "_original")
comparacao <- merge(metricas, base, by = "rede", all.x = TRUE, sort = FALSE)
for (m in c("TP", "FP", "FN", "TN", "sensibilidade", "especificidade",
            "f1", "prop_saturacao_099")) {
  comparacao[[paste0("delta_", m)]] <- comparacao[[m]] -
    comparacao[[paste0(m, "_original")]]
}
comparacao <- comparacao[order(match(comparacao$rede, REDES),
                               match(comparacao$metodo, METODOS)), ]

ref <- arestas[arestas$metodo == "original",
               c("rede", "alvo", "origem", "selecionada")]
names(ref)[4L] <- "selecionada_original"
mudancas <- merge(arestas[arestas$metodo != "original", ], ref,
                  by = c("rede", "alvo", "origem"), all.x = TRUE, sort = FALSE)
mudancas$ganhou_TP <- as.integer(mudancas$verdadeira == 1L &
                                   mudancas$selecionada == 1L &
                                   mudancas$selecionada_original == 0L)
mudancas$perdeu_TP <- as.integer(mudancas$verdadeira == 1L &
                                   mudancas$selecionada == 0L &
                                   mudancas$selecionada_original == 1L)
mudancas$adicionou_FP <- as.integer(mudancas$verdadeira == 0L &
                                      mudancas$selecionada == 1L &
                                      mudancas$selecionada_original == 0L)
mudancas$removeu_FP <- as.integer(mudancas$verdadeira == 0L &
                                    mudancas$selecionada == 0L &
                                    mudancas$selecionada_original == 1L)

pasta_global <- file.path(SAIDA, "comparacao_global")
dir.create(pasta_global, recursive = TRUE, showWarnings = FALSE)
escrever_csv(comparacao, file.path(pasta_global, "comparacao_metricas.csv"))
escrever_csv(arestas, file.path(pasta_global, "comparacao_arestas.csv"))
escrever_csv(pesos, file.path(pasta_global, "comparacao_pesos.csv"))
escrever_csv(dinamica, file.path(pasta_global, "comparacao_dinamica_por_alvo.csv"))
escrever_csv(por_alvo, file.path(pasta_global, "comparacao_recuperacao_por_alvo.csv"))
escrever_csv(mudancas, file.path(pasta_global, "mudancas_vs_original.csv"))
if (length(TRACOS)) escrever_csv(concatenar(TRACOS),
                                 file.path(pasta_global, "historico_passos_mrmr.csv"))

# Um grafico por rede com os mesmos eixos e criterios de comparacao.
for (rede in REDES) {
  tab <- comparacao[comparacao$rede == rede, , drop = FALSE]
  tab <- tab[match(METODOS, tab$metodo), , drop = FALSE]
  png(file.path(pasta_global, paste0("comparacao_", rede, ".png")),
      width = 1500, height = 1050, res = 150)
  par(mfrow = c(2, 1), mar = c(7, 5, 4, 2))
  mat <- rbind(TP = tab$TP, FP = tab$FP, FN = tab$FN)
  barplot(mat, beside = TRUE, names.arg = tab$metodo, las = 2,
          col = c("steelblue", "firebrick", "gray60"),
          ylab = "Numero de arestas", main = paste(rede, "- recuperacao de conexoes"))
  legend("topright", legend = rownames(mat),
         fill = c("steelblue", "firebrick", "gray60"), bty = "n")
  barplot(tab$prop_saturacao_099, names.arg = tab$metodo, las = 2,
          ylim = c(0, 1), col = "steelblue",
          ylab = "Fracao de P > 0,99",
          main = paste(rede, "- saturacao logistica"))
  dev.off()
}

informar("COMPARACAO FINAL CONCLUIDA. Resultados em: %s", pasta_global)
print(comparacao[, c("rede", "metodo", "TP", "FP", "FN", "TN",
                     "sensibilidade", "especificidade", "f1",
                     "prop_saturacao_099", "delta_TP", "delta_FP")],
      row.names = FALSE)
