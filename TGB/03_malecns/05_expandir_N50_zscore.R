# =============================================================================
# TCC | EXPANSAO Z-SCORE N=50 | T=40000 e T=80000
# Execucao pelo RStudio: abra o arquivo, Ctrl+A, Ctrl+Enter (ou Source).
# Execute na maquina virtual que possui /home/murilo/tcc_malecns.
# Preserva resultados historicos; checkpoints individuais por alvo.
# Nao requer pacotes adicionais (somente base R).
# =============================================================================

PROJETO <- "/home/murilo/tcc_malecns"
SCRIPT_REFERENCIA <- file.path(PROJETO, "TCC_MC_N12_N25_N50_N100", "02_monte_carlo.R")
MATRIZ_REFERENCIA <- file.path(PROJETO, "resultados", "expansao_anatomica", "W_subrede_N50.rds")
CONFIG_HISTORICA <- file.path(PROJETO, "resultados", "comparacao_tres_matrizes_pesos_T10000_B199", "N50", "zscore_positivos", "configuracao.rds")
SAIDA <- file.path(PROJETO, "resultados", "zscore_N50_T40000_T80000_v1")

TS <- c(40000L, 80000L)
ALFA <- 0.05
LAMBDA <- 1
B <- 199L
BURN_IN <- 500L
SEMENTE_CADEIA <- 2026L
SEMENTE_SURROGATES <- 170000L
VERSAO <- "zscore_N50_expansao_T_v1_20261111"

# Altere para TRUE apenas se quiser verificar sem simular ou selecionar.
SOMENTE_AUDITORIA <- FALSE

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), sprintf(...), "\n", sep = "")
validar <- function(teste, mensagem) if (!isTRUE(teste)) stop(mensagem, call. = FALSE)

# Salvamento atomico: evita checkpoints incompletos caso a sessao interrompa.
salvar_rds <- function(objeto, caminho) {
  dir.create(dirname(caminho), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".temporario_", tmpdir = dirname(caminho))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(objeto, tmp)
  validar(!file.exists(caminho), paste("Arquivo ja existe; nao sobrescrever:", caminho))
  validar(file.rename(tmp, caminho), paste("Erro ao finalizar arquivo:", caminho))
}
# Os CSVs sao produtos derivados e podem ser reconstituidos dos checkpoints.
escrever_csv_atomico <- function(dados, caminho) {
  tmp <- tempfile(".csv_temporario_", tmpdir = dirname(caminho))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  write.csv(dados, tmp, row.names = FALSE, na = "")
  # Uma substituicao de CSV derivado nao afeta os checkpoints.
  if (file.exists(caminho)) unlink(caminho)
  validar(file.rename(tmp, caminho), paste("Erro ao salvar CSV:", caminho))
}

validar(file.exists(SCRIPT_REFERENCIA), "Script 02_monte_carlo.R nao encontrado.")
validar(file.exists(MATRIZ_REFERENCIA), "Matriz N50 nao encontrada.")
validar(file.exists(CONFIG_HISTORICA), "Configuracao historica z-score nao encontrada.")

# Extrai APENAS as definicoes de funcoes, sem executar o Monte Carlo original.
nomes_funcoes <- c(
  "transformar_W", "simular_gl", "criar_variaveis", "entropia_rotulos",
  "codigos_condicionantes", "preparar_estatistica", "informacao_rapida",
  "selecionar_alvo", "entropia_legada", "informacao_legada", "testar_estatistica"
)
codigo <- parse(SCRIPT_REFERENCIA)
amb <- new.env(parent = asNamespace("stats"))
for (expressao in as.list(codigo)) {
  if (is.call(expressao) && identical(expressao[[1L]], as.name("<-")) &&
      is.symbol(expressao[[2L]]) &&
      as.character(expressao[[2L]]) %in% nomes_funcoes) {
    eval(expressao, envir = amb)
  }
}
validar(all(vapply(nomes_funcoes, exists, logical(1), envir = amb, inherits = FALSE)),
        "Faltam funcoes no script de referencia. Nao prosseguir.")

W_original <- as.matrix(readRDS(MATRIZ_REFERENCIA))
storage.mode(W_original) <- "double"
cfg_hist <- readRDS(CONFIG_HISTORICA)
validar(identical(dim(W_original), c(50L, 50L)), "Matriz nao tem dimensao 50x50.")
validar(identical(rownames(W_original), colnames(W_original)), "Ordem dos neuronios inconsistente.")
validar(!anyNA(W_original) && all(is.finite(W_original)) && all(W_original >= 0), "Pesos anatomicos invalidos.")
validar(all(diag(W_original) == 0), "Autoconexoes anatomicas inesperadas.")
validar(sum(W_original > 0) == 787L, "Numero de conexoes diferente de 787.")
validar(identical(W_original, cfg_hist$W_original), "Matriz difere da historica.")
W_sim <- amb$transformar_W(W_original, "zscore_positivos")
validar(identical(W_sim, cfg_hist$W_sim), "Z-score difere da matriz historica.")
validar(identical(cfg_hist$B, B) && identical(cfg_hist$lambda, LAMBDA) &&
          identical(cfg_hist$alfa, ALFA) && identical(cfg_hist$burn_in, BURN_IN) &&
          identical(cfg_hist$semente_cadeia, SEMENTE_CADEIA) &&
          identical(cfg_hist$semente_surrogates, SEMENTE_SURROGATES),
        "Parametros diferentes dos da configuracao historica.")
validar(all(TS > 1L), "T precisa ser maior que 1.")
validar(all(TS %in% c(40000L, 80000L)), "T inesperado.")

# Codigo de referencia congelado na configuracao para evitar mistura de versoes.
md5_referencia <- unname(tools::md5sum(SCRIPT_REFERENCIA))
log_msg("Matriz N50 validada (787 arestas). Z-score identico ao historico.")
log_msg("Codigo 02_monte_carlo.R: md5=%s", md5_referencia)
log_msg("Resultados novos em: %s", SAIDA)

if (SOMENTE_AUDITORIA) {
  log_msg("Auditoria concluida. Nenhuma simulacao executada.")
} else {
  dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)
  resultados <- list()
  for (T_AMOSTRA in TS) {
    pasta <- file.path(SAIDA, paste0("T", T_AMOSTRA), "zscore_positivos", "alfa_0_05")
    dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
    cfg <- list(
      versao = VERSAO, codigo_referencia_md5 = md5_referencia,
      rede = "N50", metodo = "zscore_positivos", W_original = W_original,
      W_sim = W_sim, T = T_AMOSTRA, burn_in = BURN_IN, B = B,
      lambda = LAMBDA, alfa = ALFA, semente_cadeia = SEMENTE_CADEIA,
      semente_surrogates = SEMENTE_SURROGATES
    )
    arquivo_cfg <- file.path(pasta, "configuracao.rds")
    if (file.exists(arquivo_cfg)) {
      validar(identical(readRDS(arquivo_cfg), cfg),
              paste("Configuracao incompativel; parada por seguranca:", pasta))
    } else {
      salvar_rds(cfg, arquivo_cfg)
    }
    
    arq_cadeia <- file.path(pasta, "cadeia.rds")
    if (file.exists(arq_cadeia)) {
      cadeia <- readRDS(arq_cadeia)
      log_msg("T=%d: cadeia recuperada", T_AMOSTRA)
    } else {
      log_msg("T=%d: simulando cadeia", T_AMOSTRA)
      set.seed(SEMENTE_CADEIA)
      cadeia <- amb$simular_gl(W_sim, T_AMOSTRA, BURN_IN)
      salvar_rds(cadeia, arq_cadeia)
      log_msg("T=%d: cadeia salva", T_AMOSTRA)
    }
    X <- cadeia$X
    validar(identical(dim(X), c(T_AMOSTRA, 50L)), "Dimensao da cadeia incorreta.")
    validar(identical(colnames(X), colnames(W_original)), "Ordem dos neuronios da cadeia incorreta.")
    validar(all(X %in% 0:1), "Valores de disparo invalidos.")
    validar(identical(dim(cadeia$P), c(T_AMOSTRA, 50L)), "Matriz P invalida.")
    validar(identical(dim(cadeia$V), c(T_AMOSTRA, 50L)), "Matriz V invalida.")
    
    W_est <- matrix(0L, 50L, 50L, dimnames = dimnames(W_original))
    tracos <- list()
    for (i in seq_len(50L)) {
      arq_alvo <- file.path(pasta, sprintf("alvo_%02d.rds", i))
      if (file.exists(arq_alvo)) {
        obj <- readRDS(arq_alvo)
        validar(identical(obj$config, cfg) && identical(obj$i, i),
                paste("Checkpoint incompativel:", arq_alvo))
        selecao <- obj$selecao
        log_msg("T=%d alvo=%02d/50: recuperado", T_AMOSTRA, i)
      } else {
        log_msg("T=%d alvo=%02d/50: iniciando", T_AMOSTRA, i)
        inicio <- Sys.time()
        set.seed(SEMENTE_SURROGATES + i)
        selecao <- amb$selecionar_alvo(X, i, B, ALFA, LAMBDA)
        validar(is.list(selecao) && all(c("pais", "traco") %in% names(selecao)),
                "Resultado de selecao inesperado.")
        validar(!anyDuplicated(selecao$pais) &&
                  all(selecao$pais %in% setdiff(seq_len(50L), i)),
                "Selecao com candidatos invalidos.")
        salvar_rds(list(config = cfg, i = i, selecao = selecao), arq_alvo)
        log_msg("T=%d alvo=%02d/50: salvo em %.1f min", T_AMOSTRA, i,
                as.numeric(difftime(Sys.time(), inicio, units = "mins")))
      }
      validar(all(selecao$pais %in% setdiff(seq_len(50L), i)), "Checkpoint com pais invalidos.")
      if (length(selecao$pais)) W_est[i, selecao$pais] <- 1L
      if (nrow(selecao$traco)) {
        tr <- selecao$traco
        tr$T <- T_AMOSTRA
        tr$alvo <- rownames(W_original)[i]
        tracos[[length(tracos) + 1L]] <- tr
      }
    }
    
    idx <- row(W_est) != col(W_est)
    v <- (W_original[idx] > 0)
    e <- (W_est[idx] == 1L)
    TP <- sum(v & e); FP <- sum(!v & e)
    FN <- sum(v & !e); TN <- sum(!v & !e)
    validar(TP + FP + FN + TN == 2450L && TP + FN == 787L && TN + FP == 1663L,
            "Matriz de confusao invalida.")
    metricas <- data.frame(
      rede = "N50", metodo = "zscore_positivos", T = T_AMOSTRA,
      burn_in = BURN_IN, B = B, lambda = LAMBDA, alfa = ALFA,
      TP = TP, FP = FP, FN = FN, TN = TN,
      sensibilidade = TP / (TP + FN),
      especificidade = TN / (TN + FP),
      precisao = if (TP + FP > 0) TP / (TP + FP) else NA_real_,
      F1 = if (2 * TP + FP + FN > 0) 2 * TP / (2 * TP + FP + FN) else NA_real_,
      n_selecionadas = TP + FP,
      prop_saturacao_099 = mean(cadeia$P > 0.99),
      prop_saturacao_0999 = mean(cadeia$P > 0.999),
      taxa_disparo_media = mean(X)
    )
    arq_estimativa <- file.path(pasta, "W_estimada.rds")
    if (file.exists(arq_estimativa)) {
      validar(identical(readRDS(arq_estimativa), W_est),
              paste("Matriz estimada anterior e checkpoints divergem:", pasta))
    } else {
      salvar_rds(W_est, arq_estimativa)
    }
    escrever_csv_atomico(metricas, file.path(pasta, "metricas.csv"))
    if (length(tracos)) escrever_csv_atomico(do.call(rbind, tracos), file.path(pasta, "passos_mrmr.csv"))
    resultados[[length(resultados) + 1L]] <- metricas
    escrever_csv_atomico(do.call(rbind, resultados), file.path(SAIDA, "resumo_parcial.csv"))
    log_msg("CONCLUIDO T=%d | TP=%d FP=%d FN=%d TN=%d | sens=%.4f | F1=%.4f",
            T_AMOSTRA, TP, FP, FN, TN, metricas$sensibilidade, metricas$F1)
  }
  escrever_csv_atomico(do.call(rbind, resultados), file.path(SAIDA, "resumo_geral.csv"))
  log_msg("TODOS OS EXPERIMENTOS Z-SCORE CONCLUIDOS.")
}
