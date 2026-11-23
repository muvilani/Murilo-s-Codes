# TCC — Secao 4.1: bateria completa, cenarios 3 -> 2 -> 1
# R base; checkpoint por replica; resultados informacionais.
# As metricas novas exigem RECALCULAR as replicacoes originais.
# NAO usa checkpoints da bateria antiga; pasta separada e assinatura de configuracao.

TESTE_RAPIDO <- FALSE
FORCAR_RECALCULO <- FALSE
# Depois de examinar os graficos dos cenarios 2 e 1, ajuste seus valores aqui.
BURNIN <- c(`3` = 2000L, `2` = 2000L, `1` = 2000L)
# 2000 para cenarios 2 e 1 e PROVISORIO, nao validado pelos pilotos.
if (TESTE_RAPIDO) {
  cenarios <- c(3L, 2L, 1L)
  tamanhos <- 300L
  repeticoes <- 2L
  B <- 19L
  lambdas <- c(0, 1)
  alfas <- 0.05
} else {
  cenarios <- c(3L, 2L, 1L)
  tamanhos <- c(1200L, 6000L, 12000L, 60000L, 120000L)
  repeticoes <- 100L
  B <- 199L
  lambdas <- c(0, 0.25, 0.50, 0.75, 1)
  alfas <- c(0.01, 0.05, 0.10)
}

PASTA <- if (TESTE_RAPIDO) "resultados_4_1_v2_TESTE" else "resultados_4_1_v2_COMPLETA"
dir.create(PASTA, recursive = TRUE, showWarnings = FALSE)
# Cada cenario processa TODOS os tamanhos antes de passar ao proximo.

# 1. Matrizes, simulador e funcoes informacionais (mesmas do original)
matriz_cenario <- function(cenario) {
  if (cenario == 1) {
    W <- matrix(c(
      0,  0,  1,  1,  1,
      0,  0,  1,  1,  1,
      1,  1,  0,  1, -4,
      1,  1,  1,  0, -4,
      1,  1, -4, -4,  0
    ), nrow = 5, byrow = TRUE)
  } else if (cenario == 2) {
    W <- matrix(c(
      0,  0,   3,   3,   3,
      0,  0,   3,   3,   3,
      3,  3,   0,   3, -12,
      3,  3,   3,   0, -12,
      3,  3, -12, -12,   0
    ), nrow = 5, byrow = TRUE)
  } else if (cenario == 3) {
    W <- matrix(c(
      0,  0,   3,  3,   3,
      0,  0,   1,  1,   1,
      3,  1,   0,  1, -12,
      3,  1,   1,  0,  -4,
      3,  1, -12, -4,   0
    ), nrow = 5, byrow = TRUE)
  } else {
    stop("Cenario invalido: escolha 1, 2 ou 3.")
  }
  return(W)
}

# =============================================================================
# 2. SIMULADOR GL COM DECAIMENTO 1/2 (CONFORME FERREIRA ET AL.)
#
# Em cada instante t:
#   potencial_i = soma_j W[i,j] * memoria[i,j]
#   P(X[t,i] = 1 | passado) = plogis(potencial_i)
#
# memoria[i,j] registra os spikes passados de j com decaimento 1/2,
# contados somente APOS o ultimo spike do neuronio-alvo i.
# Se i disparar agora, sua linha da memoria volta a zero.
# Nota: o artigo usa plogis(potencial), SEM intercepto logit(0.02).
# =============================================================================

simular_gl <- function(W, tamanho, burn_in) {
  n <- nrow(W)
  total <- tamanho + burn_in
  X <- matrix(0L, nrow = total, ncol = n)
  memoria <- matrix(0, nrow = n, ncol = n)
  
  for (t in seq_len(total)) {
    potencial <- rowSums(W * memoria)
    probabilidade <- plogis(potencial)
    X[t, ] <- rbinom(n, size = 1, prob = probabilidade)
    
    # Para o proximo instante: eventos antigos pesam metade e o spike
    # atual, se houver, passa a contar com peso 1.
    memoria <- memoria / 2 + matrix(X[t, ], nrow = n, ncol = n,
                                    byrow = TRUE)
    
    # Depois de um spike do alvo, descartamos o passado anterior a ele.
    memoria[X[t, ] == 1L, ] <- 0
  }
  
  # Descarta o periodo de aquecimento.
  return(X[(burn_in + 1):total, , drop = FALSE])
}

# =============================================================================
# 3. CONSTRUIR U: NUMERO DE SPIKES ENTRE O ULTIMO SPIKE DE i E t
#
# Para cada t aproveitado:
# U[t,j] = soma X[s,j] para ultimo_spike_i < s < t.
# Usamos soma acumulada para nao percorrer o intervalo a cada instante.
# =============================================================================

criar_variaveis <- function(X, i) {
  Tobs <- nrow(X)
  n <- ncol(X)
  
  # ultimo[t] = posicao do ultimo spike de i ANTES de t.
  ultimo <- c(0L, cummax(ifelse(X[1:(Tobs - 1), i] == 1L,
                                seq_len(Tobs - 1), 0L)))
  tempos <- which(ultimo > 0L)
  if (length(tempos) == 0L) stop("O alvo nao disparou na amostra.")
  
  U <- matrix(0L, nrow = length(tempos), ncol = n)
  for (j in seq_len(n)) {
    acumulada <- cumsum(X[, j])
    U[, j] <- acumulada[tempos - 1L] - acumulada[ultimo[tempos]]
  }
  
  return(list(y = X[tempos, i], U = U,
              tempos = tempos, ultimo = ultimo[tempos]))
}

# =============================================================================
# 4. INFORMACAO MUTUA EMPIRICA (ESTIMADOR PLUG-IN)
#
# Entropia: H(Z) = - soma_z p(z) log(p(z))
# Informacao condicional: I(U;Y|S) = H(U,S)+H(Y,S)-H(S)-H(U,Y,S)
# Com S vazio, volta a I(U;Y) = H(U)+H(Y)-H(U,Y).
# O logaritmo natural mede informacao em nats.
# =============================================================================

entropia <- function(dados) {
  if (is.null(dim(dados))) dados <- matrix(dados, ncol = 1)
  # Transforma cada linha (configuracao) em um texto, ex. "0|2|1".
  padroes <- do.call(paste, c(as.data.frame(dados), sep = "|"))
  frequencias <- as.numeric(table(padroes))
  p <- frequencias / sum(frequencias)
  return(-sum(p * log(p)))
}

informacao <- function(u, y, S) {
  if (ncol(S) == 0L) {
    valor <- entropia(u) + entropia(y) - entropia(cbind(u, y))
  } else {
    valor <- entropia(cbind(u, S)) + entropia(cbind(y, S)) -
      entropia(S) - entropia(cbind(u, y, S))
  }
  # Trata apenas possiveis erros numericos na casa de 10^-15.
  return(max(0, valor))
}

# =============================================================================
# 5. SELECAO DE UM ALVO (TODOS OS LAMBDAS E ALFAS)
#
# A cada etapa:
#   A) calcula a informacao condicional dos candidatos restantes;
#   B) cria B conjuntos surrogate, deslocando circularmente cada candidato;
#   C) usa o maximo das informacoes de cada surrogate para obter o limiar;
#   D) permite apenas candidatos com informacao condicional > limiar;
#   E) entre eles, escolhe quem tiver o maior escore Q_lambda.
#
# A secao 4.1 usa I(U;Y|S) para PASSAR NO TESTE e Q_lambda para ESCOLHER.
# Uma pequena memoria (cache) evita recalcular uma etapa ja visitada para
# outro valor de lambda ou alpha, deixando a bateria bem menos custosa.
# =============================================================================


# 2. Selecao + auditoria das medidas informacionais
selecionar_alvo <- function(X, i, B, lambdas, alfas) {
  # Caso excepcional: o alvo nao teve disparo ANTES do ultimo instante.
  # Assim, nao existe nenhum t da amostra com ultimo_spike_i < t.
  # O procedimento original parava em criar_variaveis(); agora devolve
  # conjunto vazio para todos os parametros, sem inventar informacao.
  # Interpretacao: sem evidencias para selecionar pais, NAO ausencia de arestas.
  # Inclui alvo sem disparos ou com unico disparo no ultimo instante.
  if (!any(X[seq_len(nrow(X) - 1L), i] == 1L)) {
    message(sprintf(
      "AVISO [T=%d, alvo=%d]: sem disparo anterior utilizavel; nenhum pai sera selecionado para esse alvo.",
      nrow(X), i
    ))
    conjuntos_vazios <- list()
    for (lambda in lambdas) {
      for (alfa in alfas) {
        conjuntos_vazios[[paste(lambda, alfa, sep = "_")]] <- integer(0)
      }
    }
    return(list(conjuntos = conjuntos_vazios, candidatos = data.frame(),
                selecionados = data.frame()))
  }
  
  dados <- criar_variaveis(X, i)
  U <- dados$U
  y <- dados$y
  n <- ncol(X)
  candidatos_todos <- setdiff(seq_len(n), i)
  
  vazio <- matrix(nrow = length(y), ncol = 0)
  relevancia <- sapply(candidatos_todos, function(j) {
    informacao(U[, j], y, vazio)
  })
  names(relevancia) <- as.character(candidatos_todos)
  
  cache <- new.env(parent = emptyenv())
  
  # Calcula um passo para um determinado conjunto de selecionados S.
  avaliar_passo <- function(S) {
    chave <- paste0("S_", paste(sort(S), collapse = "_"))
    if (exists(chave, envir = cache, inherits = FALSE)) {
      return(get(chave, envir = cache, inherits = FALSE))
    }
    
    candidatos <- setdiff(candidatos_todos, S)
    condicionantes <- U[, S, drop = FALSE]
    
    observado <- sapply(candidatos, function(j) {
      informacao(U[, j], y, condicionantes)
    })
    names(observado) <- as.character(candidatos)
    
    maximos_nulos <- numeric(B)
    for (b in seq_len(B)) {
      valores_nulos <- numeric(length(candidatos))
      
      for (k in seq_along(candidatos)) {
        j <- candidatos[k]
        serie <- X[, j]
        deslocamento <- sample.int(nrow(X) - 1L, size = 1L)
        
        # Circular shift: desloca j sem mexer no alvo i nem nos
        # condicionantes S. O comprimento da serie e preservado.
        serie_surrogate <- c(tail(serie, deslocamento),
                             head(serie, nrow(X) - deslocamento))
        
        # Reconstrucao de U usando o MESMO ultimo spike do alvo original.
        acumulada <- cumsum(serie_surrogate)
        u_surrogate <- acumulada[dados$tempos - 1L] -
          acumulada[dados$ultimo]
        
        valores_nulos[k] <- informacao(u_surrogate, y, condicionantes)
      }
      # Maximo entre candidatos = corrige a busca por varios candidatos.
      maximos_nulos[b] <- max(valores_nulos)
    }
    
    limites <- sapply(alfas, function(a) {
      as.numeric(quantile(maximos_nulos, probs = 1 - a, type = 1))
    })
    names(limites) <- as.character(alfas)
    
    resultado <- list(candidatos = candidatos, observado = observado,
                      limites = limites)
    assign(chave, resultado, envir = cache)
    return(resultado)
  }
  
  resultados <- list()
  linhas_candidatos <- list()
  linhas_selecionados <- list()
  for (lambda in lambdas) {
    for (alfa in alfas) {
      selecionados <- integer(0)
      etapa <- 0L
      repeat {
        if (length(selecionados) == length(candidatos_todos)) break
        etapa <- etapa + 1L
        passo <- avaliar_passo(selecionados)
        limiar <- passo$limites[as.character(alfa)]
        
        # So entram na disputa os candidatos com informacao condicional
        # MAIOR do que o limiar obtido pelos surrogates.
        aprovados <- passo$candidatos[passo$observado > limiar]
        # Armazena TODOS os candidatos avaliados na etapa, inclusive os reprovados.
        nomes <- as.character(passo$candidatos)
        rel <- as.numeric(relevancia[nomes])
        inc <- as.numeric(passo$observado[nomes])
        red <- rel - inc  # Pode ser negativa: informacao condicional > marginal.
        notas_todas <- (1 - lambda) * rel + lambda * inc
        escolhido <- NA_integer_
        if (length(aprovados) > 0L) {
          notas <- (1 - lambda) * relevancia[as.character(aprovados)] +
            lambda * passo$observado[as.character(aprovados)]
          escolhido <- aprovados[which.max(notas)]
        }
        linhas_candidatos[[length(linhas_candidatos) + 1L]] <- data.frame(
          alvo = i, lambda = lambda, alfa = alfa, etapa = etapa,
          candidato = passo$candidatos, relevancia = rel, incremental = inc,
          redundancia = red, escore = notas_todas, limiar = as.numeric(limiar),
          elegivel = passo$candidatos %in% aprovados,
          selecionado = !is.na(escolhido) & passo$candidatos == escolhido
        )
        if (length(aprovados) == 0L) break
        
        # Mesma escolha do algoritmo original; apenas registramos as medidas.
        linhas_selecionados[[length(linhas_selecionados) + 1L]] <- data.frame(
          alvo = i, lambda = lambda, alfa = alfa, etapa = etapa,
          candidato = escolhido,
          relevancia = as.numeric(relevancia[as.character(escolhido)]),
          incremental = as.numeric(passo$observado[as.character(escolhido)]),
          redundancia = as.numeric(relevancia[as.character(escolhido)] -
                                     passo$observado[as.character(escolhido)])
        )
        selecionados <- c(selecionados, escolhido)
      }
      
      chave_resultado <- paste(lambda, alfa, sep = "_")
      resultados[[chave_resultado]] <- selecionados
    }
  }
  combinar <- function(lista) if (length(lista)) do.call(rbind, lista) else data.frame()
  return(list(conjuntos = resultados, candidatos = combinar(linhas_candidatos),
              selecionados = combinar(linhas_selecionados)))
}



# 3. Checagens basicas ---------------------------------------------------------
X_pequeno <- cbind(c(1, 0, 0, 1, 0, 0), c(0, 1, 1, 0, 1, 0))
stopifnot(identical(as.integer(criar_variaveis(X_pequeno, 1)$U[, 2]),
                    c(0L, 1L, 2L, 0L, 1L)))
z <- c(0, 1, 0, 1)
stopifnot(abs(informacao(z, z, matrix(nrow = 4, ncol = 0)) - log(2)) < 1e-10)

salvar_seguro <- function(obj, caminho) {
  tmp <- paste0(caminho, ".tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, caminho)) {
    if (file.exists(caminho)) file.remove(caminho)
    if (!file.rename(tmp, caminho)) stop("Falha ao salvar ", caminho)
  }
}

# Validacao impede reutilizar checkpoints de parametros diferentes.
assinatura <- function(cenario, T) list(versao = "v2_info_20261019",
                                        cenario = cenario, T = T, repeticoes = repeticoes,
                                        burn_in = as.integer(BURNIN[as.character(cenario)]), B = B,
                                        lambdas = lambdas, alfas = alfas)
valido <- function(x, ass, replica) {
  if (!is.list(x) || !identical(x$assinatura, ass)) return(FALSE)
  d <- x$resumo
  if (!is.data.frame(d) || nrow(d) != length(lambdas)*length(alfas)) return(FALSE)
  if (!all(d$replica == replica)) return(FALSE)
  esperado <- expand.grid(lambda = lambdas, alfa = alfas)
  observado <- d[, c("lambda", "alfa")]
  identical(sort(do.call(paste, c(esperado, sep = ":"))),
            sort(do.call(paste, c(observado, sep = ":")))) &&
    is.data.frame(x$candidatos) && is.data.frame(x$selecionados)
}

# Colunas de candidatos / selecionados que permanecem quando nenhuma aresta entra.
vazio_detalhes <- function(selecionado = FALSE) {
  x <- data.frame(alvo = integer(), lambda = numeric(), alfa = numeric(),
                  etapa = integer(), candidato = integer(), relevancia = numeric(),
                  incremental = numeric(), redundancia = numeric())
  if (!selecionado) {
    x$escore <- numeric(); x$limiar <- numeric()
    x$elegivel <- logical(); x$selecionado <- logical()
  }
  x
}

# Agregacao das replicas concluidas: permite inspecionar parciais e retomar.
atualizar_resumos <- function(pasta_cenario, ass, T) {
  pasta_t <- file.path(pasta_cenario, paste0("T_", T))
  paths <- file.path(pasta_t, sprintf("rep_%03d.rds", seq_len(repeticoes)))
  todos <- lapply(seq_along(paths), function(r) {
    if (!file.exists(paths[r])) return(NULL)
    obj <- tryCatch(readRDS(paths[r]), error = function(e) NULL)
    if (valido(obj, ass, r)) obj$resumo else NULL
  })
  todos <- Filter(Negate(is.null), todos)
  if (!length(todos)) return(invisible(NULL))
  tab <- do.call(rbind, todos)
  write.csv(tab, file.path(pasta_t, "replicas_parcial.csv"), row.names = FALSE)
  met <- c("TP", "FP", "FN", "TN", "sensibilidade", "especificidade",
           "proporcao_alvos_exatos", "grafo_exato", "relevancia_media",
           "redundancia_media", "incremental_media", "n_selecionados")
  agrup <- list(cenario = tab$cenario, T = tab$T, lambda = tab$lambda, alfa = tab$alfa)
  medias <- aggregate(tab[met], by = agrup,
                      FUN = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  desvios <- aggregate(tab[met], by = agrup,
                       FUN = function(x) if (sum(is.finite(x)) < 2) NA_real_ else sd(x, na.rm = TRUE))
  names(desvios)[-(1:4)] <- paste0(names(desvios)[-(1:4)], "_dp")
  nrep <- aggregate(tab$replica, by = agrup, FUN = length)
  names(nrep)[5] <- "n_replicas"
  resumo <- merge(merge(medias, desvios, by = c("cenario", "T", "lambda", "alfa")),
                  nrep, by = c("cenario", "T", "lambda", "alfa"))
  resumo <- resumo[order(resumo$alfa, resumo$lambda), ]
  write.csv(resumo, file.path(pasta_t, "resumo_parcial.csv"), row.names = FALSE)
  invisible(resumo)
}

# 4. Bateria ------------------------------------------------------------------
for (cenario in cenarios) {
  W <- matriz_cenario(cenario)
  bc <- as.integer(BURNIN[as.character(cenario)])
  if (is.na(bc) || bc < 0) stop("Burn-in invalido para cenario ", cenario)
  pasta_c <- file.path(PASTA, paste0("cenario_", cenario))
  dir.create(pasta_c, recursive = TRUE, showWarnings = FALSE)
  cat("\n=========================================\n", "CENARIO ", cenario,
      " | burn-in = ", bc, "\n=========================================\n", sep="")
  for (tamanho in tamanhos) {
    ass <- assinatura(cenario, tamanho)
    pasta_t <- file.path(pasta_c, paste0("T_", tamanho))
    dir.create(pasta_t, recursive = TRUE, showWarnings = FALSE)
    t_inicio <- Sys.time()
    for (r in seq_len(repeticoes)) {
      arquivo <- file.path(pasta_t, sprintf("rep_%03d.rds", r))
      if (file.exists(arquivo) && !FORCAR_RECALCULO) {
        salvo <- tryCatch(readRDS(arquivo), error = function(e) NULL)
        if (valido(salvo, ass, r)) {
          cat(sprintf("C%d | T=%d | replica=%d/%d — checkpoint valido\n",
                      cenario, tamanho, r, repeticoes))
          next
        }
      }
      set.seed(170000 + cenario * 100000 + tamanho + r)
      cat(sprintf("C%d | T=%d | replica=%d/%d — calculando\n",
                  cenario, tamanho, r, repeticoes))
      X <- simular_gl(W, tamanho = tamanho, burn_in = bc)
      n <- ncol(X)
      estimados <- lapply(seq_len(n), function(i) selecionar_alvo(X, i, B, lambdas, alfas))
      candidatos <- lapply(estimados, `[[`, "candidatos")
      candidatos <- candidatos[vapply(candidatos, nrow, integer(1)) > 0]
      candidatos <- if (length(candidatos)) do.call(rbind, candidatos) else vazio_detalhes()
      selecionados <- lapply(estimados, `[[`, "selecionados")
      selecionados <- selecionados[vapply(selecionados, nrow, integer(1)) > 0]
      selecionados <- if (length(selecionados)) do.call(rbind, selecionados) else vazio_detalhes(TRUE)
      candidatos$cenario <- cenario; candidatos$T <- tamanho; candidatos$replica <- r
      selecionados$cenario <- cenario; selecionados$T <- tamanho; selecionados$replica <- r
      linhas <- list()
      for (lambda in lambdas) for (alfa in alfas) {
        chave <- paste(lambda, alfa, sep="_")
        tp <- fp <- fn <- tn <- exatos <- 0L
        for (i in seq_len(n)) {
          verdadeiro <- setdiff(which(W[i, ] != 0), i)
          estimado <- estimados[[i]]$conjuntos[[chave]]
          tp <- tp + length(intersect(verdadeiro, estimado))
          fp <- fp + length(setdiff(estimado, verdadeiro))
          fn <- fn + length(setdiff(verdadeiro, estimado))
          tn <- tn + length(setdiff(setdiff(seq_len(n), i), union(verdadeiro, estimado)))
          exatos <- exatos + as.integer(setequal(verdadeiro, estimado))
        }
        ss <- selecionados[selecionados$lambda == lambda & selecionados$alfa == alfa, , drop=FALSE]
        medio <- function(x) if (!length(x)) NA_real_ else mean(x)
        linhas[[length(linhas)+1L]] <- data.frame(
          cenario = cenario, T = tamanho, replica = r, lambda = lambda, alfa = alfa,
          TP = tp, FP = fp, FN = fn, TN = tn,
          sensibilidade = tp/(tp+fn), especificidade = tn/(tn+fp),
          proporcao_alvos_exatos = exatos/n, grafo_exato = as.integer(fp == 0 && fn == 0),
          relevancia_media = medio(ss$relevancia),
          redundancia_media = medio(ss$redundancia),
          incremental_media = medio(ss$incremental), n_selecionados = nrow(ss))
      }
      obj <- list(assinatura = ass, resumo = do.call(rbind, linhas),
                  candidatos = candidatos, selecionados = selecionados)
      stopifnot(valido(obj, ass, r))
      salvar_seguro(obj, arquivo)
      # Atualiza o CSV a cada 5 replicas e na primeira, para acompanhar andamento.
      if (r == 1L || r %% 5L == 0L || r == repeticoes)
        atualizar_resumos(pasta_c, ass, tamanho)
    }
    atualizar_resumos(pasta_c, ass, tamanho)
    cat("CONCLUIDO C", cenario, " T=", tamanho, " | duracao = ",
        round(as.numeric(difftime(Sys.time(), t_inicio, units="mins")), 1), " min\n", sep="")
  }
}

# 5. Consolidacao da tabela principal ----------------------------------------
resumos <- list()
for (cenario in cenarios) for (tamanho in tamanhos) {
  arq <- file.path(PASTA, paste0("cenario_", cenario), paste0("T_", tamanho), "resumo_parcial.csv")
  if (file.exists(arq)) resumos[[length(resumos)+1L]] <- read.csv(arq)
}
if (length(resumos)) {
  tabela <- do.call(rbind, resumos)
  tabela <- tabela[order(match(tabela$cenario, cenarios), tabela$T, tabela$alfa, tabela$lambda), ]
  write.csv(tabela, file.path(PASTA, "tabela_completa_por_configuracao.csv"), row.names = FALSE)
  # A tabela de apresentacao preserva sensibilidade, especificidade e alvos exatos.
  mostrar <- tabela[, c("cenario", "T", "alfa", "lambda", "sensibilidade",
                        "especificidade", "proporcao_alvos_exatos",
                        "relevancia_media", "redundancia_media", "n_replicas")]
  mostrar[, c("sensibilidade", "especificidade", "proporcao_alvos_exatos")] <-
    100 * mostrar[, c("sensibilidade", "especificidade", "proporcao_alvos_exatos")]
  write.csv(mostrar, file.path(PASTA, "tabela_apresentacao.csv"), row.names = FALSE)
}
cat("\nFim. Resultados em: ", normalizePath(PASTA), "\n", sep="")
cat("Detalhes por candidato/etapa: campos 'candidatos' e 'selecionados' de cada rep_XXX.rds.\n")
