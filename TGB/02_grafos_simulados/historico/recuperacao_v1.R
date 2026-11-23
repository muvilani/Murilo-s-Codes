# =============================================================================
# TCC - SECAO 4.1: RECUPERACAO DE GRAFOS NEURONAIS
# VERSAO DE RECUPERACAO PARA MAQUINA VIRTUAL - CENARIOS 1 E 2
# CHECKPOINT AUTOMATICO POR REPLICACAO
# CORRECAO: TRATAMENTO DE ALVO SEM DISPARO ANTERIOR UTILIZAVEL
# Cassiavilani e Ferreira (manuscrito em desenvolvimento)
#
# OBJETIVO
# Simular os cenarios 1, 2 e 3 de Ferreira et al. (2025), aplicar a selecao
# baseada em informacao mutua e usar surrogates para decidir QUANDO PARAR.
#
# NENHUM PACOTE E NECESSARIO: todo o codigo usa apenas o R base.
# Execute no RStudio com Source ou selecione os blocos e use Ctrl+Enter.
# =============================================================================

# =============================================================================
# 0. CONFIGURACOES: COMECE COM TESTE_RAPIDO = TRUE
# =============================================================================

TESTE_RAPIDO <- FALSE
FORCAR_RECALCULO <- FALSE
# Cada replicacao recebe uma semente propria (veja Secao 7).
set.seed(17)

if (TESTE_RAPIDO) {
  cenarios <- 1
  tamanhos <- 300
  repeticoes <- 2
  B <- 19
  lambdas <- c(0, 1)
  alfas <- 0.05
} else {
  cenarios <- 1:2
  tamanhos <- c(1200, 6000, 12000, 60000)
  repeticoes <- 100
  B <- 199  # Escolha de implementacao; a secao 4.1 nao fixa B.
  lambdas <- c(0, 0.25, 0.50, 0.75, 1)
  alfas <- c(0.01, 0.05, 0.10)
}

burn_in <- 2000  # Escolha de implementacao; fixa em todas as configuracoes.
# Mesmo diretorio da execucao original da maquina virtual: reutiliza os .rds prontos.
pasta_resultados <- if (TESTE_RAPIDO) "resultados_4_1_teste" else "resultados_4_1"
dir.create(pasta_resultados, showWarnings = FALSE)

# AVISO: a bateria completa e muito custosa por causa dos B surrogates por
# etapa, por alvo e por replica. Primeiro valide o TESTE_RAPIDO.

# =============================================================================
# 1. AS TRES MATRIZES VERDADEIRAS DE FERREIRA ET AL. (2025)
# W[i,j] e o peso da conexao j -> i.
# Zero = nao existe aresta. Peso negativo = efeito inibitorio.
# =============================================================================

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
    return(conjuntos_vazios)
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
  for (lambda in lambdas) {
    for (alfa in alfas) {
      selecionados <- integer(0)
      
      repeat {
        if (length(selecionados) == length(candidatos_todos)) break
        passo <- avaliar_passo(selecionados)
        limiar <- passo$limites[as.character(alfa)]
        
        # So entram na disputa os candidatos com informacao condicional
        # MAIOR do que o limiar obtido pelos surrogates.
        aprovados <- passo$candidatos[passo$observado > limiar]
        if (length(aprovados) == 0L) break
        
        # Entre os aprovados, escolher usando o escore mRMR generalizado.
        notas <- (1 - lambda) * relevancia[as.character(aprovados)] +
          lambda * passo$observado[as.character(aprovados)]
        escolhido <- aprovados[which.max(notas)]
        selecionados <- c(selecionados, escolhido)
      }
      
      chave_resultado <- paste(lambda, alfa, sep = "_")
      resultados[[chave_resultado]] <- selecionados
    }
  }
  return(resultados)
}

# =============================================================================
# 6. CONFERENCIAS SIMPLES ANTES DE SIMULAR
# =============================================================================

# Neste exemplo: alvo 1 dispara em t=1 e t=4; candidato 2 em t=2,3,5.
X_pequeno <- cbind(c(1, 0, 0, 1, 0, 0),
                   c(0, 1, 1, 0, 1, 0))
teste_U <- criar_variaveis(X_pequeno, i = 1)$U[, 2]
stopifnot(identical(as.integer(teste_U), c(0L, 1L, 2L, 0L, 1L)))

# Se Y=U com metade zeros e metade uns, I(U;Y)=log(2) nats.
z <- c(0, 1, 0, 1)
stopifnot(abs(informacao(z, z, matrix(nrow = 4, ncol = 0)) - log(2)) < 1e-10)

# Testes de regressao: alvo sem spikes ou com spike apenas no ultimo t.
X_zero <- cbind(rep(0L, 6), c(0L, 1L, 0L, 1L, 0L, 1L))
X_final <- cbind(c(0L, 0L, 0L, 0L, 0L, 1L), rep(0L, 6))
for (teste_X in list(X_zero, X_final)) {
  teste_vazio <- suppressMessages(
    selecionar_alvo(teste_X, i = 1, B = 1,
                    lambdas = c(0, 0.25, 0.5, 0.75, 1),
                    alfas = c(0.01, 0.05, 0.10))
  )
  stopifnot(length(teste_vazio) == 15L,
            all(vapply(teste_vazio, length, integer(1L)) == 0L))
}

cat("Testes basicos passaram!\n")

# =============================================================================
# 7. EXPERIMENTO: 1 REDE VERDADEIRA X VARIAS ESTIMATIVAS
# Avaliamos TP, FP, FN, TN, sensibilidade, especificidade e recuperacao exata.
# Os tres cenarios tem 18 arestas presentes e 2 ausentes (sem autolacos).
# =============================================================================

resultados_finais <- list()

# Testa se um conjunto de linhas esta completo para a configuracao/replica.
validar_linhas <- function(df, cenario, tamanho, replicas_esperadas) {
  chaves <- c("cenario", "T", "replica", "lambda", "alfa",
              "TP", "FP", "FN", "TN", "sensibilidade", "especificidade",
              "proporcao_alvos_exatos", "grafo_exato")
  if (!is.data.frame(df) || !all(chaves %in% names(df))) return(FALSE)
  if (nrow(df) != length(replicas_esperadas) * length(lambdas) * length(alfas))
    return(FALSE)
  if (!all(df$cenario == cenario) || !all(df$T == tamanho)) return(FALSE)
  esperado <- expand.grid(replica = replicas_esperadas,
                          lambda = lambdas, alfa = alfas)
  observado <- df[, c("replica", "lambda", "alfa"), drop = FALSE]
  chaves1 <- do.call(paste, c(esperado, sep = ":"))
  chaves2 <- do.call(paste, c(observado, sep = ":"))
  identical(sort(chaves1), sort(chaves2))
}

pasta_checkpoints <- file.path(pasta_resultados, "checkpoints")
dir.create(pasta_checkpoints, recursive = TRUE, showWarnings = FALSE)

# Escreve o RDS primeiro em arquivo temporario para reduzir risco de corrupcao.
salvar_seguro <- function(obj, caminho) {
  temporario <- paste0(caminho, ".tmp")
  saveRDS(obj, temporario)
  if (file.exists(caminho)) file.remove(caminho)
  if (!file.rename(temporario, caminho))
    stop("Falha ao salvar: ", caminho)
}

for (cenario in cenarios) {
  W <- matriz_cenario(cenario)
  
  for (tamanho in tamanhos) {
    arquivo <- file.path(pasta_resultados,
                         paste0("cenario_", cenario, "_T_", tamanho, ".rds"))
    
    if (file.exists(arquivo) && !FORCAR_RECALCULO) {
      tabela_salva <- tryCatch(readRDS(arquivo), error = function(e) NULL)
      if (validar_linhas(tabela_salva, cenario, tamanho, seq_len(repeticoes))) {
        cat("Configuracao completa: reutilizando ", arquivo, "\n", sep = "")
        resultados_finais[[length(resultados_finais) + 1L]] <- tabela_salva
        next
      }
      warning("Arquivo completo invalido: sera recalculado: ", arquivo)
    }
    
    linhas <- vector("list", repeticoes)
    for (r in seq_len(repeticoes)) {
      checkpoint <- file.path(pasta_checkpoints,
                              sprintf("cenario_%d_T_%d_rep_%03d.rds", cenario, tamanho, r))
      
      if (file.exists(checkpoint) && !FORCAR_RECALCULO) {
        salvo <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
        if (validar_linhas(salvo, cenario, tamanho, r)) {
          linhas[[r]] <- salvo
          cat("Cenario", cenario, "| T =", tamanho,
              "| Replica", r, "de", repeticoes, "-- salva, pulando\n")
          next
        }
      }
      
      # Uma nova semente para cada replica permite reiniciar com identicos dados.
      # Nao altera o metodo de simulacao, selecao ou teste.
      set.seed(170000 + cenario * 100000 + tamanho + r)
      cat("Cenario", cenario, "| T =", tamanho,
          "| Replica", r, "de", repeticoes, "-- calculando\n")
      
      linhas_replica <- list()
      X <- simular_gl(W, tamanho = tamanho, burn_in = burn_in)
      n <- ncol(X)
      
      # A selecao e feita separadamente para cada neuronio-alvo.
      estimados <- lapply(seq_len(n), function(i) {
        selecionar_alvo(X, i, B = B, lambdas = lambdas, alfas = alfas)
      })
      
      for (lambda in lambdas) {
        for (alfa in alfas) {
          tp <- fp <- fn <- tn <- 0L
          alvos_exatos <- 0L
          chave <- paste(lambda, alfa, sep = "_")
          
          for (i in seq_len(n)) {
            verdadeiro <- setdiff(which(W[i, ] != 0), i)
            estimado <- estimados[[i]][[chave]]
            
            tp <- tp + length(intersect(verdadeiro, estimado))
            fp <- fp + length(setdiff(estimado, verdadeiro))
            fn <- fn + length(setdiff(verdadeiro, estimado))
            tn <- tn + length(setdiff(setdiff(seq_len(n), i),
                                      union(verdadeiro, estimado)))
            alvos_exatos <- alvos_exatos + as.integer(setequal(verdadeiro, estimado))
          }
          
          linhas_replica[[length(linhas_replica) + 1]] <- data.frame(
            cenario = cenario, T = tamanho, replica = r,
            lambda = lambda, alfa = alfa,
            TP = tp, FP = fp, FN = fn, TN = tn,
            sensibilidade = tp / (tp + fn),
            especificidade = tn / (tn + fp),
            proporcao_alvos_exatos = alvos_exatos / n,
            grafo_exato = as.integer(fp == 0 && fn == 0)
          )
        }
      }
      
      # Somente apos a replica inteira terminar, salva seu checkpoint.
      tabela_replica <- do.call(rbind, linhas_replica)
      if (!validar_linhas(tabela_replica, cenario, tamanho, r))
        stop("Replica incompleta. Cenario: ", cenario, " T: ", tamanho,
             " Replica: ", r)
      salvar_seguro(tabela_replica, checkpoint)
      linhas[[r]] <- tabela_replica
    }
    
    tabela <- do.call(rbind, linhas)
    if (!validar_linhas(tabela, cenario, tamanho, seq_len(repeticoes)))
      stop("Configuracao incompleta: cenario ", cenario, ", T = ", tamanho)
    salvar_seguro(tabela, arquivo)
    cat("CONFIGURACAO FINALIZADA E SALVA: ", arquivo, "\n", sep = "")
    resultados_finais[[length(resultados_finais) + 1L]] <- tabela
  }
}

resultados <- do.call(rbind, resultados_finais)
write.csv(resultados, file.path(pasta_resultados, "replicas.csv"),
          row.names = FALSE)

# Medias das metricas ao longo das repeticoes de Monte Carlo.
resumo <- aggregate(
  cbind(TP, FP, FN, TN, sensibilidade, especificidade,
        proporcao_alvos_exatos, grafo_exato) ~ cenario + T + lambda + alfa,
  data = resultados, FUN = mean
)

write.csv(resumo, file.path(pasta_resultados, "resumo.csv"), row.names = FALSE)
cat("\nRESUMO DAS SIMULACOES:\n")
print(resumo)
cat("\nArquivos salvos em:", pasta_resultados, "\n")
