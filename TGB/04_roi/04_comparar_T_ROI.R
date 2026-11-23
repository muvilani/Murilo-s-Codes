# =============================================================================
# TCC | ROI LOP_L_col_09_28 | COMPARACAO T=40000 E T=80000
# CORRIGIDO: funcoes copiadas literalmente do codigo fornecido pelo autor.
# GL com decaimento; selecao com CONTAGENS sem decaimento.
# CMI condicionada no conjunto COMPLETO; surrogates independentes por candidato.
# Nao utiliza resultados da implementacao pareada anterior.
# Salve este arquivo e execute com source() no RStudio.
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




# =============================================================================
# ADAPTACAO PARA A REGIAO ANATOMICA — EXECUCAO
# =============================================================================
executar_roi_corrigida <- function() {
  PROJETO <- "/home/murilo/tcc_malecns"
  ROI <- "LOP_L_col_09_28"
  TAMANHOS <- c(40000L, 80000L)
  ALFAS <- c(0.01, 0.05, 0.10)
  LAMBDA <- 1
  B <- 199L
  BURN_IN <- 500L
  REPLICAS <- 1L
  SEED_GL <- 20261118L
  SEED_SUR <- 170000L
  
  SAIDA <- file.path(PROJETO, "comparacao_T_ROI_metodo_original_v1")
  dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)
  
  log_msg <- function(...) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
    flush.console()
  }
  salvar <- function(obj, caminho) {
    tmp <- tempfile(".temporario_", tmpdir = dirname(caminho))
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(obj, tmp)
    if (!file.rename(tmp, caminho)) stop("Falha ao salvar: ", caminho)
  }
  salvar_csv <- function(obj, caminho) {
    tmp <- tempfile(".temporario_", tmpdir = dirname(caminho))
    on.exit(unlink(tmp), add = TRUE)
    write.csv(obj, tmp, row.names = FALSE, na = "")
    if (!file.rename(tmp, caminho)) stop("Falha ao salvar CSV: ", caminho)
  }
  
  # Checagens do codigo de referencia, sem executar sua bateria antiga.
  X_pequeno <- cbind(c(1,0,0,1,0,0), c(0,1,1,0,1,0))
  stopifnot(identical(
    as.integer(criar_variaveis(X_pequeno, 1)$U[,2]),
    c(0L,1L,2L,0L,1L)
  ))
  z <- c(0,1,0,1)
  stopifnot(abs(informacao(z,z,matrix(nrow=4,ncol=0))-log(2)) < 1e-10)
  
  W <- as.matrix(readRDS(file.path(
    PROJETO, "levantamento_regioes", "matrizes", ROI, "W_regiao.rds"
  )))
  stopifnot(
    identical(dim(W), c(64L,64L)),
    all(is.finite(W)), all(W >= 0), all(diag(W) == 0),
    sum(W > 0) == 268L, sum(W) == 437,
    !is.null(rownames(W)), identical(rownames(W),colnames(W)),
    !anyDuplicated(rownames(W)),
    identical(LAMBDA, 1), all(ALFAS > 0 & ALFAS < 1),
    all(TAMANHOS > 1L), BURN_IN >= 0L
  )
  N <- nrow(W)
  IDS <- rownames(W)
  W_sim <- log1p(W) / log1p(max(W))
  diag(W_sim) <- 0
  
  nomes <- c("simular_gl", "criar_variaveis", "entropia", "informacao", "selecionar_alvo")
  codigo_funcoes <- setNames(lapply(nomes, function(nm) {
    f <- get(nm, envir = environment(simular_gl))
    list(argumentos = formals(f), corpo = deparse(body(f)))
  }), nomes)
  
  cfg <- list(
    versao = "ROI_contagens_CMI_conjunta_v1",
    roi = ROI, metodo = "log_maximo", W_original = W, W_sim = W_sim,
    tamanhos = TAMANHOS, alfas = ALFAS, lambda = LAMBDA,
    B = B, burn_in = BURN_IN, replicas = REPLICAS,
    seed_gl = SEED_GL, seed_sur = SEED_SUR,
    codigo_funcoes = codigo_funcoes,
    desenho = "uma trajetoria longa por replica; prefixos apos burn-in"
  )
  arquivo_cfg <- file.path(SAIDA, "configuracao.rds")
  if (file.exists(arquivo_cfg)) {
    if (!identical(readRDS(arquivo_cfg), cfg)) {
      stop("Configuracao ou funcoes diferentes. Use OUTRA pasta de saida.")
    }
  } else salvar(cfg, arquivo_cfg)
  
  writeLines(capture.output(sessionInfo()), file.path(SAIDA,"sessionInfo.txt"))
  
  avaliar <- function(selecao, i, tt, rep, segundos) {
    verdadeiros <- setdiff(which(W[i,] > 0), i)
    do.call(rbind, lapply(ALFAS, function(a) {
      chave <- paste(LAMBDA, a, sep = "_")
      estimados <- selecao$conjuntos[[chave]]
      if (is.null(estimados)) stop("Chave de selecao ausente: ", chave)
      stopifnot(!anyDuplicated(estimados),
                all(estimados %in% setdiff(seq_len(N),i)))
      tp <- length(intersect(verdadeiros,estimados))
      fp <- length(setdiff(estimados,verdadeiros))
      fn <- length(setdiff(verdadeiros,estimados))
      tn <- N-1L-tp-fp-fn
      data.frame(roi=ROI,metodo="log_maximo",T=tt,replica=rep,
                 lambda=LAMBDA,alfa=a,alvo=IDS[i],
                 n_pais_verdadeiros=length(verdadeiros),
                 n_selecionados=length(estimados),TP=tp,FP=fp,FN=fn,TN=tn,
                 exato=as.integer(setequal(verdadeiros,estimados)),
                 segundos_alvo=segundos)
    }))
  }
  resumir <- function(d) {
    grupos <- split(d,interaction(d$T,d$replica,d$alfa,drop=TRUE))
    tab <- do.call(rbind,lapply(grupos,function(q) {
      tp<-sum(q$TP); fp<-sum(q$FP); fn<-sum(q$FN); tn<-sum(q$TN)
      data.frame(T=q$T[1],replica=q$replica[1],alfa=q$alfa[1],lambda=LAMBDA,
                 n_alvos=nrow(q),completo=length(unique(q$alvo))==N,
                 TP=tp,FP=fp,FN=fn,TN=tn,
                 precisao=if(tp+fp>0) tp/(tp+fp) else NA_real_,
                 sensibilidade=if(tp+fn>0) tp/(tp+fn) else NA_real_,
                 especificidade=if(tn+fp>0) tn/(tn+fp) else NA_real_,
                 F1=if(2*tp+fp+fn>0) 2*tp/(2*tp+fp+fn) else NA_real_,
                 n_exatos=sum(q$exato),fracao_exatos=mean(q$exato))
    }))
    rownames(tab)<-NULL
    tab[order(tab$replica,tab$T,tab$alfa),]
  }
  
  # Reconstroi os CSVs a partir dos checkpoints presentes, inclusive parciais.
  atualizar <- function() {
    arquivos <- list.files(SAIDA,pattern="^alvo_[0-9]+\\.rds$",
                           recursive=TRUE,full.names=TRUE)
    if (!length(arquivos)) return(invisible(NULL))
    metricas <- lapply(arquivos,function(f) {
      obj<-readRDS(f)
      if (!identical(obj$config,cfg)) stop("Checkpoint incompativel: ",f)
      obj$metricas
    })
    d<-do.call(rbind,metricas)
    stopifnot(!anyDuplicated(d[c("T","replica","alfa","alvo")]),
              all(d$TP+d$FN==d$n_pais_verdadeiros),
              all(d$TP+d$FP==d$n_selecionados),
              all(d$TP+d$FP+d$FN+d$TN==N-1L))
    s<-resumir(d)
    salvar_csv(d,file.path(SAIDA,"metricas_por_alvo.csv"))
    salvar_csv(s,file.path(SAIDA,"resumo_comparacao_T.csv"))
    invisible(s)
  }
  
  atualizar()
  for (rep in REPLICAS) {
    pasta_rep <- file.path(SAIDA,sprintf("replica_%02d",rep))
    dir.create(pasta_rep,showWarnings=FALSE)
    arquivo_X <- file.path(pasta_rep,"cadeia_longa.rds")
    if (file.exists(arquivo_X)) {
      cadeia<-readRDS(arquivo_X)
      if (!identical(cadeia$config,cfg)) stop("Cadeia com configuracao diferente.")
      X_longo<-cadeia$X
      log_msg("Cadeia recuperada: replica ",rep)
    } else {
      log_msg("Simulando GL: replica ",rep,"; T maximo = ",max(TAMANHOS))
      set.seed(SEED_GL+rep)
      X_longo<-simular_gl(W_sim,max(TAMANHOS),BURN_IN)
      colnames(X_longo)<-IDS
      salvar(list(config=cfg,X=X_longo),arquivo_X)
    }
    stopifnot(identical(dim(X_longo),c(max(TAMANHOS),N)),
              identical(colnames(X_longo),IDS),all(X_longo %in% 0:1))
    
    for (tt in TAMANHOS) {
      pasta_T<-file.path(pasta_rep,paste0("T",tt))
      dir.create(pasta_T,showWarnings=FALSE)
      X<-X_longo[seq_len(tt),,drop=FALSE]
      for (i in seq_len(N)) {
        arquivo<-file.path(pasta_T,sprintf("alvo_%03d.rds",i))
        if (file.exists(arquivo)) {
          obj<-readRDS(arquivo)
          if (!identical(obj$config,cfg) || !identical(obj$i,i) ||
              !identical(obj$T,tt) || !identical(obj$replica,rep)) {
            stop("Checkpoint incompativel: ",arquivo)
          }
          log_msg("T=",tt," | alvo=",i,"/",N," | recuperado")
          next
        }
        log_msg("T=",tt," | replica=",rep," | alvo=",i,"/",N,
                " | ID=",IDS[i]," | iniciando")
        inicio<-proc.time()[["elapsed"]]
        set.seed(SEED_SUR+10000L*rep+i)
        selecao<-selecionar_alvo(X,i,B,lambdas=LAMBDA,alfas=ALFAS)
        segundos<-proc.time()[["elapsed"]]-inicio
        met<-avaliar(selecao,i,tt,rep,segundos)
        salvar(list(config=cfg,T=tt,replica=rep,i=i,selecao=selecao,
                    metricas=met,segundos=segundos),arquivo)
        atualizar()
        log_msg("T=",tt," | alvo=",i," | salvo em ",round(segundos/60,2)," min")
        print(met[,c("alfa","TP","FP","FN","TN","exato")],row.names=FALSE)
        flush.console()
      }
      s<-atualizar()
      print(s,row.names=FALSE)
      flush.console()
    }
  }
  log_msg("CONCLUIDO. Resultados: ",SAIDA)
}

executar_roi_corrigida()
