# =============================================================================
# ENTREGAVEL B -- ESTUDO DE SENSIBILIDADE DA REGRA DE PARADA
# SCRIPT AUTONOMO -- FUNCOES-NUCLEO INCORPORADAS
#
# Este arquivo NAO depende de:
#   funcoes_core_regra_parada.R
#
# As funcoes de simulacao, informacao mutua, reconstrucao de U,
# selecao forward, regra de parada, topologia e geracao dos cenarios
# estao incorporadas neste proprio arquivo.
# =============================================================================

rm(list = ls())

# =============================================================================
# FUNCOES-NUCLEO -- incorporadas neste script autonomo
#
# Este arquivo e' fonte comum para os dois entregaveis:
#   - entregavel_A_regra_de_parada.R   (pedido do orientador)
#   - entregavel_B_sensibilidade.R     (estudo de sensibilidade, iniciativa propria)
#
# Reaproveita a base MAIS RAPIDA ja validada no pipeline da Tarefa 21 (F):
#   simulador vetorizado (runif + comparacao, em vez de rbinom) e informacao
#   mutua por contagens (tabulate, em vez de table+paste). A logica matematica
#   e' identica a`s Tarefas 19/20 -- so a implementacao e' mais rapida.
#
# NOVIDADE deste arquivo (nao existia nas Tarefas 19-21): o simulador agora
# tambem devolve os disparos BRUTOS X de TODOS os neuronios (nao so do alvo),
# porque a regra de parada (Secao 3 do manuscrito) precisa deslocar
# circularmente o trem de disparos de cada candidato para gerar os surrogates,
# e reconstruir U^(j->alvo) a partir do disparo deslocado + idade real do alvo.
# Isso e' feito por soma cumulativa (reconstruir_U_shift), O(T) por candidato,
# SEM precisar re-simular a rede inteira -- e' o que torna o teste de B
# surrogates por passo computacionalmente viavel.
#
# Validacao da reconstrucao (feita durante o desenvolvimento, nao repetida
# aqui): reconstruir_U_shift(X_cand_sem_deslocar, A_alvo) reproduz EXATAMENTE
# (diferenca maxima = 0) a coluna U_alvo produzida pela simulacao original.
# =============================================================================

# -----------------------------------------------------------------------
# 1. SIMULADOR GL (vetorizado, identico em logica a`s Tarefas 19-21)
# -----------------------------------------------------------------------
# Devolve, apos remover o burn-in:
#   X       -- matriz T x N de disparos brutos de TODOS os neuronios
#   U_alvo  -- matriz T x N, coluna j = U_t^(j -> alvo) (so' faz sentido nas
#              colunas dos candidatos; a coluna do proprio alvo nao e' usada)
#   A_alvo  -- vetor T, idade do alvo (A_t^(alvo) = t - L_t^(alvo) - 1), SEM
#              truncar em A_cap (a truncagem, quando necessaria, e' aplicada
#              depois, so' nas quantidades diagnosticas condicionadas em idade)
simulate_GL <- function(W, T_total, burn_in, mu, alvo, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(W); T_full <- T_total + burn_in
  Unif <- matrix(runif(T_full * N), nrow = T_full, ncol = N)
  X      <- matrix(0L, T_full, N)
  U_alvo <- matrix(0,  T_full, N)
  A_alvo <- integer(T_full)
  U <- matrix(0, N, N); A <- integer(N)
  for (t in seq_len(T_full)) {
    p <- 1 / (1 + exp(-(mu + rowSums(W * U))))
    f <- as.numeric(Unif[t, ] < p)
    X[t, ]      <- f
    U_alvo[t, ] <- U[alvo, ]
    A_alvo[t]   <- A[alvo]
    U <- (U + matrix(f, N, N, byrow = TRUE)) * (1 - f)
    A <- (A + 1) * (1 - f)
  }
  idx <- (burn_in + 1):T_full
  list(X = X[idx, , drop = FALSE], U_alvo = U_alvo[idx, , drop = FALSE], A_alvo = A_alvo[idx])
}

# -----------------------------------------------------------------------
# 2. INFORMACAO MUTUA EMPIRICA (marginal se z = NULL, condicional c.c.)
#    Identica, em logica e em implementacao, a` versao do pipeline F.
# -----------------------------------------------------------------------
info_mutua <- function(x, y, z = NULL) {
  n <- length(x); xc <- match(x, sort(unique(x))); yc <- match(y, sort(unique(y)))
  kx <- max(xc); ky <- max(yc)
  if (is.null(z)) {
    zc <- rep(1L, n)
  } else {
    z <- as.matrix(z)
    zc <- match(z[, 1], sort(unique(z[, 1])))
    if (ncol(z) > 1) for (l in 2:ncol(z)) {
      cl  <- match(z[, l], sort(unique(z[, l])))
      comb <- (zc - 1) * max(cl) + cl
      zc  <- match(comb, sort(unique(comb)))
    }
  }
  kz <- max(zc)
  # as.double(): tabulate() devolve inteiro (32 bits); em T grande (>~46000
  # observacoes concentradas numa mesma celula, o que ocorre por ex. em
  # T_total = 100000 nos cenarios GL) o PRODUTO de duas contagens pode
  # ultrapassar .Machine$integer.max e estourar silenciosamente (virando NA
  # com um warning de "integer overflow"). Forcar double aqui evita isso.
  cxyz <- as.double(tabulate((zc - 1) * kx * ky + (xc - 1) * ky + yc, kz * kx * ky))
  cxz  <- as.double(tabulate((zc - 1) * kx + xc, kz * kx))
  cyz  <- as.double(tabulate((zc - 1) * ky + yc, kz * ky))
  cz   <- as.double(tabulate(zc, kz))
  i  <- which(cxyz > 0); zi <- (i - 1) %/% (kx * ky) + 1; r <- (i - 1) %% (kx * ky)
  xi <- r %/% ky + 1; yi <- r %% ky + 1
  sum(cxyz[i] / n * log(cxyz[i] * cz[zi] / (cxz[(zi - 1) * kx + xi] * cyz[(zi - 1) * ky + yi])))
}

# -----------------------------------------------------------------------
# 3. RECONSTRUCAO DE U^(j->alvo) A PARTIR DO DISPARO BRUTO (via soma cumulativa)
#    U_t^(j->alvo) = # disparos de j em (L_t^(alvo), t-1], L_t^(alvo) = t-A_t-1
# -----------------------------------------------------------------------
reconstruir_U_shift <- function(X_cand, A_alvo) {
  Tn <- length(X_cand)
  L  <- pmax(seq_len(Tn) - A_alvo - 1, 0)
  Ccum  <- c(0, cumsum(X_cand))
  upper <- pmax(seq_len(Tn) - 1, 0)
  Ccum[upper + 1] - Ccum[L + 1]
}

# Deslocamento circular de um vetor por 'lag' posicoes (usado para os
# surrogates da regra de parada, Secao 3, passo 4 do manuscrito).
deslocar_circular <- function(x, lag) {
  n <- length(x)
  if (lag %% n == 0) return(x)
  x[((seq_len(n) - lag - 1) %% n) + 1]
}

# -----------------------------------------------------------------------
# 4. SELECAO FORWARD (SEM regra de parada) -- ranking completo, como nas
#    Tarefas 19-21. Usada aqui so' para comparacao (coluna "ranking normal"
#    nos dois entregaveis), reaproveitando o tamanho verdadeiro |V_i|.
# -----------------------------------------------------------------------
selecao_forward <- function(X_alvo, U_candidatos, lambda = 1) {
  restantes <- names(U_candidatos); S <- character(0)
  while (length(restantes) > 0) {
    Q <- sapply(restantes, function(j) {
      Rel <- info_mutua(U_candidatos[[j]], X_alvo)
      if (length(S) == 0) return(Rel)
      Inc <- info_mutua(U_candidatos[[j]], X_alvo, do.call(cbind, U_candidatos[S]))
      (1 - lambda) * Rel + lambda * Inc
    })
    prox <- names(which.max(Q)); S <- c(S, prox); restantes <- setdiff(restantes, prox)
  }
  S
}

# -----------------------------------------------------------------------
# 5. SELECAO FORWARD **COM REGRA DE PARADA** (Secao 3 do manuscrito)
#
#   Passos 1-9 da Secao 3, generalizados para lambda arbitrario (o texto do
#   artigo usa apenas a informacao incremental pura, o que equivale a
#   lambda = 1 na formula Q = (1-lambda)*Rel + lambda*Inc; deixamos lambda
#   configuravel para o estudo de sensibilidade do Entregavel B).
#
# Argumentos:
#   X_alvo, A_alvo        -- disparos e idade do neuronio alvo (saida do simulador)
#   U_candidatos          -- lista nomeada de vetores U^(j->alvo), UM POR
#                            candidato, computados a partir dos dados reais
#                            (nao deslocados)
#   X_candidatos          -- lista nomeada dos disparos BRUTOS de cada
#                            candidato (usados para gerar os surrogates)
#   lambda                -- peso de Inc na formula Q (1 = artigo original)
#   B                     -- numero de replicas surrogate por passo
#   alpha                 -- nivel de significancia do teste de maximo
#   resample_cada_passo   -- se TRUE, sorteia um novo lag para cada candidato
#                            a cada passo m (interpretacao mais literal do
#                            texto); se FALSE, sorteia um unico lag por
#                            candidato no inicio e reutiliza em todos os passos
#
# Retorna:
#   V_hat      -- vizinhanca estimada (vetor de rotulos, na ordem de selecao)
#   historico  -- data.frame com uma linha por passo: M_obs, limiar, se parou,
#                 quantos candidatos restavam e qual seria o proximo escolhido
# -----------------------------------------------------------------------
selecao_forward_com_parada <- function(X_alvo, A_alvo, U_candidatos, X_candidatos,
                                       lambda = 1, B = 200, alpha = 0.05,
                                       resample_cada_passo = TRUE) {
  candidatos <- names(U_candidatos)
  Tn <- length(X_alvo)
  S  <- character(0)
  restantes <- candidatos
  historico <- vector("list", length(candidatos))
  passo <- 0
  
  score_j <- function(Uj, S_atual) {
    Rel <- info_mutua(Uj, X_alvo)
    if (length(S_atual) == 0) return(Rel)          # Inc(j | vazio) = Rel(j), eq. apos (1.13)
    Inc <- info_mutua(Uj, X_alvo, do.call(cbind, U_candidatos[S_atual]))
    (1 - lambda) * Rel + lambda * Inc
  }
  
  lags_fixos <- NULL
  if (!resample_cada_passo) {
    lags_fixos <- setNames(sample.int(Tn - 1, length(candidatos), replace = TRUE), candidatos)
  }
  
  while (length(restantes) > 0) {
    passo <- passo + 1
    
    # --- passo 2-3 da Secao 3: escore observado e maximo observado ---
    scores_obs <- sapply(restantes, function(j) score_j(U_candidatos[[j]], S))
    M_obs <- max(scores_obs)
    
    # --- passos 4-5: B replicas surrogate (shift circular dos candidatos
    #     restantes; alvo e candidatos ja selecionados ficam intocados) ---
    M_surr <- numeric(B)
    for (b in seq_len(B)) {
      scores_surr <- sapply(restantes, function(j) {
        lag <- if (resample_cada_passo) sample.int(Tn - 1, 1) else lags_fixos[[j]]
        Xs  <- deslocar_circular(X_candidatos[[j]], lag)
        Us  <- reconstruir_U_shift(Xs, A_alvo)
        score_j(Us, S)
      })
      M_surr[b] <- max(scores_surr)
    }
    
    # --- passo 6: limiar nulo (quantil empirico 1-alpha dos maximos surrogate) ---
    limiar <- as.numeric(quantile(M_surr, probs = 1 - alpha, type = 7, names = FALSE))
    parou  <- M_obs <= limiar
    
    historico[[passo]] <- data.frame(
      passo = passo, M_obs = M_obs, limiar = limiar, parou = parou,
      n_restantes = length(restantes), candidato_top = names(which.max(scores_obs)),
      stringsAsFactors = FALSE
    )
    
    # --- passo 7: aplica o criterio de parada ---
    if (parou) break
    
    # --- passo 8: seleciona o melhor candidato restante e atualiza S ---
    prox <- names(which.max(scores_obs))
    S <- c(S, prox)
    restantes <- setdiff(restantes, prox)
  }
  
  list(V_hat = S, historico = do.call(rbind, historico[seq_len(passo)]))
}

# -----------------------------------------------------------------------
# 6. UTILITARIOS DE TOPOLOGIA (identicos, em logica, a`s Tarefas 19-21)
# -----------------------------------------------------------------------
construir_W <- function(N, arestas) {
  W <- matrix(0, N, N)
  for (e in arestas) W[e["to"], e["from"]] <- e["peso"]
  W
}

pais_verdadeiros <- function(arestas, alvo_idx, rotulos) {
  if (length(arestas) == 0) return(character(0))
  origens <- unlist(lapply(arestas, function(e) if (e["to"] == alvo_idx) e["from"] else NA))
  origens <- sort(unique(origens[!is.na(origens)]))
  if (length(origens) == 0) return(character(0))
  rotulos[origens]
}

resumir_estatisticas <- function(x) {
  c(media = mean(x), dp = sd(x), mediana = median(x),
    q05 = unname(quantile(x, 0.05)), q95 = unname(quantile(x, 0.95)))
}

# -----------------------------------------------------------------------
# 7. DEFINICAO UNIFICADA DOS CENARIOS (GL, B2, F) -- MESMOS parametros
#    (topologias, pesos, T, burn-in, n_mc) usados nas Tarefas 19-21, para
#    que os resultados COM regra de parada sejam diretamente comparaveis
#    aos resultados SEM regra de parada ja obtidos.
#
#    Cada job e' uma lista com: id, familia ("GL"/"B2"/"F"), cenario, N,
#    rotulos, alvo_idx, arestas (lista de c(from=,to=,peso=)), T_total,
#    burn_in, n_mc, seed_base, parametros (para registro/legenda).
# -----------------------------------------------------------------------

p0            <- 0.02
mu_intercepto <- log(p0 / (1 - p0))
A_cap         <- 10

rotulos_3 <- c("j", "k", "i")
rotulos_5 <- c("ell", "j", "k", "m", "i")

# --- 7.1 GL-0 a GL-4 (Tarefa 19 + extensao da Tarefa 19b: alvo generico) ---
# Topologias-base (mesmas do doc "Extensao das simulacoes diagnosticas GL").
.cenarios_GL_base <- list(
  "GL-0" = list(arestas = list()),
  "GL-1" = list(arestas = list(c(from = 2, to = 3, peso = NA))),                       # k -> i
  "GL-2" = list(arestas = list(c(from = 1, to = 3, peso = NA))),                       # j -> i
  "GL-3" = list(arestas = list(c(from = 1, to = 2, peso = NA), c(from = 2, to = 3, peso = NA))),  # j->k->i
  "GL-4" = list(arestas = list(c(from = 2, to = 1, peso = NA), c(from = 2, to = 3, peso = NA)))   # k->j, k->i
)

gerar_jobs_GL <- function(T_total = 100000, burn_in = 2000, n_mc = 10,
                          omega_vals = c(1.5, 4, 5), alvos = c(k = 2, i = 3)) {
  jobs <- list(); contador <- 0
  for (nome_cenario in names(.cenarios_GL_base)) {
    arestas_base <- .cenarios_GL_base[[nome_cenario]]$arestas
    for (omega in omega_vals) {
      arestas_omega <- lapply(arestas_base, function(e) { e["peso"] <- omega; e })
      for (nome_alvo in names(alvos)) {
        contador <- contador + 1
        alvo_idx <- alvos[[nome_alvo]]
        jobs[[length(jobs) + 1]] <- list(
          id = sprintf("GL_%s_omega%s_alvo%s", nome_cenario, omega, nome_alvo),
          familia = "GL", cenario = nome_cenario, N = 3, rotulos = rotulos_3,
          alvo_idx = alvo_idx, arestas = arestas_omega,
          T_total = T_total, burn_in = burn_in, n_mc = n_mc,
          seed_base = 10000 * contador + 1,
          parametros = list(omega = omega, alvo = nome_alvo)
        )
      }
    }
  }
  jobs
}

# --- 7.2 B2-0 a B2-4 (Tarefa 20) ---
gerar_jobs_B2 <- function(n_mc = 100, T_vals = c(5000, 10000, 40000)) {
  definir_burn_in <- function(T_total) max(500, min(2000, round(0.1 * T_total)))
  omega_vals_B2_1 <- c(0.5, 1.5, 5)
  ab_vals_B2_2    <- list(c(alpha = 1.5, beta = 5), c(alpha = 4, beta = 4), c(alpha = 5, beta = 1.5))
  ab_vals_B2_3    <- list(c(alpha = 1.5, beta = 1.5), c(alpha = 5, beta = 1.5), c(alpha = 1.5, beta = 5))
  abg_vals_B2_4   <- list(c(alpha = 4, beta = 4, gamma = 4),
                          c(alpha = 4, beta = 4, gamma = 1.5),
                          c(alpha = 5, beta = 1.5, gamma = 4))
  jobs <- list(); contador <- 0
  
  for (T_total in T_vals) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("B2-0_T%d", T_total), familia = "B2", cenario = "B2-0",
      N = 3, rotulos = rotulos_3, alvo_idx = 3, arestas = list(),
      T_total = T_total, burn_in = definir_burn_in(T_total), n_mc = n_mc,
      seed_base = 10000 * contador + 1, parametros = list())
  }
  for (T_total in T_vals) for (omega in omega_vals_B2_1) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("B2-1_T%d_omega%s", T_total, omega), familia = "B2", cenario = "B2-1",
      N = 3, rotulos = rotulos_3, alvo_idx = 3,
      arestas = list(c(from = 2, to = 3, peso = omega)),
      T_total = T_total, burn_in = definir_burn_in(T_total), n_mc = n_mc,
      seed_base = 10000 * contador + 1, parametros = list(omega = omega))
  }
  for (T_total in T_vals) for (par in ab_vals_B2_2) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("B2-2_T%d_a%s_b%s", T_total, par["alpha"], par["beta"]),
      familia = "B2", cenario = "B2-2", N = 3, rotulos = rotulos_3, alvo_idx = 3,
      arestas = list(c(from = 1, to = 2, peso = unname(par["alpha"])),
                     c(from = 2, to = 3, peso = unname(par["beta"]))),
      T_total = T_total, burn_in = definir_burn_in(T_total), n_mc = n_mc,
      seed_base = 10000 * contador + 1,
      parametros = list(alpha = unname(par["alpha"]), beta = unname(par["beta"])))
  }
  for (T_total in T_vals) for (par in ab_vals_B2_3) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("B2-3_T%d_a%s_b%s", T_total, par["alpha"], par["beta"]),
      familia = "B2", cenario = "B2-3", N = 3, rotulos = rotulos_3, alvo_idx = 3,
      arestas = list(c(from = 2, to = 1, peso = unname(par["alpha"])),
                     c(from = 2, to = 3, peso = unname(par["beta"]))),
      T_total = T_total, burn_in = definir_burn_in(T_total), n_mc = n_mc,
      seed_base = 10000 * contador + 1,
      parametros = list(alpha = unname(par["alpha"]), beta = unname(par["beta"])))
  }
  for (T_total in T_vals) for (par in abg_vals_B2_4) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("B2-4_T%d_a%s_b%s_g%s", T_total, par["alpha"], par["beta"], par["gamma"]),
      familia = "B2", cenario = "B2-4", N = 5, rotulos = rotulos_5, alvo_idx = 5,
      arestas = list(c(from = 1, to = 2, peso = unname(par["alpha"])),
                     c(from = 2, to = 5, peso = unname(par["beta"])),
                     c(from = 3, to = 5, peso = unname(par["gamma"]))),
      T_total = T_total, burn_in = definir_burn_in(T_total), n_mc = n_mc,
      seed_base = 10000 * contador + 1,
      parametros = list(alpha = unname(par["alpha"]), beta = unname(par["beta"]),
                        gamma = unname(par["gamma"])))
  }
  jobs
}

# --- 7.3 F-0 a F-2 (Tarefa 21) ---
gerar_jobs_F <- function(n_mc = 100, T_vals = c(5000, 10000, 40000), burn_in_fixo = 2000) {
  ab_F1  <- list(c(alpha = 5, beta = 1), c(alpha = 10, beta = 1), c(alpha = 10, beta = 0.5))
  abg_F2 <- list(c(alpha = 5, beta = 0.5, gamma = 4), c(alpha = 5, beta = 4, gamma = 0.5),
                 c(alpha = 5, beta = 0.5, gamma = 0.5))
  jobs <- list(); contador <- 0
  
  for (T_total in T_vals) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("F-0-null_T%d", T_total), familia = "F", cenario = "F-0-null",
      N = 3, rotulos = rotulos_3, alvo_idx = 3, arestas = list(),
      T_total = T_total, burn_in = burn_in_fixo, n_mc = n_mc,
      seed_base = 100000 + 1000 * contador, parametros = list())
  }
  for (T_total in T_vals) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("F-0-dir_T%d", T_total), familia = "F", cenario = "F-0-dir",
      N = 3, rotulos = rotulos_3, alvo_idx = 3,
      arestas = list(c(from = 2, to = 3, peso = 0.5)),
      T_total = T_total, burn_in = burn_in_fixo, n_mc = n_mc,
      seed_base = 100000 + 1000 * contador, parametros = list(omega = 0.5))
  }
  for (T_total in T_vals) for (par in ab_F1) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("F-1_T%d_a%s_b%s", T_total, par["alpha"], par["beta"]),
      familia = "F", cenario = "F-1", N = 3, rotulos = rotulos_3, alvo_idx = 3,
      arestas = list(c(from = 1, to = 2, peso = unname(par["alpha"])),
                     c(from = 2, to = 3, peso = unname(par["beta"]))),
      T_total = T_total, burn_in = burn_in_fixo, n_mc = n_mc,
      seed_base = 100000 + 1000 * contador,
      parametros = list(alpha = unname(par["alpha"]), beta = unname(par["beta"])))
  }
  for (T_total in T_vals) for (par in abg_F2) {
    contador <- contador + 1
    jobs[[length(jobs) + 1]] <- list(
      id = sprintf("F-2_T%d_a%s_b%s_g%s", T_total, par["alpha"], par["beta"], par["gamma"]),
      familia = "F", cenario = "F-2", N = 5, rotulos = rotulos_5, alvo_idx = 5,
      arestas = list(c(from = 1, to = 2, peso = unname(par["alpha"])),
                     c(from = 2, to = 5, peso = unname(par["beta"])),
                     c(from = 3, to = 5, peso = unname(par["gamma"]))),
      T_total = T_total, burn_in = burn_in_fixo, n_mc = n_mc,
      seed_base = 100000 + 1000 * contador,
      parametros = list(alpha = unname(par["alpha"]), beta = unname(par["beta"]),
                        gamma = unname(par["gamma"])))
  }
  jobs
}

gerar_todos_os_jobs <- function(n_mc_GL = 10, n_mc_B2 = 100, n_mc_F = 100,
                                T_total_GL = 100000, burn_in_GL = 2000) {
  jobs <- c(gerar_jobs_GL(T_total = T_total_GL, burn_in = burn_in_GL, n_mc = n_mc_GL),
            gerar_jobs_B2(n_mc = n_mc_B2),
            gerar_jobs_F(n_mc = n_mc_F))
  for (idx in seq_along(jobs)) jobs[[idx]]$job_idx <- idx
  jobs
}

# -----------------------------------------------------------------------
# 8. RODAR UMA REPLICA DE UM JOB: devolve tudo que os dois entregaveis
#    precisam (dados simulados + candidatos), para nao repetir a simulacao
#    entre diferentes combinacoes de hiperparametro da regra de parada.
# -----------------------------------------------------------------------
simular_replica_job <- function(job, seed) {
  W <- construir_W(job$N, job$arestas)
  candidatos_idx <- setdiff(seq_len(job$N), job$alvo_idx)
  rot_cand <- job$rotulos[candidatos_idx]
  
  s <- simulate_GL(W, job$T_total, job$burn_in, mu_intercepto, job$alvo_idx, seed = seed)
  
  list(
    X_alvo       = s$X[, job$alvo_idx],
    A_alvo       = s$A_alvo,
    U_candidatos = setNames(lapply(candidatos_idx, function(idx) s$U_alvo[, idx]), rot_cand),
    X_candidatos = setNames(lapply(candidatos_idx, function(idx) s$X[, idx]), rot_cand),
    pais_verdadeiros = pais_verdadeiros(job$arestas, job$alvo_idx, job$rotulos)
  )
}

# =============================================================================
# ENTREGAVEL B -- Estudo de sensibilidade da regra de parada (iniciativa propria,
# para discutir com o orientador; usa OS MESMOS 93 cenarios do Entregavel A)
#
# O QUE ESTE SCRIPT FAZ:
#   Para cada uma das 93 configuracoes ja usadas nas Tarefas 19-21 (e no
#   Entregavel A), varia um GRID de hiperparametros da regra de parada e mede
#   o impacto de cada um no resultado (tamanho de V_hat, taxa de recuperacao
#   exata, falsos positivos/negativos, passo de parada, custo computacional):
#
#     B (num. de surrogates por passo) .... grade em GRADE_B
#     alpha (nivel do teste de maximo) .... grade em GRADE_ALPHA
#     lambda (peso de Inc na formula Q) ... grade em GRADE_LAMBDA
#     estrategia de shift ................. resorteado a cada passo vs. fixo
#
#   IMPORTANTE -- OTIMIZACAO DE CUSTO (nao muda a logica, so' a ordem das
#   contas): para uma mesma replica de um mesmo cenario, os dados simulados
#   (X, U, A) NAO dependem dos hiperparametros da regra de parada. Por isso,
#   cada replica e' simulada UMA UNICA VEZ e reutilizada para TODAS as
#   combinacoes do grid daquele cenario -- em vez de re-simular a rede a cada
#   combinacao (o que seria = numero_de_combinacoes x mais lento sem
#   nenhum ganho de informacao).
#
# CUSTO COMPUTACIONAL -- LEIA ANTES DE RODAR A BATERIA COMPLETA:
#   Mesmo com a otimizacao acima, o custo da REGRA DE PARADA em si (nao da
#   simulacao) escala linearmente com o tamanho do grid. Com a grade padrao
#   deste script (4 x 3 x 3 x 2 = 72 combinacoes) x 93 configuracoes x
#   N_MC_SENSIBILIDADE replicas, o tempo total estimado (ver benchmark
#   impresso automaticamente no console ao iniciar RODAR_BATERIA_COMPLETA)
#   pode chegar a DEZENAS DE HORAS em 1 nucleo. Por isso:
#     (a) N_MC_SENSIBILIDADE por padrao e' MENOR que o n_mc usado no
#         Entregavel A (default: 20 em vez de ate 100) -- isso e' uma
#         escolha documentada aqui, nao pedida explicitamente por voce;
#         ajuste se quiser mais precisao ao custo de mais tempo.
#     (b) a execucao e' idempotente por CONFIGURACAO (nao por combinacao de
#         hiperparametro): cada configuracao roda o grid inteiro de uma vez
#         e salva um unico .rds; internamente interrompida, retomar so'
#         refaz configuracoes ainda nao salvas.
#     (c) use N_CORES > 1 (Linux/Mac) para paralelizar entre replicas.
#   Se o tempo estimado for inviavel, reduza GRADE_B/GRADE_ALPHA/GRADE_LAMBDA
#   ou N_MC_SENSIBILIDADE no Painel de Controle antes de rodar a bateria toda.
# =============================================================================


# =============================================================================
# SECAO 0 -- PAINEL DE CONTROLE
# =============================================================================
RODAR_TESTE_DE_SANIDADE <- FALSE    # grid pequeno, T pequeno, minutos
RODAR_BATERIA_COMPLETA  <- TRUE   # grid completo x 93 configs -- LEIA O AVISO ACIMA
FORCAR_RECALCULO        <- FALSE
RODAR_AGREGACAO         <- TRUE

PASTA_RESULTADOS     <- "resultados_entregavel_B"
PASTA_SAIDA_AGREGADA <- "saida_agregada_entregavel_B"

SEMENTE_GLOBAL <- 17
N_CORES        <- max(1, parallel::detectCores() - 1)

# --- grade de hiperparametros da regra de parada (ESTUDO DE SENSIBILIDADE) ---
GRADE_B       <- c(50, 100, 200, 500)
GRADE_ALPHA   <- c(0.01, 0.05, 0.10)
GRADE_LAMBDA  <- c(0, 0.5, 1)          # 1 = regra literal do manuscrito (Secao 3)
GRADE_RESAMPLE <- c(TRUE, FALSE)       # TRUE = resorteia lag a cada passo m

# --- escala das simulacoes (REDUZIDA em relacao ao Entregavel A -- ver aviso) ---
N_MC_SENSIBILIDADE <- 20   # Entregavel A usa 10 (GL) / 100 (B2, F); aqui, 20 para todos
T_TOTAL_GL <- 100000
BURN_IN_GL <- 2000

set.seed(SEMENTE_GLOBAL)


# =============================================================================
# SECAO 1 -- GRADE DE HIPERPARAMETROS
# =============================================================================
gerar_grade_hp <- function(B_vals = GRADE_B, alpha_vals = GRADE_ALPHA,
                           lambda_vals = GRADE_LAMBDA, resample_vals = GRADE_RESAMPLE) {
  grade <- expand.grid(B = B_vals, alpha = alpha_vals, lambda = lambda_vals,
                       resample = resample_vals, KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)
  grade$hp_id <- sprintf("B%d_a%s_l%s_r%s", grade$B, grade$alpha, grade$lambda, grade$resample)
  grade
}


# =============================================================================
# SECAO 2 -- EXECUCAO DE UM JOB (todas as replicas x toda a grade de hp)
# =============================================================================
executar_job_B <- function(job, grade_hp, n_mc, n_cores = 1) {
  pais_verd <- pais_verdadeiros(job$arestas, job$alvo_idx, job$rotulos)
  
  rodar_uma_replica <- function(r) {
    seed_r <- job$seed_base + r
    d <- simular_replica_job(job, seed = seed_r)   # simulado UMA vez, reusado em toda a grade
    
    linhas_hp <- lapply(seq_len(nrow(grade_hp)), function(g) {
      hp <- grade_hp[g, ]
      t0 <- Sys.time()
      res <- selecao_forward_com_parada(
        d$X_alvo, d$A_alvo, d$U_candidatos, d$X_candidatos,
        lambda = hp$lambda, B = hp$B, alpha = hp$alpha,
        resample_cada_passo = hp$resample
      )
      tempo <- as.numeric(Sys.time() - t0, units = "secs")
      
      V_hat <- res$V_hat
      data.frame(
        replicacao = r, hp_id = hp$hp_id, B = hp$B, alpha = hp$alpha,
        lambda = hp$lambda, resample = hp$resample,
        tamanho_Vhat = length(V_hat),
        exato = setequal(V_hat, pais_verd),
        falso_positivo = length(setdiff(V_hat, pais_verd)) > 0,
        falso_negativo = length(setdiff(pais_verd, V_hat)) > 0,
        passo_parada = nrow(res$historico),
        tempo_seg = tempo,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, linhas_hp)
  }
  
  linhas <- if (n_cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_len(n_mc), rodar_uma_replica, mc.cores = n_cores)
  } else {
    lapply(seq_len(n_mc), rodar_uma_replica)
  }
  
  resultado <- do.call(rbind, linhas)
  attr(resultado, "pais_verdadeiros") <- pais_verd
  attr(resultado, "job") <- job
  resultado
}


# =============================================================================
# SECAO 3 -- BENCHMARK RAPIDO (estima tempo total antes de rodar tudo)
# =============================================================================
estimar_tempo_total <- function(jobs, grade_hp, n_mc) {
  # usa o job MAIS CARO (maior T_total) como referencia pessimista
  job_ref <- jobs[[which.max(sapply(jobs, function(j) j$T_total))]]
  hp_ref  <- grade_hp[which.max(grade_hp$B), ]   # combinacao mais cara do grid
  
  d <- simular_replica_job(job_ref, seed = 999999)
  t0 <- Sys.time()
  invisible(selecao_forward_com_parada(d$X_alvo, d$A_alvo, d$U_candidatos, d$X_candidatos,
                                       lambda = hp_ref$lambda, B = hp_ref$B, alpha = hp_ref$alpha,
                                       resample_cada_passo = hp_ref$resample))
  t_max <- as.numeric(Sys.time() - t0, units = "secs")
  
  # estimativa grosseira: usa a MEDIA do grid de B como proxy do custo tipico
  # (custo escala ~linear em B), aplicada a todas as configuracoes/replicas
  fator_B_medio <- mean(grade_hp$B) / max(grade_hp$B)
  custo_medio_estimado <- t_max * fator_B_medio
  
  total_execucoes <- length(jobs) * nrow(grade_hp) * n_mc
  tempo_serial_h  <- total_execucoes * custo_medio_estimado / 3600
  
  cat(sprintf("\n[benchmark] 1 execucao no pior caso (T=%d, B=%d): %.2fs\n",
              job_ref$T_total, hp_ref$B, t_max))
  cat(sprintf("[benchmark] execucoes totais previstas: %d configuracoes x %d combinacoes x %d replicas = %d\n",
              length(jobs), nrow(grade_hp), n_mc, total_execucoes))
  cat(sprintf("[benchmark] tempo estimado (1 nucleo, MUITO aproximado): %.1f horas\n", tempo_serial_h))
  cat(sprintf("[benchmark] com N_CORES=%d: ~%.1f horas\n\n", N_CORES, tempo_serial_h / N_CORES))
  invisible(tempo_serial_h)
}


# =============================================================================
# SECAO 4 -- TESTE DE SANIDADE (grid e T pequenos)
# =============================================================================
rodar_teste_de_sanidade <- function() {
  grade_pequena <- gerar_grade_hp(B_vals = c(30, 100), alpha_vals = c(0.05, 0.10),
                                  lambda_vals = c(0.5, 1), resample_vals = c(TRUE, FALSE))
  cat("Grade de teste:", nrow(grade_pequena), "combinacoes\n")
  
  jobs <- gerar_todos_os_jobs(n_mc_GL = 2, n_mc_B2 = 2, n_mc_F = 2, T_total_GL = 3000, burn_in_GL = 500)
  jobs <- lapply(jobs, function(j) { j$T_total <- min(j$T_total, 3000); j$burn_in <- min(j$burn_in, 500); j })
  escolhidos <- list(
    Filter(function(j) j$familia == "GL" && j$cenario == "GL-0", jobs)[[1]],
    Filter(function(j) j$familia == "B2" && j$cenario == "B2-1", jobs)[[1]],
    Filter(function(j) j$familia == "F"  && j$cenario == "F-2",  jobs)[[1]]
  )
  
  for (job in escolhidos) {
    t0 <- Sys.time()
    res <- executar_job_B(job, grade_pequena, n_mc = 3, n_cores = 1)
    dt <- as.numeric(Sys.time() - t0, units = "secs")
    cat(sprintf("\n=== %s (%.1fs, %d linhas) ===\n", job$id, dt, nrow(res)))
    agr <- aggregate(cbind(exato, tamanho_Vhat) ~ B + alpha, data = res, FUN = mean)
    print(agr)
  }
}


# =============================================================================
# SECAO 5 -- BATERIA COMPLETA (idempotente por configuracao)
# =============================================================================
rodar_bateria_completa <- function() {
  if (!dir.exists(PASTA_RESULTADOS)) dir.create(PASTA_RESULTADOS, recursive = TRUE)
  
  grade_hp <- gerar_grade_hp()
  jobs <- gerar_todos_os_jobs(n_mc_GL = N_MC_SENSIBILIDADE, n_mc_B2 = N_MC_SENSIBILIDADE,
                              n_mc_F = N_MC_SENSIBILIDADE, T_total_GL = T_TOTAL_GL, burn_in_GL = BURN_IN_GL)
  n_jobs <- length(jobs)
  
  estimar_tempo_total(jobs, grade_hp, N_MC_SENSIBILIDADE)
  
  cat("*******************************************************************\n")
  cat(sprintf("*** Grade de hiperparametros: %d combinacoes (B x alpha x lambda x shift)\n", nrow(grade_hp)))
  cat(sprintf("*** n_mc por configuracao: %d | configuracoes: %d\n", N_MC_SENSIBILIDADE, n_jobs))
  cat("*** Execucao idempotente por CONFIGURACAO (nao por combinacao de hp).\n")
  cat("*******************************************************************\n\n")
  
  for (idx in seq_along(jobs)) {
    job <- jobs[[idx]]
    arq <- file.path(PASTA_RESULTADOS, sprintf("resultado_%s.rds", job$id))
    if (file.exists(arq) && !FORCAR_RECALCULO) {
      cat(sprintf("[%03d/%d] %-40s -- ja existe, pulando\n", idx, n_jobs, job$id))
      next
    }
    t0 <- Sys.time()
    res <- executar_job_B(job, grade_hp, n_mc = N_MC_SENSIBILIDADE, n_cores = N_CORES)
    dt <- as.numeric(Sys.time() - t0, units = "secs")
    saveRDS(res, arq)
    cat(sprintf("[%03d/%d] %-40s concluido em %7.1fs (%d linhas)\n",
                idx, n_jobs, job$id, dt, nrow(res)))
  }
}


# =============================================================================
# SECAO 6 -- AGREGACAO: efeito marginal de cada hiperparametro
# =============================================================================
agregar_entregavel_B <- function() {
  if (!dir.exists(PASTA_SAIDA_AGREGADA)) dir.create(PASTA_SAIDA_AGREGADA, recursive = TRUE)
  arqs <- list.files(PASTA_RESULTADOS, pattern = "\\.rds$", full.names = TRUE)
  if (length(arqs) == 0) { warning("Nenhum .rds encontrado em ", PASTA_RESULTADOS); return(invisible(NULL)) }
  
  todos <- lapply(arqs, function(arq) {
    res <- readRDS(arq); job <- attr(res, "job")
    res$familia <- job$familia; res$cenario <- job$cenario; res$id <- job$id
    res$T_total <- job$T_total
    res$n_pais_verdadeiros <- length(attr(res, "pais_verdadeiros"))
    res
  })
  base <- do.call(rbind, todos)
  saveRDS(base, file.path(PASTA_SAIDA_AGREGADA, "entregavel_B_bruto_completo.rds"))
  
  # --- 1) tabela completa por configuracao x combinacao de hp (media sobre replicas) ---
  por_config_hp <- aggregate(
    cbind(tamanho_Vhat, exato, falso_positivo, falso_negativo, passo_parada, tempo_seg) ~
      familia + cenario + id + T_total + n_pais_verdadeiros + B + alpha + lambda + resample,
    data = base, FUN = mean)
  write.csv(por_config_hp, file.path(PASTA_SAIDA_AGREGADA, "tabela_completa_por_config_e_hp.csv"),
            row.names = FALSE)
  
  # --- 2) efeito MARGINAL de cada hiperparametro (media sobre tudo o mais) ---
  efeito_B <- aggregate(cbind(exato, tamanho_Vhat, falso_positivo, falso_negativo, tempo_seg) ~ B,
                        data = base, FUN = mean)
  efeito_alpha <- aggregate(cbind(exato, tamanho_Vhat, falso_positivo, falso_negativo, tempo_seg) ~ alpha,
                            data = base, FUN = mean)
  efeito_lambda <- aggregate(cbind(exato, tamanho_Vhat, falso_positivo, falso_negativo, tempo_seg) ~ lambda,
                             data = base, FUN = mean)
  efeito_resample <- aggregate(cbind(exato, tamanho_Vhat, falso_positivo, falso_negativo, tempo_seg) ~ resample,
                               data = base, FUN = mean)
  
  write.csv(efeito_B,        file.path(PASTA_SAIDA_AGREGADA, "efeito_B.csv"),        row.names = FALSE)
  write.csv(efeito_alpha,    file.path(PASTA_SAIDA_AGREGADA, "efeito_alpha.csv"),    row.names = FALSE)
  write.csv(efeito_lambda,   file.path(PASTA_SAIDA_AGREGADA, "efeito_lambda.csv"),   row.names = FALSE)
  write.csv(efeito_resample, file.path(PASTA_SAIDA_AGREGADA, "efeito_resample.csv"), row.names = FALSE)
  
  # --- 3) efeito marginal separado por familia de cenario (GL / B2 / F) ---
  efeito_B_por_familia <- aggregate(cbind(exato, tamanho_Vhat) ~ familia + B, data = base, FUN = mean)
  efeito_alpha_por_familia <- aggregate(cbind(exato, tamanho_Vhat) ~ familia + alpha, data = base, FUN = mean)
  efeito_lambda_por_familia <- aggregate(cbind(exato, tamanho_Vhat) ~ familia + lambda, data = base, FUN = mean)
  write.csv(efeito_B_por_familia, file.path(PASTA_SAIDA_AGREGADA, "efeito_B_por_familia.csv"), row.names = FALSE)
  write.csv(efeito_alpha_por_familia, file.path(PASTA_SAIDA_AGREGADA, "efeito_alpha_por_familia.csv"), row.names = FALSE)
  write.csv(efeito_lambda_por_familia, file.path(PASTA_SAIDA_AGREGADA, "efeito_lambda_por_familia.csv"), row.names = FALSE)
  
  cat(sprintf("\nAgregacao concluida: %d configuracoes, %d linhas brutas.\n", length(arqs), nrow(base)))
  cat("Arquivos gerados em", PASTA_SAIDA_AGREGADA, ":\n")
  cat("  tabela_completa_por_config_e_hp.csv  (uma linha por configuracao x combinacao de hp)\n")
  cat("  efeito_B.csv / efeito_alpha.csv / efeito_lambda.csv / efeito_resample.csv (efeito marginal)\n")
  cat("  efeito_*_por_familia.csv (efeito marginal, separado por GL / B2 / F)\n")
  
  list(por_config_hp = por_config_hp, efeito_B = efeito_B, efeito_alpha = efeito_alpha,
       efeito_lambda = efeito_lambda, efeito_resample = efeito_resample)
}


# =============================================================================
# EXECUCAO CONTROLADA (depende apenas das flags da Secao 0)
# =============================================================================
if (RODAR_TESTE_DE_SANIDADE) { cat(">>> Teste de sanidade\n"); rodar_teste_de_sanidade() }
if (RODAR_BATERIA_COMPLETA)  { rodar_bateria_completa() }
if (RODAR_AGREGACAO) {
  agregacao_B <- agregar_entregavel_B()
  if (!is.null(agregacao_B)) {
    cat("\n--- Efeito marginal de alpha ---\n"); print(agregacao_B$efeito_alpha)
    cat("\n--- Efeito marginal de lambda ---\n"); print(agregacao_B$efeito_lambda)
  }
}