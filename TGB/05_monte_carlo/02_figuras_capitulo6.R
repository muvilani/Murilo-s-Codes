
# ============================================================
# TCC - CAPITULO 6
# Matriz anatomica original e frequencia de selecao mRMR
# ROI: LOP_L_col_09_28
#
# Gera:
# 1. Matriz recuperada (PNG e PDF)
# 2. Comparacao original x recuperada (PNG e PDF)
# 3. Matrizes de frequencias e contagens (CSV)
#
# Nao executa novas simulacoes.
# ============================================================

local({
  
  # ----------------------------------------------------------
  # 1. CONFIGURACAO
  # ----------------------------------------------------------
  
  ROOT <- normalizePath(
    getwd(), winslash = "/", mustWork = TRUE
  )
  
  DIR_OUT <- file.path(ROOT, "resultados_MC_rapido")
  
  ARQ_W <- file.path(
    ROOT,
    "levantamento_regioes",
    "matrizes",
    "LOP_L_col_09_28",
    "W_regiao.rds"
  )
  
  if (!file.exists(ARQ_W)) {
    stop("Arquivo W_regiao.rds nao encontrado.")
  }
  
  W <- as.matrix(readRDS(ARQ_W))
  N <- nrow(W)
  
  stopifnot(
    N == 64L,
    ncol(W) == 64L,
    sum(W > 0) == 268L
  )
  
  Wbin <- 1.0 * (W > 0)
  diag(Wbin) <- NA_real_
  
  # ----------------------------------------------------------
  # 2. VALIDAR REPLICAS
  # ----------------------------------------------------------
  
  pastas <- file.path(
    DIR_OUT,
    sprintf("rep_%03d", 1:100)
  )
  
  validar_replica <- function(p, numero_replica) {
    
    if (!dir.exists(p)) return(FALSE)
    
    arquivos <- file.path(
      p, sprintf("alvo_%03d.rds", 1:N)
    )
    
    if (!all(file.exists(arquivos))) return(FALSE)
    
    for (i in seq_len(N)) {
      
      valido <- tryCatch({
        
        obj <- readRDS(arquivos[i])
        
        is.list(obj) &&
          length(obj$alvo) == 1L &&
          !is.na(obj$alvo) &&
          obj$alvo == i &&
          length(obj$replica) == 1L &&
          !is.na(obj$replica) &&
          obj$replica == numero_replica &&
          !is.null(obj$pais) &&
          is.numeric(obj$pais) &&
          !anyNA(obj$pais) &&
          all(obj$pais == as.integer(obj$pais)) &&
          !anyDuplicated(obj$pais) &&
          all(obj$pais %in% setdiff(1:N, i))
        
      }, error = function(e) FALSE)
      
      if (!valido) return(FALSE)
    }
    
    TRUE
  }
  
  validas <- vapply(
    seq_along(pastas),
    function(r) validar_replica(pastas[r], r),
    logical(1)
  )
  
  completas <- pastas[validas]
  R <- length(completas)
  
  cat("\n====================================\n")
  cat("REPLICAS VALIDAS:", R, "/100\n")
  
  if (R < 100) {
    cat("REPLICAS QUE PRECISAM DE VERIFICACAO:\n")
    print(which(!validas))
  }
  
  if (R == 0) {
    stop("Nenhuma replica completa e valida encontrada.")
  }
  
  # ----------------------------------------------------------
  # 3. MATRIZ DE FREQUENCIAS
  # ----------------------------------------------------------
  
  contagens <- matrix(
    0L, N, N, dimnames = dimnames(W)
  )
  
  for (p in completas) {
    
    for (i in seq_len(N)) {
      
      arq <- file.path(
        p, sprintf("alvo_%03d.rds", i)
      )
      
      obj <- readRDS(arq)
      
      if (length(obj$pais) > 0) {
        contagens[i, obj$pais] <-
          contagens[i, obj$pais] + 1L
      }
    }
  }
  
  freq <- 100 * contagens / R
  diag(freq) <- NA_real_
  
  write.csv(
    freq,
    file.path(DIR_OUT, "frequencia_selecao_64x64.csv"),
    na = ""
  )
  
  write.csv(
    contagens,
    file.path(DIR_OUT, "contagens_selecao_64x64.csv")
  )
  
  writeLines(
    c(
      paste("Replicas validas:", R),
      paste("Replicas excluidas:", paste(which(!validas),
                                         collapse = ", ")),
      "Frequencias em porcentagem de selecao.",
      "Matriz: linhas = alvos; colunas = origens.",
      "T = 40000; B = 199; alpha = 0.05.",
      "Lambda = 1; transformacao = log_maximo.",
      "Diagonal cinza = autoconexao excluida."
    ),
    file.path(DIR_OUT, "LEGENDA_METODOLOGICA.txt")
  )
  
  # ----------------------------------------------------------
  # 4. PALETAS
  # ----------------------------------------------------------
  
  # Tons de azul com contraste elevado
  cores_freq <- grDevices::colorRampPalette(
    c(
      "#FFFFFF",
      "#D5EDF6",
      "#8AC5E3",
      "#3182BD",
      "#08519C",
      "#001B44"
    )
  )(101)
  
  cor_conexao <- "#08306B"
  cor_ausencia <- "#FFFFFF"
  cor_diagonal <- "#E3E6EA"
  
  # ----------------------------------------------------------
  # 5. DESENHAR MATRIZ
  # ----------------------------------------------------------
  
  desenhar_matriz <- function(mat, titulo, original = FALSE) {
    
    graphics::par(
      mar = c(4.2, 4.5, 3.0, 1.0),
      family = "sans",
      bg = "white"
    )
    
    graphics::plot.new()
    
    graphics::plot.window(
      xlim = c(0.5, N + 0.5),
      ylim = c(0.5, N + 0.5),
      asp = 1
    )
    
    # Fundo cinza para a diagonal
    graphics::rect(
      0.5, 0.5, N + 0.5, N + 0.5,
      col = cor_diagonal,
      border = NA
    )
    
    if (original) {
      
      graphics::image(
        x = seq_len(N),
        y = seq_len(N),
        z = t(mat[N:1, , drop = FALSE]),
        col = c(cor_ausencia, cor_conexao),
        zlim = c(0, 1),
        add = TRUE,
        useRaster = TRUE
      )
      
    } else {
      
      graphics::image(
        x = seq_len(N),
        y = seq_len(N),
        z = t(mat[N:1, , drop = FALSE]) / 100,
        col = cores_freq,
        zlim = c(0, 1),
        add = TRUE,
        useRaster = TRUE
      )
    }
    
    ticks <- unique(c(
      1L, seq.int(8L, N, by = 8L), N
    ))
    
    graphics::axis(
      1,
      at = ticks,
      labels = ticks,
      cex.axis = 0.8,
      lwd = 0,
      lwd.ticks = 0
    )
    
    graphics::axis(
      2,
      at = N + 1L - ticks,
      labels = ticks,
      las = 1,
      cex.axis = 0.8,
      lwd = 0,
      lwd.ticks = 0
    )
    
    graphics::mtext(
      "Origem (pré-sináptico)",
      side = 1, line = 2.6, cex = 0.9
    )
    
    graphics::mtext(
      "Alvo (pós-sináptico)",
      side = 2, line = 2.8, cex = 0.9
    )
    
    graphics::title(
      main = titulo,
      cex.main = 1.1,
      font.main = 2
    )
    
    graphics::box(col = "#A1ACB5")
  }
  
  # ----------------------------------------------------------
  # 6. LEGENDAS HORIZONTAIS
  # ----------------------------------------------------------
  
  desenhar_legenda_freq <- function() {
    
    graphics::par(
      mar = c(2.5, 4.5, 0.1, 4.5),
      family = "sans"
    )
    
    graphics::plot.new()
    
    graphics::plot.window(
      xlim = c(0, 100),
      ylim = c(0, 1)
    )
    
    for (v in 0:99) {
      graphics::rect(
        v, 0.38, v + 1, 0.78,
        col = cores_freq[v + 1],
        border = NA
      )
    }
    
    graphics::axis(
      1,
      at = c(0, 25, 50, 75, 100),
      labels = paste0(c(0, 25, 50, 75, 100), "%"),
      cex.axis = 0.85,
      pos = 0.32,
      lwd = 0,
      lwd.ticks = 0
    )
    
    graphics::mtext(
      "Frequência de seleção",
      side = 1,
      line = 1.6,
      cex = 0.9
    )
  }
  
  desenhar_legenda_original <- function() {
    
    graphics::par(
      mar = c(2.5, 4.5, 0.1, 4.5),
      family = "sans"
    )
    
    graphics::plot.new()
    
    graphics::plot.window(
      xlim = c(0, 100),
      ylim = c(0, 1)
    )
    
    graphics::rect(
      10, 0.4, 16, 0.8,
      col = cor_ausencia,
      border = "#A1ACB5"
    )
    
    graphics::text(
      19, 0.6,
      "Ausência",
      adj = 0,
      cex = 0.92
    )
    
    graphics::rect(
      48, 0.4, 54, 0.8,
      col = cor_conexao,
      border = NA
    )
    
    graphics::text(
      57, 0.6,
      "Conexão",
      adj = 0,
      cex = 0.92
    )
  }
  
  # ----------------------------------------------------------
  # 7. FIGURA COMPARATIVA
  # ----------------------------------------------------------
  
  desenhar_comparacao <- function() {
    
    graphics::par(
      oma = c(0.2, 0.2, 2.2, 0.2),
      family = "sans",
      bg = "white"
    )
    
    graphics::layout(
      matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        byrow = TRUE
      ),
      heights = c(12, 1.6),
      widths = c(1, 1)
    )
    
    desenhar_matriz(
      Wbin,
      "Matriz anatômica original",
      original = TRUE
    )
    
    desenhar_matriz(
      freq,
      sprintf("Matriz recuperada", R),
      original = FALSE
    )
    
    desenhar_legenda_original()
    desenhar_legenda_freq()
    
    graphics::mtext(
      "Conectividade anatômica e frequência de seleção",
      side = 3,
      outer = TRUE,
      line = 0.6,
      cex = 1.25,
      font = 2
    )
    
    graphics::layout(1)
  }
  
  # ----------------------------------------------------------
  # 8. FIGURA RECUPERADA ISOLADA
  # ----------------------------------------------------------
  
  desenhar_recuperada <- function() {
    
    graphics::par(
      oma = c(0.2, 0.2, 0.5, 0.2),
      family = "sans",
      bg = "white"
    )
    
    graphics::layout(
      matrix(c(1, 2), nrow = 2),
      heights = c(12, 1.5)
    )
    
    desenhar_matriz(
      freq,
      sprintf(
        "Frequência de seleção por conexão (R = %d)",
        R
      ),
      original = FALSE
    )
    
    desenhar_legenda_freq()
    
    graphics::layout(1)
  }
  
  # ----------------------------------------------------------
  # 9. SALVAR COMPARACAO
  # ----------------------------------------------------------
  
  grDevices::png(
    file.path(
      DIR_OUT,
      "comparacao_original_recuperada_64x64.png"
    ),
    width = 3800,
    height = 2100,
    res = 260,
    type = "cairo"
  )
  
  desenhar_comparacao()
  grDevices::dev.off()
  
  grDevices::cairo_pdf(
    file.path(
      DIR_OUT,
      "comparacao_original_recuperada_64x64.pdf"
    ),
    width = 15.5,
    height = 8.8,
    onefile = FALSE
  )
  
  desenhar_comparacao()
  grDevices::dev.off()
  
  # ----------------------------------------------------------
  # 10. SALVAR RECUPERADA ISOLADA
  # ----------------------------------------------------------
  
  grDevices::png(
    file.path(
      DIR_OUT,
      "matriz_frequencia_64x64.png"
    ),
    width = 2800,
    height = 2900,
    res = 300,
    type = "cairo"
  )
  
  desenhar_recuperada()
  grDevices::dev.off()
  
  grDevices::cairo_pdf(
    file.path(
      DIR_OUT,
      "matriz_frequencia_64x64.pdf"
    ),
    width = 9.5,
    height = 10,
    onefile = FALSE
  )
  
  desenhar_recuperada()
  grDevices::dev.off()
  
  # ----------------------------------------------------------
  # 11. CONCLUSAO
  # ----------------------------------------------------------
  
  cat("\n====================================\n")
  cat("FIGURAS GERADAS COM SUCESSO!\n")
  cat("Replicas utilizadas:", R, "\n")
  cat("Diretorio:", DIR_OUT, "\n")
  cat("====================================\n")
  
})
