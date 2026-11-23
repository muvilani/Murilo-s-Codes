# CARREGANDO PACOTES ------------------------------------------------------
library(glmnet)
library(ggplot2)
library(infotheo)
library(philentropy)
library(fastmit)
library(RTransferEntropy)

#Função faz uma matriz onde cada linha representa o neurônio e cada coluna a proporção de 0 até aquele instante no tempo
matriz_prop_1 <- function(matriz,X,n,t){
  for (i in 1:n){
    for (j in 1:t+4){
      matriz[i,j] <- sum(X[i,][1:j])/j
    }
  }
  return(matriz)
}

# Atualização da função matriz_x
matriz_x <- function(m, indice) {
  matriz_x <- m[indice, , drop = FALSE]
  return(matriz_x)
}

#Função taxa de disparo logísstica
phi = function(x){
  p = 1 / (1 + exp(-0.35*(x-0.5)))
  #p = exp(x)/(1+exp(x))
  return(p)
}

# Atualização da função matriz_pesos
matriz_pesos <- function(m, indice) {
  matriz_pesos <- m[!(1:nrow(m) %in% indice), , drop = FALSE]
  return(matriz_pesos)
}

V_novo_fixa <- function(X, indice_neuronio, tempo) {
  x <- matriz_x(X, indice_neuronio)
  k <- 4  # Número fixo de passos anteriores
  inicio <- tempo - k + 1
  adjusted_w <- matriz_pesos(w, indice_neuronio)[, rep(indice_neuronio, k)]
  potencial <- sum((1 / (2^k)) * matriz_pesos(X, indice_neuronio)[,(inicio:tempo)] %*% adjusted_w)
  return(potencial)
}

##Função que gera a cadeia X sem influencia externa
gera_cadeia_2_fixa <- function(X){
  for (j in 5:ncol(X)){
    for (i in 1:n){
      X[i,j] = ifelse(runif(1,0,1) <= phi(V_novo_fixa(X,i,j-1)), 1, 0)
    }
  }
  return(X)
}

#executar novamente a função gera_cadeia_2 com apenas 2 neurônios
n <- 2
set.seed(1)
t = 20
w <- matrix(c(0, 2, 
              2, 0), nrow = n, ncol = n, byrow = TRUE)

X <- matrix(c(1, 0, 0, 1,
              0, 1, 1, 0), nrow = n, byrow = TRUE)

X <- cbind(X, matrix(0, nrow = n, ncol = t))

X1 <- gera_cadeia_2_fixa(X)

#criando a matriz de proporções de 1 nas sequências e fazendo os gráficos para o burn in ---------------
proporcao_1 <- matrix(0, nrow = n, ncol =t+4) #gero a matriz so pra usar na função primeiro
proporcao_1 <- matriz_prop_1(proporcao_1,X,n,t)

# Criar um data frame com os dados para o gráfico
dados <- data.frame(Tempo = 1:ncol(X), 
                    Proporcao_1_Neuronio1 = proporcao_1[1,], 
                    Proporcao_1_Neuronio2 = proporcao_1[2,]
)

# Transformar o data frame em formato longo
dados_longo <- tidyr::pivot_longer(dados, cols = -Tempo, names_to = "Neuronio", values_to = "Proporcao_1")

# Plotar o gráfico combinado
ggplot(dados_longo, aes(x = Tempo, y = Proporcao_1, color = Neuronio)) +
  geom_line() +
  labs(x = "Tempo", y = "Proporção de 1s") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

#burn in da cadeia -----------------------------
X1 <- matrix(0, nrow = n, ncol = t+4-2500)
for (i in 1:n){
  X1[i,] <- X[i,][2501:(t+4)]
}

entropia_transferencia(X1[1,],X1[2,],3)
transfer_entropy(X1[2,],X1[1,], lx = 4, ly = 4) 
calc_te(X1[2,],X1[1,], lx = 4, ly = 4)
0.0068*42500*2
