# CARREGANDO PACOTES ------------------------------------------------------
library(glmnet)
library(ggplot2)
library(infotheo)
library(philentropy)
library(fastmit)

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
  #p = 1 / (1 + exp(-0.35*(x-0.5)))
  #p = exp(3*x) / (1 + exp(3*x))
  p = exp(x)/(1+exp(x))
  return(p)
}

#Função verifica indice do ultimo disparo que ocorreu na sequência de disparos x, que á a linha do neurônio especificada
ultimo = function(x,tempo_inicial,tempo_final){
  u = tail(which(x[tempo_inicial:tempo_final] == 1), n = 1)
  return(u)
}

# Atualização da função matriz_pesos
matriz_pesos <- function(m, indice) {
  matriz_pesos <- m[!(1:nrow(m) %in% indice), , drop = FALSE]
  return(matriz_pesos)
}

V_novo <- function(X, indice_neuronio, tempo) {
  x <- matriz_x(X, indice_neuronio)
  indice <- ultimo(x, 1, tempo)
  dim_matriz_pesos_X <- dim(matriz_pesos(X, indice_neuronio))
  dim_matriz_pesos_w <- dim(matriz_pesos(w, indice_neuronio))
  potencial <- if (indice == tempo) {
    0
  } else if (indice != (tempo - 1)) {
    adjusted_w <- matriz_pesos(w, indice_neuronio)[, rep(indice_neuronio, tempo - indice)]
    sum((1 / (2^(tempo - indice))) * matriz_pesos(X, indice_neuronio)[,(indice + 1):tempo] %*% adjusted_w)
  } else {
    adjusted_w <- matriz_pesos(w, indice_neuronio)[, rep(indice_neuronio, tempo - indice + 1)]
    sum((1 / (2^(tempo - indice))) * matriz_pesos(X, indice_neuronio)[,(indice + 1):tempo] %*% adjusted_w)
  }
  return(potencial)
}

##Função que gera a cadeia X sem influencia externa
gera_cadeia_2 <- function(X){
  for (j in 5:ncol(X)){
    for (i in 1:n){
      X[i,j] = ifelse(runif(1,0,1) <= phi(alpha[i] + V_novo(X,i,j-1)), 1, 0)
    }
  }
  return(X)
}

#executar novamente a função gera_cadeia_2 com apenas 2 neurônios
n <- 2
set.seed(1)
t = 4996
w <- matrix(c(0, 10, 
              0, 0), nrow = n, ncol = n, byrow = TRUE)
alpha <- c(0.0,0.0)

X <- matrix(c(1, 0, 0, 1,
              0, 1, 1, 0), nrow = n, byrow = TRUE)

X <- cbind(X, matrix(0, nrow = n, ncol = t))

X2 <- gera_cadeia_2(X)
table(X2[1,])
table(X2[2,])
16067/40000
31700/40000
#criando a matriz de proporções de 1 nas sequências e fazendo os gráficos para o burn in ---------------
proporcao_1 <- matrix(0, nrow = n, ncol =t+4) #gero a matriz so pra usar na função primeiro
proporcao_1 <- matriz_prop_1(proporcao_1,X2,n,t)
#proporcao_1[1,]
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

#aplicar o burn-in nos dados e salvar as 3 cadeias diferentes

#burn-in de 2500

#burn in da cadeia -----------------------------
X1 <- matrix(0, nrow = n, ncol = t+4-2500)
for (i in 1:n){
  X1[i,] <- X[i,][2501:(t+4)]
}

###############################################################################
####################### CHAMA A FUNÇÃO PARA B AMOSTRAS ##################################
B = 100
valores_obs <- rep(1,B)
valores_obs2 <- rep(1,B)
valores_obs3 <- rep(1,B)
set.seed(1)
n=2
w = matrix(c(     0,     0,
                  0,     0), nrow=n,ncol=n, byrow = T)
t= 4996
system.time(print(estima_qschiq(B,t,valores_obs,valores_obs2)))

qchisq(0.05,15,lower.tail = TRUE)
qchisq(0.05,1,lower.tail = TRUE)

#####################################################################################
############################## TESTES PARA FUNÇÃO PHI  ##############################
# Definindo a função p(x) mudar b para deixar ela mais inclinada sem tirar a taxa de 0 pra 0.5
p <- function(x) {
  exp('b'x) / (1 + exp('b'x))
}

'p <- function(x) {
  1 / (1 + exp(-0.35*(x-0.5)))
}'


# Plot da função
curve(p, from = -10, to = 10, n = 1000, col = "blue", xlab = "x", ylab = "p(x)",
      main = "Gráfico da função p(x)", ylim = c(0, 1))
abline(v = 0, lty = 2)
abline(h = 0.5, lty = 2)
phi(0)
