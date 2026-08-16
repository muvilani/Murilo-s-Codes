# CARREGANDO PACOTES ------------------------------------------------------
library(glmnet)
library(ggplot2)

# Gerando a Matriz de pesos sinápticos ----------------------------------------------
set.seed(15)
n = 4
w = matrix(c(     0,     4,   0, 0, #
                  0,     0,   6, 0, #
                  0,     0,   0, 0,
                  0,     0,   0, 0)     , nrow=n,ncol=n, byrow = T) #


# Funções necessárias para gerar X ----------------------------------------

#Função faz uma matriz onde cada linha representa o neurônio e cada coluna a proporção de 0 até aquele instante no tempo
matriz_prop_1 <- function(matriz,X,n,t){
  for (i in 1:n){
    for (j in 1:t+4){
      matriz[i,j] <- sum(X[i,][1:j])/j
    }
  }
  return(matriz)
}

#Função que seleciona todas as linhas da matriz "m" menos a do "indice" indicado 
matriz_pesos = function(m,indice){
  matriz_pesos = matrix(c(m[!(row(m) %in% c(indice))]), nrow=nrow(m)-1)
  return(matriz_pesos)
}

#Função que seleciona a linha "indice" da matriz "m" indicada
matriz_x = function(m,indice){
  matriz_x = c(m[(row(m) %in% c(indice))])
  return(matriz_x)
}

#Função taxa de disparo logísstica
phi = function(x){
  p = exp(x)/(1+exp(x))
  return(p)
}

#Função verifica indice do ultimo disparo que ocorreu na sequência de disparos x, que á a linha do neurônio especificada
ultimo = function(x,tempo_inicial,tempo_final){
  u = tail(which(x[tempo_inicial:tempo_final] == 1), n = 1)
  return(u)
}

#Função que calcula o potencial de membrana do neurônio "indice_neurônio" no tempo t 
V_novo = function(X,indice_neuronio,tempo){
  potencial = if(ultimo(matriz_x(X,indice_neuronio),1,tempo) == tempo) 0  
  else if(ultimo(matriz_x(X,indice_neuronio),1,tempo) != (tempo-1))    sum((1/(2^(tempo - ultimo(matriz_x(X,indice_neuronio),1,tempo))))*rowSums(matriz_pesos(X,indice_neuronio)[,(ultimo(matriz_x(X,indice_neuronio),1,tempo)+1):tempo])%*%matriz_pesos(w,indice_neuronio)[,indice_neuronio]) 
  else  sum((1/(2^(tempo - ultimo(matriz_x(X,indice_neuronio),1,tempo))))*matriz_pesos(X,indice_neuronio)[,(ultimo(matriz_x(X,indice_neuronio),1,tempo)+1):tempo]%*%matriz_pesos(w,indice_neuronio)[,indice_neuronio]) 
  return(potencial)
}


##Função que gera a cadeia X
#gera_cadeia <- function(X){
#  for (j in 5:ncol(X)){
#    for (i in 1:n){
#      X[i,j] = if (sum(w[1,]) > 0) ifelse(runif(1,0,1) <= phi(1 + V_novo(X,i,j-1)), 1, 0)
#      else ifelse(runif(1,0,1) <= phi(-1 + V_novo(X,i,j-1)), 1, 0)
#    }
#  }
#  return(X)
#}

##Função que gera a cadeia X sem influencia externa
gera_cadeia_2 <- function(X){
  for (j in 5:ncol(X)){
    for (i in 1:n){
      X[i,j] = ifelse(runif(1,0,1) <= phi(V_novo(X,i,j-1)), 1, 0)
    }
  }
  return(X)
}

t = 19996 #normalmente 7496

X = matrix(c(1,0,0,0,
             0,1,1,0,
             0,0,1,1,
             0,1,1,0), nrow=n, byrow=T)


X = cbind(X, matrix(0,nrow=n,ncol=t))

#X <- gera_cadeia(X)
X <- gera_cadeia_2(X)

#criando a matriz de proporções de 1 nas sequências e fazendo os gráficos para o burn in ---------------
proporcao_1 <- matrix(0, nrow =5, ncol =t+4) #gero a matriz so pra usar na função primeiro
proporcao_1 <- matriz_prop_1(proporcao_1,X,n,t)

# Criar um data frame com os dados para o gráfico
dados <- data.frame(Tempo = 1:ncol(X), 
                    Proporcao_1_Neuronio1 = proporcao_1[1,], 
                    Proporcao_1_Neuronio2 = proporcao_1[2,],
                    Proporcao_1_Neuronio3 = proporcao_1[3,],
                    Proporcao_1_Neuronio4 = proporcao_1[4,],
                    Proporcao_1_Neuronio5 = proporcao_1[5,])

# Transformar o data frame em formato longo
dados_longo <- tidyr::pivot_longer(dados, cols = -Tempo, names_to = "Neuronio", values_to = "Proporcao_1")

# Plotar o gráfico combinado
ggplot(dados_longo, aes(x = Tempo, y = Proporcao_1, color = Neuronio)) +
  geom_line() +
  labs(x = "Tempo", y = "Proporção de 1s") +
  scale_color_manual(values = c("blue", "red", "green", "purple", "orange")) +
  theme_minimal()

#aplicar o burn-in nos dados e salvar as 3 cadeias diferentes

#burn-in de 2500

#burn in da cadeia -----------------------------
X13 <- matrix(0, nrow = 5, ncol = t+4-2500)
for (i in 1:5){
  X13[i,] <- X[i,][2501:(t+4)]
}

##############
#X1 cadeia com pesos 4 entre os neurônios 3 e 4
#X2 cadeia com pesos -4 entre os neurônios 3 e 4
#X3 cadeia com pesos 4 e -4 respectivamente entre os neurônios 3 e 4
#X10 cadeia com 17500 observações pesos 4 e 4
#X11 cadeia com 17500 observações com neurônio 3 (peso 4) e neurônio 4 (peso -2 no neurônio 3)
#X12 cadeia com 2 peso -2 no 3, 3 e 4 peso 4 mutuamente entre eles.
#X12 cadeia com 2 peso -3 no 3, 3 e 4 peso 4 mutuamente entre eles.
save(X1,X2,X3,X4,X5,X6,X7,X1_pot,X2_pot,X3_pot, file = "cadeias")
load("cadeias")


