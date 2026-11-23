# CARREGANDO PACOTES ------------------------------------------------------
library(glmnet)
library(ggplot2)
library(infotheo)
library(philentropy)
library(fastmit)

# Function to atualize X matrix
matriz_x <- function(m, indice) {
  matriz_x <- m[indice, , drop = FALSE]
  return(matriz_x)
}

# logistic function
phi = function(x){
  p = exp(x)/(1+exp(x))
  return(p)
}

#Function to verify where the last 1 appeared in the time series sequence
ultimo = function(x,tempo_inicial,tempo_final){
  u = tail(which(x[tempo_inicial:tempo_final] == 1), n = 1)
  return(u)
}

# Function to atualize the W matrix
matriz_pesos <- function(m, indice) {
  matriz_pesos <- m[!(1:nrow(m) %in% indice), , drop = FALSE]
  return(matriz_pesos)
}

#Function to atualize the potential of the neuron in time t 
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

##Function to generate the chain of spikes from the neurons
gera_cadeia_2 <- function(X){
  for (j in 5:ncol(X)){
    for (i in 1:n){
      X[i,j] = ifelse(runif(1,0,1) <= phi(alpha[i] + V_novo(X,i,j-1)), 1, 0)
    }
  }
  return(X)
}

