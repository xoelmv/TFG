
#----------------------------#
#---  CARGAMOS LIBRERIAS  ---#
#----------------------------#

library(fclust)
library(TSclust)
library(stats)




#----------------------------#
#---  AJUSTAR PARAMETROS  ---#
#----------------------------#

# valores de m para fuzzy
m_values <- c(seq(from =1.2, to=1.8, by=0.3))
#m_values <- c(1.2, 1.8)
#----------------------------#
#---  AUXILIAR FUNCTIONS  ---#
#----------------------------#



# Función que calcula las autocorrelaciones:
Acf=function(x, lag.max = 10){
  
  if (dim(x)[2]>10){
    coef_ACF_series = matrix(0,dim(x)[1] , ncol = lag.max)
  } else {
    coef_ACF_series = matrix(0,dim(x)[1] , ncol = dim(x)[2]-1)
  }
  
  for (i in 1:dim(x)[1]){
    coeficientes <- stats::acf(x[i,], lag.max = lag.max, plot = FALSE)$acf[-1]
    coef_ACF_series[i,] <- coeficientes
  }
  
  row.names(coef_ACF_series) <- row.names(x)
  return(coef_ACF_series)
}

# Función que calcula las autocorrelaciones parciales:
Pacf=function(x, lag.max = 10){
  
  if (dim(x)[2]>10){
    coef_PACF_series = matrix(0,dim(x)[1] , ncol = lag.max)
  } else {
    coef_PACF_series = matrix(0,dim(x)[1] , ncol = dim(x)[2]-1)
  }
  
  
  for (i in 1:dim(x)[1]){
    coeficientes <- stats::pacf(x[i,], lag.max = lag.max, plot = FALSE)$acf
    coef_PACF_series[i,] <- coeficientes
  }
  
  row.names(coef_PACF_series) <- row.names(x)
  return(coef_PACF_series)
}


# Función que calcula los coeficientes autoregresivos (Piccolo): 
pic=function(x){
  
  autoregressive_coef = list()
  max.order = 0
  
  for (i in 1:dim(x)[1]){
    output <- stats::ar(x[i,])
    ar <- output$ar
    autoregressive_coef[[i]] <- ar
    
    order <- output$order
    if (order > max.order){
      max.order <- order
    }
    if (order == 0){
      autoregressive_coef[[i]] <- arima(x[i,], order = c(1, 0, 0))$coef[1]
    }
    
  }
  
  X_pic = matrix(0,dim(x)[1] , ncol = max.order)
  for (i in 1:dim(x)[1]){
    
    ar <- autoregressive_coef[[i]]
    if (length(ar) > 0){
      X_pic[i,] <- c(ar, rep(0, max.order - length(autoregressive_coef[[i]]) ) )
    }
  }
  
  row.names(X_pic) <- row.names(x)
  return(X_pic)
}


# Función que calcula la matriz de autocovarianzas cuantil:
qaf=function(x,tau,L){
  
  x=x[!is.na(x)]
  n=length(tau)
  
  # Quantiles:
  quant=quantile(x,probs=tau,na.rm=TRUE)
  
  xx=matrix(0,ncol=length(quant),nrow=length(x))
  for (i in 1:length(quant)) xx[,i]=((x<=quant[i])*1)-tau[i]
  N=dim(xx)[1]
  
  comb=expand.grid(1:n,1:n)
  I=dim(comb)[1]
  
  list.cov.mat=list()
  auto.cov=NA; name=NA
  for (l in 1:L){
    cov.mat=matrix(0,ncol=n,nrow=n)
    for(i in 1:I){
      cov.mat[comb[i,1],comb[i,2]]=cov(xx[-((N-l+1):N),comb[i,1]],xx[-(1:l),comb[i,2]])
      name=c(name,paste(l,"-",comb[i,2],"-",comb[i,1],sep=""))
    }
    colnames(cov.mat)=names(quant); rownames(cov.mat)=names(quant)
    list.cov.mat[[l]]=cov.mat
    auto.cov=c(auto.cov,as.vector(t(cov.mat)))
    
  }
  auto.cov=auto.cov[-1]
  names(auto.cov)=name[-1]
  names(list.cov.mat)=paste("lag",1:L,sep="")
  return(auto.cov)
}
# Función que calcula los medoides iniciales para arrancar el algoritmo con el método PAM:
init.part = function(X,K){
  diss1 <- dist(X)^2 
  clust1 <- cutree(hclust(diss1, method = "complete"), k = K)
  names(clust1) <- rownames(X)
  diss1 <- as.matrix(diss1)
  colnames(diss1) = rownames(diss1) = rownames(X)
  
  # Compute medoids:
  medoids <- numeric()
  for(i in 1:K){
    if(sum(clust1 == i) == 1){
      medoids[i] <- names(which(clust1 == i))
    }else{
      medoids[i] <- names(which.min(colSums(diss1[clust1 == i, clust1 == i])))
    }
  }
  
  C.index <- sort(match(medoids, rownames(X)))
  return(C.index)
}


# Función que calcula los "membership degrees":
Uic.FCMdC = function(dic, m){
  
  dic_m <- dic^(1/(m-1))
  
  U <- 1/(dic_m*rowSums(1/dic_m))
  
  U[U == "NaN"] <- 1
  
  return(U)
  
}


#-------------------#
#---  QAF-FCMdC  ---#
#-------------------#



# FUZZY C-MEDOIDS: Función que ejecuta todo el algoritmo
# Y: Matriz que contiene el conjunto de series a clasificar (por fila)
# m: parámetro que controla lo 'fuzzy' que será la clasificación
# C: número de clusters
# max.iter: número máximo de iteraciones
# aleat: si su valor es TRUE, la partición inicial se selecciona de forma aleatoria, en caso contrario se computa el método PAM
FCMdC_QAF=function(Y, D, X, m, C, max.iter, aleat){
  
  # Transforma Y en un data.frame y se renombran las filas
  Y <- data.frame(Y)
  rownames(Y) <- paste("x",1:dim(Y)[1],sep="")
  
  # Almacena en n el número de series a clasificar
  n <- dim(Y)[1]
  
  X <- X                 # features matrix
  distance <- D          # dissimilarity matrix 
  
  # Se decide entre elegir la partición de forma aleatoria (si aleat = TRUE)
  # o calcular la partición inicial mediante algoritmo PAM basado en la métrica QAF (si aleat = FALSE)
  # En este punto, en medoids se almacena el número de fila de cada uno de los medoides iniciales.
  if(aleat){
    medoids <- sort(sample(1:n, size = C))
  }else{
    medoids <- init.part(X, C)
  }
  
  # Se inicializa el valor de iteración a 1
  iter <- 1
  
  # Se define la matriz en la que se iran almacenando los valores de la función objetivo para cada medoide
  fobj_c <- matrix(0, nrow = max.iter, ncol = C)
  
  # Bucle que computa el algoritmo propiamente dicho hasta que para, ya sea porque converge
  # o porque se alcance el número máximo de iteraciones
  for(j in 1:max.iter){
    
    # A partir de la segunda iteración, se almacena el valor de los medoides de dos iteraciones antes
    # Se explica más adelante su uso en un criterio de convergencia del algoritmo
    if(iter > 1) medoids_old2 <- medoids_old
    
    # Se almacena el valor de los medoides de la iteración anterior
    # Se explica más adelante su uso en un criterio de convergencia del algoritmo
    medoids_old <- medoids
    
    # Distancias de cada serie a los medoides 
    d_ic <- distance[, medoids, drop = FALSE]
    
    # Cálculo de los nuevos "membership degrees"
    U <- Uic.FCMdC(d_ic, m)
    
    # Definición de la matriz de distancias que se actualiza en cada iteración
    fobj <- matrix(Inf, ncol = C, nrow = n)
    
    # Bucle que se encarga de actualizar el valor de los medoides
    for(c in 1:C){
      
      # Se calcula el valor de la fuinción objetivo para cada una de las series que no son medoides
      for (i in 1:n) fobj[i, c] <- U[-medoids, c]^m %*% distance[-medoids, i]
      
      # Se le asocia valor infinito a los medoides para evitar que pueda seleccionarse la misma serie como dos medoides
      fobj[medoids, c] <- Inf
      
      # Calcula la función objetivo para el medoide
      fobj_c[j, c] <- U[-medoids, c]^m %*% d_ic[-medoids,c]
      
      # En caso de que el valor de la función objetivo de alguna de las series sea menor que el de el c-esimo medoide
      # Se actualiza el nuevo medoide
      if(min(fobj[, c]) < fobj_c[j, c]){
        
        # Se identifica la serie que hace mínima la función objetivo y se actualiza el valor del c-ésimo medoide
        medoids[c] <- which.min(fobj[, c])
        
      }
    }
    
    # Comprueba si en dos iteraciones consecutivas los medoides no cambian
    # En caso afirmativo, el algoritmo para
    if(all(sort(medoids) == sort(medoids_old)) | iter == max.iter) break
    
    # En caso contrario, se realiza una nueva iteración hasta llegar al número máximo de iteraciones
    # Puede suceder que se entre en bucle entre dos particiones distintas, en ese caso se define el siguiente criterio de parada.
    # En caso de qe en dos iteraciones no consecuitivas se repita la misma partición, se seleccionara de las dos que entran en conflicto,
    # aquella que minimice la función objetivo
    if((iter>1) && (all(sort(medoids) == sort(medoids_old2)))){     
      
      if(min(fobj_c[iter, ]) < min(fobj_c[iter-1, ])){ 
        
        medoids <- sort(medoids) 
        
      }else{
        
        medoids <- sort(medoids_old) 
        
      }
      
      break
      
    }
    
    # Se actualiza el valor de la iteración
    iter <- iter+1
    
  }
  
  # Una cez que el bucle termina, ya sea por alcanzar el número máximo de iteraciones o porque el algoritmo converja
  # se calcula la matriz de "membership degrees", los medoides y la partición final.
  
  # Cálculo de la matriz de "membership degrees":
  medoids_opt <- medoids
  d_ic_opt <- distance[, medoids_opt, drop = FALSE]
  U_opt <- Uic.FCMdC(d_ic_opt, m)                 
  
  
  # Cálculo de la partición "hard" final seleccionando el cluster con el mayor "memberhip degree" para cada serie
  clust <- apply(U_opt, 1, which.max)
  names(clust) <- rownames(U_opt)
  
  # Valores finales devueltos por la función
  results <- list()
  results$U <- U_opt
  results$clustering <- clust
  results$medoids <- medoids_opt
  results$n_iter <- iter
  results$k <- C
  results$m <- m
  results$Y <- Y
  results$d_ic <- d_ic_opt
  results$distance <- distance
  return(results)
}



#--------------------------------#
#---  CREAMOS LOS ESCENARIOS  ---#
#--------------------------------#


f_escenario.1 = function(R=10, S=10, nclust=2, T=200, ar.min_c1, ar.max_c1,ar.min_c2, ar.max_c2){
  results = list()
  R_series=list()
  for (r in 1:R){
    series <- matrix(0, nrow = S * nclust, ncol = T)
    serie_names <- character()
    
    for (i in 1:(S * nclust)) {
      if (i <= S) {
        ar_params <- runif(1, min = ar.min_c1, max = ar.max_c1)
        serie_names <- c(serie_names, paste0("serie", i, ".1"))
      } else {
        ar_params <- runif(1, min = ar.min_c2, max = ar.max_c2)
        serie_names <- c(serie_names, paste0("serie", i, ".2"))
      }
      serie <- arima.sim(list(order = c(1, 0, 0), ar = ar_params), n = T)
      #output <- stats::ar(serie)
      #print(paste("Serie:", i, ", AR:", ar_params, ", AR estimado:", output$ar))
      series[i, ] <- serie
      
    }
    rownames(series) <- serie_names
    R_series[[r]] <- series
  }
  true_labels <- c(rep(1, S), rep(2, S))
  results$datos <- R_series
  results$R <- R   
  results$C <- nclust
  results$true_labels <- true_labels
  return(results)
}

f_escenario.2 = function(R=10, S=10, nclust=2, T=200, ma.min_c1, ma.max_c1,ma.min_c2, ma.max_c2){
  results = list()
  R_series=list()
  for (r in 1:R){
    series <- matrix(0, nrow = S * nclust, ncol = T)
    serie_names <- character()
    
    for (i in 1:(S * nclust)) {
      if (i <= S) {
        ma_params <- runif(1, min = ma.min_c1, max = ma.max_c1)
        serie_names <- c(serie_names, paste0("serie", i, ".1"))
      } else {
        ma_params <- runif(1, min = ma.min_c2, max = ma.max_c2)
        serie_names <- c(serie_names, paste0("serie", i, ".2"))
      }
      serie <- arima.sim(list(order = c(0, 0, 1), ma = ma_params), n = T)
      #output <- stats::ar(serie)
      #print(paste("Serie:", i, ", AR:", ar_params, ", AR estimado:", output$ar))
      series[i, ] <- serie
      
    }
    rownames(series) <- serie_names
    R_series[[r]] <- series
  }
  true_labels <- c(rep(1, S), rep(2, S))
  results$datos <- R_series
  results$R <- R   
  results$C <- nclust
  results$true_labels <- true_labels
  return(results)
}

f_escenario.3 = function(R=10, S=10, nclust=3, T=200, ar.min_c1, ar.max_c1,ar.min_c2, ar.max_c2, ar.min_c3, ar.max_c3){
  results = list()
  R_series=list()
  for (r in 1:R){
    series <- matrix(0, nrow = S * nclust, ncol = T)
    serie_names <- character()
    
    for (i in 1:(S * nclust)) {
      if (i <= S) {
        ar_params <- runif(1, min = ar.min_c1, max = ar.max_c1)
        serie_names <- c(serie_names, paste0("serie", i, ".1"))
      } else if (i > S & i <= 2*S) {
        ar_params <- runif(1, min = ar.min_c2, max = ar.max_c2)
        serie_names <- c(serie_names, paste0("serie", i, ".2"))
      } else {
        ar_params <- runif(1, min = ar.min_c3, max = ar.max_c3)
        serie_names <- c(serie_names, paste0("serie", i, ".3"))
      }
      serie <- arima.sim(list(order = c(1, 0, 0), ar = ar_params), n = T)
      #output <- stats::ar(serie)
      #print(paste("Serie:", i, ", AR:", ar_params, ", AR estimado:", output$ar))
      series[i, ] <- serie
      
    }
    rownames(series) <- serie_names
    R_series[[r]] <- series
  }
  true_labels <- c(rep(1, S), rep(2, S), rep(3, S))
  results$datos <- R_series
  results$R <- R   
  results$C <- nclust
  results$true_labels <- true_labels
  return(results)
}


# sea phi.1 el coeficiente autorregresivo de las series del grupo 1, y phi.2 el correspondiente
# a las series del grupo 2, se crea una serie cuya primera mitad (T/2) se corresponde a un AR(phi.1)
# y el siguente tramo a un AR(phi.2) 

serie_2phi = function(phi.1,phi.2,T){
  burn = 100
  a <- rnorm(T+burn)
  x.aux <- numeric(T+burn)
  
  
  x.aux[1] <- a[1]
  for (i in 2:(burn+T/2)){
    x.aux[i]<- x.aux[i-1]*phi.1 + a[i]
  }
  for (j in 1:(T/2)){
    x.aux[i+j]<- x.aux[i+j-1]*phi.2 + a[i+j]
  }
  x <- x.aux[101:(T+burn)]
  
  return(x)
}

f_escenario.sin_uso = function(R=10, S=10, nclust=2, T=200, ar.min_c1, ar.max_c1,ar.min_c2, ar.max_c2){
  results = list()
  R_series=list()
  for (r in 1:R){
    series <- matrix(0, nrow = S * nclust+1, ncol = T)
    serie_names <- character()
    
    for (i in 1:(S * nclust)) {
      if (i <= S) {
        ar_params <- runif(1, min = ar.min_c1, max = ar.max_c1)
        ar_params.1 <- ar_params
        serie_names <- c(serie_names, paste0("serie", i, ".1"))
      } else {
        ar_params <- runif(1, min = ar.min_c2, max = ar.max_c2)
        ar_params.2 <- ar_params
        serie_names <- c(serie_names, paste0("serie", i, ".2"))
      }
      serie <- arima.sim(list(order = c(1, 0, 0), ar = ar_params), n = T)
      #output <- stats::ar(serie)
      #print(paste("Serie:", i, ", AR:", ar_params, ", AR estimado:", output$ar))
      series[i, ] <- serie
      
    }
    serie <- serie_2phi(ar_params.1, ar_params.2, T=200)
    #cat(ar_params.1, ar_params.2)
    serie_names <- c(serie_names, paste0("serie", i+1, ".2phi(",ar_params.1, ",", ar_params.2, ")"))
    series[i+1, ]<- serie
    rownames(series) <- serie_names
    R_series[[r]] <- series
  }
  true_labels <- c(rep(1, S), rep(2, S), 1)
  results$datos <- R_series
  results$R <- R   
  results$C <- nclust
  results$true_labels <- true_labels
  return(results)
}
# se crea una serie cuyo coeficiente autoregresivo (phi.3) está entre los dos grupos.
# phi.1 el coeficiente autorregresivo de las series del grupo 1, min = 0, max = 0.4
# phi.2 el coeficiente autorregresivo de las series del grupo 2, min = 0.6, max = 1)
# phi.3 = 0.5

f_escenario.4 = function(R=10, S=10, nclust=2, T=200, ar.min_c1, ar.max_c1,ar.min_c2, ar.max_c2, ar.extra){
  results = list()
  R_series=list()
  for (r in 1:R){
    series <- matrix(0, nrow = S * nclust+1, ncol = T)
    serie_names <- character()
    
    for (i in 1:(S * nclust)) {
      if (i <= S) {
        ar_params <- runif(1, min = ar.min_c1, max = ar.max_c1)
        ar_params.1 <- ar_params
        serie_names <- c(serie_names, paste0("serie", i, ".1"))
      } else {
        ar_params <- runif(1, min = ar.min_c2, max = ar.max_c2)
        ar_params.2 <- ar_params
        serie_names <- c(serie_names, paste0("serie", i, ".2"))
      }
      serie <- arima.sim(list(order = c(1, 0, 0), ar = ar_params), n = T)
      #output <- stats::ar(serie)
      #print(paste("Serie:", i, ", AR:", ar_params, ", AR estimado:", output$ar))
      series[i, ] <- serie
      
    }
  
    serie <- arima.sim(list(order = c(1, 0, 0), ar = ar.extra), n = T)
    serie_names <- c(serie_names, paste0("serie", i+1, ".phi_", ar.extra))
    series[i+1, ]<- serie
    rownames(series) <- serie_names
    R_series[[r]] <- series
  }
  true_labels <- c(rep(1, S), rep(2, S), 1)
  results$datos <- R_series
  results$R <- R
  results$C <- nclust
  results$true_labels <- true_labels
  return(results)
}

# Escenario con series no lineales
f_escenario.5 <- function(R=10, S=4, nclust=3, T=100,
                          C1.p1.min, C1.p1.max, C1.p2.min, C1.p2.max,
                          C2.p1.min, C2.p1.max, C2.p2.min, C2.p2.max,
                          C3.p1.min, C3.p1.max, C3.p2.min, C3.p2.max) {
  results <- list()
  R_series <- list()
  
  # Función para simular serie bilineal
  simulate_bilinear <- function(T, p1.min, p1.max, p2.min, p2.max) {
    T2 <- T + 50
    p1 <- runif(1,p1.min, p1.max)
    p2 <- runif(1,p2.min, p2.max)
    X <- numeric(T2)
    epsilon <- rnorm(T2)
    X[1] <- epsilon[1]
    for (t in 2:T2) {
      X[t] <- p1 * X[t-1] - p2 * epsilon[t-1] * X[t-1] + epsilon[t]
    }
    return(X[51:T2])
  }
  
  
  # Bucle para cada repetición
  for (r in 1:R) {
    series <- matrix(0, nrow = S * nclust, ncol = T)
    serie_names <- character()
    
    # Bucle para cada serie
    for (i in 1:(S * nclust)) {
      if (i <= S) {
        serie_names <- c(serie_names, paste0("serie", i, ".1"))
        series[i, ] <- simulate_bilinear(T, C1.p1.min, C1.p1.max, C1.p2.min, C1.p2.max)
      } else if (i > S & i <= 2*S) {
        serie_names <- c(serie_names, paste0("serie", i, ".2"))
        series[i, ] <- simulate_bilinear(T, C2.p1.min, C2.p1.max, C2.p2.min, C2.p2.max)
      } else {
        serie_names <- c(serie_names, paste0("serie", i, ".3"))
        series[i, ] <- simulate_bilinear(T, C3.p1.min, C3.p1.max, C3.p2.min, C3.p2.max)
      }
    }
    
    # Asigna nombres a las filas y almacena la serie simulada
    rownames(series) <- serie_names
    R_series[[r]] <- series
  }
  
  true_labels <- c(rep(1, S), rep(2, S), rep(3, S))
  results$datos <- R_series
  results$R <- R
  results$C <- nclust
  results$true_labels <- true_labels
  
  return(results)
}

## con escenario.5 base más una serie de tiempo outlier
f_escenario.6 <- function(R=10, S=4, nclust=2, T=100, C1.p1.min, C1.p1.max, C1.p2.min, C1.p2.max,
                          C2.p1.min, C2.p1.max, C2.p2.min, C2.p2.max,
                          Extra.p1.min, Extra.p1.max, Extra.p2.min, Extra.p2.max) {
  results <- list()
  R_series <- list()
  
  # Función para simular serie bilineal
  simulate_bilinear <- function(T, p1.min, p1.max, p2.min, p2.max) {
    T2 <- T + 50
    p1 <- runif(1,p1.min, p1.max)
    p2 <- runif(1,p2.min, p2.max)
    X <- numeric(T2)
    epsilon <- rnorm(T2)
    X[1] <- epsilon[1]
    for (t in 2:T2) {
      X[t] <- p1 * X[t-1] - p2 * epsilon[t-1] * X[t-1] + epsilon[t]
    }
    return(X[51:T2])
  }
  
  # Bucle para cada repetición
  for (r in 1:R) {
    series <- matrix(0, nrow = (S * nclust) + 1, ncol = T)  # Añadimos una fila más para la serie outlier
    serie_names <- character()
    
    # Bucle para cada serie
    for (i in 1:((S * nclust) + 1)) {
      if (i <= S) {
        serie_names <- c(serie_names, paste0("serie", i, ".1"))
        series[i, ] <- simulate_bilinear(T, C1.p1.min, C1.p1.max, C1.p2.min, C1.p2.max)
      } else if (i <= (S * nclust)) {
        serie_names <- c(serie_names, paste0("serie", i, ".2"))
        series[i, ] <- simulate_bilinear(T, C2.p1.min, C2.p1.max, C2.p2.min, C2.p2.max)
      } else {
        serie_names <- c(serie_names, "outlier.1")
        series[i, ] <- simulate_bilinear(T, Extra.p1.min, Extra.p1.max, Extra.p2.min, Extra.p2.max)
      }
    }
    
    # Asigna nombres a las filas y almacena la serie simulada
    rownames(series) <- serie_names
    R_series[[r]] <- series
  }
  
  true_labels <- c(rep(1, S), rep(2, S), 1)
  results$datos <- R_series
  results$R <- R
  results$C <- nclust 
  results$true_labels <- true_labels
  
  return(results)
}


## tres clusters con series no lineales
f_escenario.7 <- function(R=10, S=4, nclust=3, T=100, C1.p1.min, C1.p1.max, 
                          C2.p1.min, C2.p1.max, C2.p2.min, C2.p2.max,
                          C3.p1.min, C3.p1.max) {
  results <- list()
  R_series <- list()
  
  # Función para simular serie exponencial autoregresiva (cluster C1)
  simulate_exp_ar <- function(T, p1.min, p1.max) {
    T2 <- T + 50
    p1 <- runif(1,p1.min, p1.max)
    X <- numeric(T2)
    epsilon <- rnorm(T2)
    for (t in 2:T2) {
      X[t] <- (p1 - 10 * exp(-X[t-1]^2)) * X[t-1] + epsilon[t]
    }
    return(X[51:T2])
  }
  
  # Función para simular serie bilineal (cluster C2)
  simulate_bilinear <- function(T, p1.min, p1.max, p2.min, p2.max) {
    T2 <- T + 50
    p1 <- runif(1,p1.min, p1.max)
    p2 <- runif(1,p2.min, p2.max)
    X <- numeric(T2)
    epsilon <- rnorm(T2)
    X[1] <- epsilon[1]
    for (t in 2:T2) {
      X[t] <- p1 * X[t-1] - p2 * epsilon[t-1] * X[t-1] + epsilon[t]
    }
    return(X[51:T2])
  }
  
  
  # Función para simular las series del tercer cluster (C3)
  simulate_outlier <- function(T, p1.min, p1.max) {
    T2 <- T + 50
    p1 <- runif(1,p1.min, p1.max)
    X <- numeric(T2)
    epsilon <- rnorm(T2)
    for (t in 2:T2) {
      X[t] <- p1 * abs(X[t-1]) * (3 + abs(X[t-1]))^-1 + epsilon[t]
    }
    return(X[51:T2])
  }
  
  # Bucle para cada repetición
  for (r in 1:R) {
    series <- matrix(0, nrow = (S * nclust), ncol = T)  # Añadimos una fila más para la serie outlier
    serie_names <- character()
    
    # Bucle para cada serie
    for (i in 1:(S * nclust)) {
      if (i <= S) {
        serie_names <- c(serie_names, paste0("serie", i, ".1"))
        series[i, ] <- simulate_exp_ar(T, C1.p1.min, C1.p1.max)
      } else if (i > S & i <= 2*S) {
        serie_names <- c(serie_names, paste0("serie", i, ".2"))
        series[i, ] <- simulate_bilinear(T, C2.p1.min, C2.p1.max, C2.p2.min, C2.p2.max)
      } else {
        serie_names <- c(serie_names, paste0("serie", i, ".3"))
        series[i, ] <- simulate_outlier(T, C3.p1.min, C3.p1.max)
      }
    }
    
    # Asigna nombres a las filas y almacena la serie simulada
    rownames(series) <- serie_names
    R_series[[r]] <- series
  }
  
  true_labels <- c(rep(1, S), rep(2, S), rep(3, S))  # Etiquetas verdaderas, 3 para el outlier
  results$datos <- R_series
  results$R <- R
  results$C <- nclust 
  results$true_labels <- true_labels
  
  return(results)
}


#----------------------------------------#
#---  Simulamos las series de tiempo  ---#
#----------------------------------------#

escenarios <- list()

set.seed(123)

# Dos cluster de series AR. 
# coeficiente autorregresivo C1 moviendose entre 0.1 y 0.4
# coeficiente autorregresivo C2 moviendose entre 0.6 y 0.9
R_series1 <- f_escenario.1(R=1, S=10, nclust=2, T=50,
                           ar.min_c1=0.1, ar.max_c1=0.4,
                           ar.min_c2=0.6, ar.max_c2=0.9)

escenarios[[1]] <- list()
escenarios[[1]]$S <- 10
escenarios[[1]]$nclust <- 2


# Dos cluster de series MA. 
# coeficiente de medias moviles C1 moviendose entre 0.3 y 0.8
# coeficiente de medias moviles C2 moviendose entre -0.8 y -0.3
#R_series2 <- f_escenario.2(R=1, S=10, nclust=2, T=200, ma.min_c1=0.3, ma.max_c1=0.8,ma.min_c2=-0.8, ma.max_c2=-0.3)




# Tres cluster de series AR. 
# coeficiente autorregresivo C1 moviendose entre 0.4 y 0.8
# coeficiente autorregresivo C2 moviendose entre -0.8 y -0.4
# coeficiente autorregresivo C3 moviendose entre -0.2 y 0.2
R_series2 <- f_escenario.3(R=1, S=10, nclust=3, T=50,
                           ar.min_c1=0.4, ar.max_c1=0.8,
                           ar.min_c2=-0.8, ar.max_c2=-0.4,
                           ar.min_c3=-0.2, ar.max_c3=0.2)

escenarios[[2]] <- list()
escenarios[[2]]$S <- 10
escenarios[[2]]$nclust <- 3

# Mismo escenario que el anterior pero con series de longitud T = 100
#R_series3 <- f_escenario.3(R=1, S=10, nclust=3, T=100,
#                           ar.min_c1=0.4, ar.max_c1=0.8,
#                           ar.min_c2=-0.8, ar.max_c2=-0.4,
#                           ar.min_c3=-0.2, ar.max_c3=0.2)

R_series3 <- f_escenario.4(R=1, S=10, nclust=2, T=200,
              ar.min_c1=0.1, ar.max_c1=0.4,
              ar.min_c2=0.6, ar.max_c2=0.9,
              ar.extra = 0.5)

escenarios[[3]] <- list()
escenarios[[3]]$S <- 10
escenarios[[3]]$nclust <- 2 #3



# Dos cluster de series AR. Y una serie extra a caballo de los dos clusters
# coeficiente autorregresivo C1 moviendose entre 0.1 y 0.4
# coeficiente autorregresivo C2 moviendose entre 0.6 y 0.9
# coeficiente autorregresivo C3 igual a 0.5
R_series4 <- f_escenario.4(R=1, S=10, nclust=2, T=50,
                           ar.min_c1=0.1, ar.max_c1=0.4,
                           ar.min_c2=0.6, ar.max_c2=0.9,
                           ar.extra = 0.5)

escenarios[[4]] <- list()
escenarios[[4]]$S <- 10
escenarios[[4]]$nclust <- 2


# Tres cluster con series bilineares
# primer cluster: X[t] <- p1(0.1 <-> 0.3) * X[t-1] - p2(0.1 <-> 0.3) * epsilon[t-1] * X[t-1] + epsilon[t]
# segundo cluster: X[t] <- p1(0.4 <-> 0.6) * X[t-1] - p2(0.4 <-> 0.6) * epsilon[t-1] * X[t-1] + epsilon[t]
# segundo cluster: X[t] <- p1(0.7 <-> 0.9) * X[t-1] - p2(0.7 <-> 0.9) * epsilon[t-1] * X[t-1] + epsilon[t]
# T=200
R_series5 <- f_escenario.5(R=1, S=10, nclust=3, T=50, 
                           C1.p1.min=0.1, C1.p1.max=0.2, C1.p2.min=0.1, C1.p2.max=0.2,
                           C2.p1.min=0.4, C2.p1.max=0.5, C2.p2.min=0.4, C2.p2.max=0.5, 
                           C3.p1.min=0.7, C3.p1.max=0.8, C3.p2.min=0.7, C3.p2.max=0.8)


escenarios[[5]] <- list()
escenarios[[5]]$S <- 10
escenarios[[5]]$nclust <- 3


# escenario 5 pero con T= 100
#R_series8 <- f_escenario.5(R=1, S=10, nclust=3, T=500, 
#                           C1.p1.min=0.1, C1.p1.max=0.2, C1.p2.min=0.1, C1.p2.max=0.2,
#                           C2.p1.min=0.4, C2.p1.max=0.5, C2.p2.min=0.4, C2.p2.max=0.5, 
#                           C3.p1.min=0.7, C3.p1.max=0.8, C3.p2.min=0.7, C3.p2.max=0.8)

R_series8 <- f_escenario.6(R=1, S=10, nclust=2, T=50,
                           C1.p1.min=0.1, C1.p1.max=0.2, C1.p2.min=0.1, C1.p2.max=0.2,
                           C2.p1.min=0.4, C2.p1.max=0.5, C2.p2.min=0.4, C2.p2.max=0.5,
                           Extra.p1.min=0.3, Extra.p1.max=0.3, Extra.p2.min=0.3, Extra.p2.max=0.3)

escenarios[[8]] <- list()
escenarios[[8]]$S <- 10
escenarios[[8]]$nclust <- 2 #3



# Dos cluster con series bilineares + una serie a caballo
# primer cluster: X[t] <- p1(0.1 <-> 0.4) * X[t-1] - p2(0.1 <-> 0.4) * epsilon[t-1] * X[t-1] + epsilon[t]
# segundo cluster: X[t] <- p1(0.6 <-> 0.9) * X[t-1] - p2(0.6 <-> 0.9) * epsilon[t-1] * X[t-1] + epsilon[t]
# serie a caballo : X[t] <- p1=0.5 * X[t-1] - p2=0.5 * epsilon[t-1] * X[t-1] + epsilon[t]
# T=200
R_series6 <- f_escenario.6(R=100, S=10, nclust=2, T=200,
                           C1.p1.min=0.4, C1.p1.max=0.5, C1.p2.min=0.4, C1.p2.max=0.5,
                           C2.p1.min=0.7, C2.p1.max=0.8, C2.p2.min=0.7, C2.p2.max=0.8,
                           Extra.p1.min=0.6, Extra.p1.max=0.6, Extra.p2.min=0.6, Extra.p2.max=0.6)


escenarios[[6]] <- list()
escenarios[[6]]$S <- 10
escenarios[[6]]$nclust <- 2

# tres custers de series no lineales (diferentes modelos)
# primer cluster: X[t] <- (p1(0.1 <-> 0.4) - 10 * exp(-X[t-1]^2)) * X[t-1] + epsilon[t] 
# segundo cluster: X[t] <- p1(0.6 <-> 0.9) * X[t-1] - p2(0.6 <-> 0.9) * epsilon[t-1] * X[t-1] + epsilon[t]
# tercer cluster: X[t] <- p1(0.1 <-> 0.4) * abs(X[t-1]) * (3 + abs(X[t-1]))^-1 + epsilon[t]
# T=200
R_series7 <- f_escenario.7(R=1, S=10, nclust=3, T=50,
                           C1.p1.min=0.1, C1.p1.max=0.4, 
                           C2.p1.min=0.4, C2.p1.max=0.6, C2.p2.min=0.4, C2.p2.max=0.6, 
                           C3.p1.min=0.1, C3.p1.max=0.4)


escenarios[[7]] <- list()
escenarios[[7]]$S <- 10
escenarios[[7]]$nclust <- 3


test <- list()

test[[1]] <- list()
test[[1]][[1]] <- R_series1$datos
test[[1]][[2]] <- R_series1$R
test[[1]][[3]] <- R_series1$true_labels
test[[1]][[4]] <- R_series1$C


test[[2]] <- list()
test[[2]][[1]] <- R_series2$datos
test[[2]][[2]] <- R_series2$R
test[[2]][[3]] <- R_series2$true_labels
test[[2]][[4]] <- R_series2$C

test[[3]] <- list()
test[[3]][[1]] <- R_series3$datos
test[[3]][[2]] <- R_series3$R
test[[3]][[3]] <- R_series3$true_labels
test[[3]][[4]] <- R_series3$C


test[[4]] <- list()
test[[4]][[1]] <- R_series4$datos
test[[4]][[2]] <- R_series4$R
test[[4]][[3]] <- R_series4$true_labels
test[[4]][[4]] <- R_series4$C

test[[5]] <- list()
test[[5]][[1]] <- R_series5$datos
test[[5]][[2]] <- R_series5$R
test[[5]][[3]] <- R_series5$true_labels
test[[5]][[4]] <- R_series5$C


test[[6]] <- list()
test[[6]][[1]] <- R_series6$datos
test[[6]][[2]] <- R_series6$R
test[[6]][[3]] <- R_series6$true_labels
test[[6]][[4]] <- R_series6$C


test[[7]] <- list()
test[[7]][[1]] <- R_series7$datos
test[[7]][[2]] <- R_series7$R
test[[7]][[3]] <- R_series7$true_labels
test[[7]][[4]] <- R_series7$C

test[[8]] <- list()
test[[8]][[1]] <- R_series8$datos
test[[8]][[2]] <- R_series8$R
test[[8]][[3]] <- R_series8$true_labels
test[[8]][[4]] <- R_series8$C

results <- list()
times <- list()

for (test_i in 1:length(test)){
  R_series <- test[[test_i]][[1]]
  R <- test[[test_i]][[2]]
  true_labels <- test[[test_i]][[3]]
  nclust <- test[[test_i]][[4]]
  
  #-----------------------------------------------------------#
  #---  CALCULAMOS MATRICES DE DISTANCIAS-CARACTERISTICAS  ---#
  #-----------------------------------------------------------#
  
  # las metricas que trabajan con los datos en crudo tienen 
  # como matriz de características los propios datos ( en el código: x <- Y )

  R_X_eucl = list()    # free_model
  R_D_eucl = list()
  for (i in 1:R){
    Y <- R_series[[i]]
    x <- Y                               #features
    D <- as.matrix(diss(x, "EUCL"))      #dissimilaritys
    R_X_eucl[[i]] <- x
    R_D_eucl[[i]] <- D
  }
  
  
  R_X_dtw = list()    # free_model
  R_D_dtw = list()
  for (i in 1:R){
    Y <- R_series[[i]]
    x <- Y                               #features
    D <- as.matrix(diss(x, "DTW"))       #dissimilaritys
    R_X_dtw[[i]] <- x
    R_D_dtw[[i]] <- D
  }
  
  #R_X_cort_k2 = list()    # free_model
  #R_D_cort_k2 = list()
  ###for (i in 1:R){
  #  Y <- R_series[[i]]
  #  x <- Y                                #features
  #  D <- as.matrix(diss(x, "CORT", k=2))       #dissimilaritys
  #  R_X_cort_k2[[i]] <- x
  #  R_D_cort_k2[[i]] <- D
  #}
  
  #R_X_cort_k4 = list()    # free_model
  #R_D_cort_k4 = list()
  #for (i in 1:R){
  #  Y <- R_series[[i]]
  #  x <- Y                                #features
  #  D <- as.matrix(diss(x, "CORT", k=4))       #dissimilaritys
  #  R_X_cort_k4[[i]] <- x
  #  R_D_cort_k4[[i]] <- D
  #}
  
  
  # matriz de caracteristicas de ACF
  R_X_acf = list()
  R_D_acf = list()
  for (i in 1:R){
    Y <- R_series[[i]]
    x <- Acf(Y)                          #features
    D <- as.matrix(diss(x, "EUCL"))      #dissimilaritys
    rownames(x) <- rownames(Y)  # asignar nombres de filas
    R_X_acf[[i]] <- x
    R_D_acf[[i]] <- D
  }
  
  
  # matriz de caracteristicas de PACF 
  R_X_pacf = list()
  R_D_pacf = list()
  for (i in 1:R){
    Y <- R_series[[i]]
    x <- Pacf(Y)                         #features
    D <- as.matrix(diss(x, "EUCL"))      #dissimilaritys
    rownames(x) <- rownames(Y)  # asignar nombres de filas
    R_X_pacf[[i]] <- x
    R_D_pacf[[i]] <- D
  }
  
  
  # matriz de caracteristicas de Piccolo 
  R_X_picc = list()
  R_D_picc = list()
  for (i in 1:R){
    Y <- R_series[[i]]
    x <- pic(Y)                          #features
    D <- as.matrix(diss(x, "EUCL"))      #dissimilaritys
    rownames(x) <- rownames(Y)  # asignar nombres de filas
    R_X_picc[[i]] <- x
    R_D_picc[[i]] <- D
  }
  
  # matriz de caracteristicas de QAF 
  R_X_qaf = list()
  R_D_qaf = list()
  for (i in 1:R){
    Y <- R_series[[i]]
    x <- t(apply(Y,1,qaf,tau=c(0.1,0.5,0.9),L=1))   #features
    D <- as.matrix(diss(x, "EUCL"))                 #dissimilaritys
    R_X_qaf[[i]] <- x
    R_D_qaf[[i]] <- D
  }
  
  
  
  
  #-------------------------#
  #---  SIMULACION FINAL ---#
  #-------------------------#
  
  #metrics <- c("EUCL", "DTW", "CORT_k2", "CORT_k4", "ACF", "PACF", "PICCOLO", "QAF")
  metrics <- c("EUCL", "DTW", "ACF", "PACF", "PICCOLO", "QAF")
  
  D_matrixs = list()
  D_matrixs[[1]] <- R_D_eucl
  D_matrixs[[2]] <- R_D_dtw
  #D_matrixs[[3]] <- R_D_cort_k2
  #D_matrixs[[4]] <- R_D_cort_k4
  D_matrixs[[3]] <- R_D_acf
  D_matrixs[[4]] <- R_D_pacf
  D_matrixs[[5]] <- R_D_picc
  D_matrixs[[6]] <- R_D_qaf
  
  X_matrixs = list()
  X_matrixs[[1]] <- R_X_eucl
  X_matrixs[[2]] <- R_X_dtw
  #X_matrixs[[3]] <- R_X_cort_k2
  #X_matrixs[[4]] <- R_X_cort_k4
  X_matrixs[[3]] <- R_X_acf
  X_matrixs[[4]] <- R_X_pacf
  X_matrixs[[5]] <- R_X_picc
  X_matrixs[[6]] <- R_X_qaf
  
  
  
  fcmd.results <- vector("list", length = length(m_values))
  fcmd.time <- vector("list", length = length(m_values))
  
  tiempo_i <- proc.time()
  for (i in 1:length(m_values)) {
    m_value <- m_values[i]
    m_value_str <- paste('m=',as.character(m_value))
    fcmd.results[[m_value_str]] <- list()
    fcmd.time[[m_value_str]] <- list()
    
    for (m in 1:length(metrics)) {
      metric <- metrics[m]
      fcmd.results[[m_value_str]][[metric]] <- list()
      fcmd.time[[m_value_str]][[metric]] <- list()
      
      tiempo_inicio <- proc.time()
      for (r in 1:R) {
        res <- FCMdC_QAF(Y = R_series[[r]],D = D_matrixs[[m]][[r]] ,X = X_matrixs[[m]][[r]],C=nclust,m=m_value,max.iter=20,aleat=FALSE)
        res$metric <- metric
        res$true_labs <- true_labels
        
        membership <- res$U
        
        res$Ari.f <- ARI.F(VC= true_labels,U=membership)
        res$Jaccard.f <- JACCARD.F(VC=true_labels,U=membership)
        
        fcmd.results[[m_value_str]][[metric]][[r]] <- res
      }
      tiempo_final <- proc.time()
      fcmd.time[[m_value_str]][[metric]]$tiempo <- tiempo_final - tiempo_inicio
    }
  }
  results[[test_i]] <- fcmd.results
  times[[test_i]] <- fcmd.time
  
  tiempo_f <- proc.time()
  tiempo_f - tiempo_i
}


