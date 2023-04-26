# Rachid Benabdelouahed
# Baptiste Gorteau

install.packages("FactoMineR")
install.packages("plotrix")
install.packages("dplyr")
library(dplyr)
library(FactoMineR)
library(plotrix)

# paramètres graphiques

par(font.lab = 2)

# Importation des données
alcool <- read.csv("alcool.csv",sep=";",header = TRUE, row.names = 1)
alcool

# Question 1 

# fonction cr prof


ecartype <- function(x){
  sqrt(mean((x-mean(x))^2))
}

cr <- function(matrice){
  n <- nrow(matrice)
  moy <- apply(matrice, FUN = mean, MARGIN = 2)
  moy <- matrix(rep(moy, n), nrow = n, byrow = T)
  
  e_type <- apply(matrice, FUN = ecartype, MARGIN = 2)
  e_type <- matrix(rep(e_type,n), nrow = n, byrow = T)
  return((matrice-moy)/e_type)
}

# Question 2

# Création fonction ACP

acp <- function(tablecr){
  mat_cor <- cor(tablecr)
  mat_cov <- cov(tablecr)
  #calcul des inerties
  inertie <- (eigen(mat_cov)$values/sum(eigen(mat_cov)$values))*100
  inertie <- t(as.matrix(inertie))
  colnames(inertie) <- c(1:length(inertie))
  row.names(inertie) <- "Inertie des axes (en %)"
  # Coordonnées des variables
  mat_eig <- eigen(mat_cor)
  mat_val_propres <- NULL
  for(j in 1:ncol(tablecr)){
    rf <- sqrt(mat_eig$values[j])*mat_eig$vectors[,j]
    mat_val_propres <- cbind(mat_val_propres,rf)
  }
  row.names(mat_val_propres) <- colnames(tablecr)
  coord_var <- mat_val_propres[,1:3]
  colnames(coord_var) <- c("axe1","axe2","axe3")
  # Coordonnées des individus
  axe1 <- 0
  axe2 <- 0
  axe3 <- 0
  for (i in 1:ncol(tablecr)){
    axe1 <- axe1+tablecr[,i]*eigen(mat_cor)$vectors[i,1]
    axe2 <- axe2+tablecr[,i]*eigen(mat_cor)$vectors[i,2]
    axe3 <- axe3+tablecr[,i]*eigen(mat_cor)$vectors[i,3]
  }
  coord_indiv <- data.frame(axe1,axe2,axe3)
  row.names(coord_indiv) <- rownames(tablecr)
  return(list("inertie" =inertie ,"var"=coord_var,"indiv"=coord_indiv))
}

barplot(acp(cr(alcool))$inertie)
cr(alcool)

#Normé à la valeur propre et diviser par n du côté des variables

# Question 3 

plotacp <- function(X,axe1,axe2,type){
  if(type=="ind"){
    plot(X$indiv[,axe1],X$indiv[,axe2],pch=20,main="Graphique de l'ACP \n pour les individus",xlab="",ylab="",xlim=c(min(X$ind[,axe1])-0.5,max(X$ind[,axe1])+0.5),ylim=c(min(X$ind[,axe2])-0.5,max(X$ind[,axe2])+0.5))
    title(xlab = paste("Dim",axe1,"(",round(X$inertie[1,axe1],2),"%)"),ylab =paste("Dim",axe2,"(",round(X$inertie[1,axe2],2),"%)"),adj=1,line=2.2)
    grid(NULL,NULL,lty=2,lwd=1)
    abline(h=0)
    abline(v=0)
    text(X$indiv[,axe1],X$indiv[,axe2],labels=row.names(X$indiv),cex=0.9,pos=1)
  }
  if(type=="var"){
    plot(X$var[,axe1],X$var[,axe2],xlim=c(-1,1),ylim=c(-1,1),type="n",main="Graphique de l'ACP \n pour les variables",xlab="",ylab="")
    title(xlab = paste("Dim",axe1,"(",round(X$inertie[1,axe1],2),"%)"),ylab =paste("Dim",axe2,"(",round(X$inertie[1,axe2],2),"%)"),adj=1,line=2.2)
    grid(NULL,NULL,lty=2,lwd=1)
    segments(-1,0,1,0,lty=2)
    segments(0,-1,0,1,lty=2)
    draw.circle(0,0,1)
    text(X$var[,axe1],X$var[,axe2],labels=row.names(X$var),cex=0.9,pos=1)
    for(i in 1:nrow(X$var)){
      arrows(0,0,X$var[i,axe1],X$var[i,axe2],length = 0.1)
    }
  }
}

plotacp(acp(cr(alcool)),1,2,"var")

# Question 4

indsup <- function(data_frame,ind_sup){
  table_ind_sup <- data_frame[ind_sup,]
  n <- nrow(table_ind_sup)
  data_frame_2 <- data_frame[-ind_sup,]
  moy_table <- apply(data_frame_2,FUN=mean,MARGIN = 2)
  moy_table <-  matrix(rep(moy_table, n), nrow = n, byrow = T)
  ec_type_table <- apply(data_frame_2,FUN=ecartype, MARGIN = 2)
  ec_type_table <- matrix(rep(ec_type_table,n), nrow = n, byrow = T)
  ind_sup_cr <- (table_ind_sup-moy_table)/ec_type_table
  axe1 <- 0
  axe2 <- 0
  axe3 <- 0
  mat_vect_propre <- eigen(cov(cr(data_frame_2)))$vectors
  for (i in 1:ncol(table_ind_sup)){
    axe1 <- axe1+ind_sup_cr[,i]*mat_vect_propre[i,1]
    axe2 <- axe2+ind_sup_cr[,i]*mat_vect_propre[i,2]
    axe3 <- axe3+ind_sup_cr[,i]*mat_vect_propre[i,3]
  }
  coord_indiv <- data.frame(axe1,axe2,axe3)
  row.names(coord_indiv) <- row.names(data_frame[ind_sup,])
  return(coord_indiv)
}

indsup(alcool,c(1,2))
