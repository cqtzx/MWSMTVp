library(compiler) 
library(imager) 
library(ThreeWay)
library(jpeg)
library(png)
library(pracma) 
library(signal) 
library(bmp)
library(pixmap) 
library(rpca) 
library(redR) 
library(caTools) 
library(rmatio) 
library(ggplot2)
library(fBasics)
library(R.matlab)

F2norm<-cmpfun(function(M)sqrt(sum(M^2)))

solve_Lp <- function( y, lambda, p ){
  J     =   2;
  tau   =  (2*lambda*(1-p))^(1/(2-p))+p*lambda*(2*(1-p)*lambda)^((p-1)/(2-p));
  x     =   matrix(0,NROW(y),ncol(y))
  i0    =   which( abs(y)>tau );
  if (length(i0)>=1){
    y0    =   y[i0];
    t     =   abs(y0);
    for  (j in 1:J){
      t    =  abs(y0) - p*lambda*(t)^(p-1);
    }
    x[i0]   =  sign(y0)*t;
  }
  x
}


Dg <- function(x){
  m <- NROW(x);n <- NCOL(x)
  grad <- array(0,c(m,n,2))
  grad[,,1] = x - circshift(x,c(-1,0))
  grad[m,,1] = 0
  grad[,,2] = x - circshift(x,c(0,-1))
  grad[,n,2] = 0
  grad
}

dtwg <- function(grad,p){
  g1 <- grad[,,1];g2 <- grad[,,2]
  m <- NROW(g1);n <- NCOL(g1)
  f <- matrix(0,m+1,n+1);g1 <- cbind(0,rbind(0,g1));g2 <- cbind(0,rbind(0,g2))
  for(i in 2:(m+1)){
    for(j in 2:(n+1)){
      g1[(m+1),j] <- g2[i,(n+1)] <- 0
      f[i,j] <- ((abs(g1[i,j])+0.001)^(p-1))*g1[i,j]-((abs(g1[i-1,j])+0.001)^(p-1))*g1[i-1,j]+((abs(g2[i,j])+0.001)^(p-1))*g2[i,j]-((abs(g2[i,j-1])+0.001)^(p-1))*g2[i,j-1]
      
    }
  }
  f[-1,-1]
}

WFGP_2D <- function(y,lambda,n_iters,p){
  n1 <- NROW(y);n2 <- NCOL(y)
  grad_next <- array(0,c(n1,n2,2))
  grad_prev <- array(0,c(n1,n2,2))
  u <- array(0,c(n1,n2,2))
  t_prev <- 1
  for(i in 1:n_iters){
    grad_next = u + 1/8/lambda*Dg(y - lambda*dtwg(u,p)) 
    deno = array(0,c(n1,n2,2))
    deno[,,1] = pmax(1,abs(grad_next[,,1]))
    deno[,,2] = pmax(1,abs(grad_next[,,2]))
    grad_next = grad_next/deno
    t_next = (1+sqrt(1+4*t_prev^2))/2
    u = grad_next + (t_prev-1)/t_next*(grad_next-grad_prev)
    grad_prev = grad_next
    t_prev = t_next
  }
  x = y - lambda*dtwg(grad_next,p)   
  x
}

solpw <- function(y,lambda,p){
  K <- 2
  tau <- (2*lambda*(1-p))^(1/(2-p))+p*lambda*(2*(1-p)*lambda)^((p-1)/(2-p))
  x <- rep(0,length(y))
  i0 <- which(abs(y)>tau)
  svp <- length(i0)
  if(length(i0)>=1){
    y0 <- y[i0]
    t <- abs(y0)
    lambda0 <- lambda[i0]
    for(j in 1:K){
      t <- abs(y0)-p*lambda0*t^(p-1)
    }
    x[i0] <- sign(y0)*t
  }
  A <- list(x,svp)
  A
}

wscp <- function(G,p){
  n <- dim(G)[2]
  svdg <- svd(G)
  Q <- svdg$u; sigma <- svdg$d; R <- svdg$v
  NSig <- sigma[n]
  C <- 2*sqrt(2)*(NSig^2)
  Temp <- sqrt(pmax((sigma)^2-n*NSig^2,0))
  for(i in 1:3){
    w_Vec <- (C*sqrt(n)*NSig^2)/( Temp^(1/p) + 1e-6 )
    out <- solpw(sigma,w_Vec,p)
    s1 <- out[[1]]
  }
  Delta <- diag(s1)
  Q%*%Delta%*%t(R) 
}

Nee <- function(D,lambda,eta,mu,l1,l2,p){
  di <- dim(D);n1 <- di[1];n2 <- di[2];n3 <- di[3];
  Y1 <- Y2 <- E  <- H <- array(0,di)
  X <- D-E
  L <- orth(matrix(rnorm(n1*l1),n1,l1))
  R <- orth(matrix(rnorm(n2*l2),n2,l2))
  M <- array(0,c(NCOL(L),NCOL(R),n3))
  P <- array(0,c(NROW(X),NCOL(M[,,1]),n3))
  W <- array(0,c(NCOL(X),NCOL(M[,,1]),n3))
  for(i in 1:n3){
    M[,,i] <- t(L)%*%D[,,i]%*%R
  }
  l2.error <- sc <- sc1 <- 1
  tao<-1.1;k<-1;re<-0
  while(sc > 0.00001){
    X.old <- X
    for(j in 1:n3){
      P[,,j] <- (X[,,j]+Y2[,,j]/mu)%*%R%*%t(M[,,j])
    }
    P0 <- matrix(0,n1,l2)
    for(i in 1:n3){
      P0 <- P0+P[,,i]
    }
    P0.svd <- svd(P0)
    L <- P0.svd[["u"]]%*%t(P0.svd[["v"]])
    for(j in 1:n3){
      W[,,j] <- t(t(M[,,j])%*%t(L)%*%(X[,,j]+Y2[,,j]/mu))
    }
    W0 <- matrix(0,n2,l2)
    for(i in 1:n3){
      W0 <- W0+(W[,,i])
    }
    W0.svd <- svd(W0)
    R <- W0.svd[["u"]]%*%t(W0.svd[["v"]])
    imu <- 1/mu
    for(j in 1:n3){
      M[,,j] <- wscp(t(L)%*%(X[,,j]+Y2[,,j]/mu)%*%R,p)
    }
    for(j in 1:n3){
      H[,,j] <- (D[,,j]-E[,,j]+L%*%M[,,j]%*%t(R)+Y1[,,j]/mu-Y2[,,j]/mu)/2
      X[,,j] <- WFGP_2D(H[,,j],eta/(2*mu),3,p)
    }
    limu <- lambda/mu
    for(j in 1:n3){
      E[,,j] <- solve_Lp(D[,,j]-X[,,j]+Y1[,,j]/mu,limu,p)
    }
    for(j in 1:n3){
      Y1[,,j] <- Y1[,,j]+mu*(D[,,j]-X[,,j]-E[,,j])
      Y2[,,j] <- Y2[,,j]+mu*(X[,,j]-L%*%M[,,j]%*%t(R))
    }
    re <- rbind(re,F2norm(X.old-X)/F2norm(X.old))
    l2.error <- rbind(l2.error,re)
    sc <- 0
    for(i in 1:n3){
      sc<-sc+(F2norm(D[,,i]-L%*%M[,,i]%*%t(R)-E[,,i])+F2norm(X[,,i]-L%*%M[,,i]%*%t(R)))/F2norm(D[,,i])
    }
    sc1 <- rbind(sc1,sc)
    mu <- mu*tao
    k<-k+1
    if(k>200) break
  } 
  X
}




