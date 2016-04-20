## funcao para calcular os autovalores da jacobiana
## no equilibrio a partir de uma funcao generica de ode, como a usada pelo odesolve
## neq= vetor de abundancias no equilibrio
autov <- function(func,neq,pars){
  f1 <- function(t=0, y, parms=NULL, ...){
    func(time=0,init=y,parms=pars)
  }
  eigen(jacobian.full(y=neq,func=f1),only.values=TRUE)$values
}

## Calcula a matriz de comunidades para um sistema LV de competidores
jacob.lv <- function(r1,K1,r2,K2,alfa,beta){
  f1 <- function(t=0, y, parms=NULL, ...){
    LVC(time=0,init=y,parms=pars)
  }
  pars <- c(r1=r1,r2=r2,k1=K1,k2=K2,a=alfa,b=beta)
  neq <- lv.neq(r1=r1,K1=K1,r2=r2,K2=K2,alfa=alfa,beta=beta)
  jacobian.full(y=neq,func=f1)
}

## funcao para pertubar uma equacao logistica perto de um ponto e plotar o grafico
plota.logist <- function(n=1,r=0.1,K=50,time=150,perturb=0,t.perturb=80,int,...){
  tempo <- seq(0,time,int)
  init <- c(n=n)
  parms <- c(r=r,k=K)
  if(perturb==0)out3 <- ode(func=logist,y=init,parms=parms,time=tempo)
  else{
    out <- ode(func=logist,y=init,parms=parms,time=tempo[tempo<t.perturb])
    out2 <- ode(func=logist,y=c(n=max(out[,2])+perturb),parms=parms,time=seq(0,max(tempo-t.perturb),by=min(diff(tempo))))
    out2[,1] <- out2[,1]+max(out[,1])
    out3 <- rbind(out,out2)
  }
  plot(out3,main="", xlab="Tempo", ylab="N", type="l", ylim=c(0,max(out3[,2])),lwd=1.5, col="blue")
}

## Plota e perturba a logistica com efeito allee
plota.allee <- function(n=1,r=0.1,K=50,a=10,time=150,perturb=0,t.perturb=80,int,pp=TRUE,...){
  tempo <- seq(0,time,int)
  init <- c(n=n)
  parms <- c(r=r,k=K,a=a)
  if(perturb==0)out3 <- ode(func=allee,y=init,parms=parms,time=tempo)
  else{
    out <- ode(func=allee,y=init,parms=parms,time=tempo[tempo<t.perturb])
    out2 <- ode(func=allee,y=c(n=max(out[,2])+perturb),parms=parms,time=seq(0,max(tempo-t.perturb),by=min(diff(tempo))))
    out2[,1] <- out2[,1]+max(out[,1])
    out3 <- rbind(out,out2)
  }
  if(pp){
  plot(out3,main="", xlab="Tempo", ylab="N", type="l", ylim=c(0,max(out3[,2])),lwd=1.5, col="blue")
}
  if(!pp){
  plot(out3,...)
}
}


##funcao que retorna os valores em equilibrio para um sistema LV de competidores
lv.neq <- function(r1,K1,r2,K2,alfa,beta){
  n1 <- (K1-alfa*K2)/(1-alfa*beta)
  n2 <- (K2-beta*K1)/(1-alfa*beta)
  c(n1=n1,n2=n2)
  }

##Funcao para plotar sistema de competidores LV
plota.LV <- function(n=c(n1=1,n2=1),r1=0.2,K1=150,r2=0.2,K2=100,alfa=0.2,beta=0.1,time=150,perturb=c(0,0),t.perturb=80,...){
  tempo <- seq(0,time,0.5)
  init <- n
  parms <- c(r1=r1,r2=r2,k1=K1,k2=K2,a=alfa,b=beta)
  if(all(perturb==0)) out3 <- ode(func=LVC,y=init,parms=parms,time=tempo)
  else{
    out <- ode(func=LVC,y=init,parms=parms,time=tempo[tempo<t.perturb])
    init2 <- out[nrow(out),2:3]+perturb
    out2 <- ode(func=LVC,y=init2,parms=parms,time=seq(0,max(tempo-t.perturb),by=min(diff(tempo))))
    out2[,1] <- out2[,1]+max(out[,1])
    out3 <- rbind(out,out2)
  }
    plot.LV(out3,parms,main="",lwd=1.5,...)
}

## funcao para plotar o resultado da integracao numerica do sistema LV
plot.LV <- function(out,parms,legend=TRUE,...){
matplot(out[,1],out[,2:3],type="l",lty=1,col=c("blue","red"), ylim=c(min(out[,2:3]),max(parms[3:4])),xlab="tempo",ylab="N",...)
if(!legend)legend("topleft",c("N1","N2"),lty=1,col=c("blue","red"))
}


## funcao para pertubar um sistema LV perto do equilibrio e plotar o grafico
perturba.LV <- function(deltas,func=LVC,parms,tempo=seq(0,500,by=0.5),...){
  y <- neq(parms)
  out <- ode(func=func,y=y,parms=parms,time=seq(0,50,by=0.5))
  out2 <- ode(func=func,y=y+deltas,parms=parms,time=tempo)
  out2[,1] <- out2[,1]+max(out[,1])
  out3 <- rbind(out,out2)
  plot.LV(out3,parms,...)
}

## equacao logistica
logist <- function(time, init, parms){
  with(as.list(c(init,parms)),
       {
         dn <- r*n*((k-n)/k)
         list(c(dn))
       })
}

## equacao logistica com efeito Alee
allee <- function(time, init, parms){
  with(as.list(c(init,parms)),
       {
         dn <- r*n*(1-n/k)*((n-a)/k)
         list(c(dn))
       })
}


## Funcao com sistema competidores Lotka-Volterra
LVC <- function(time, init, parms){
  with(as.list(c(init,parms)),
       {
         dn1 <- r1*n1*((k1-n1-a*n2)/k1)
         dn2 <- r2*n2*((k2-n2-b*n1)/k2)
         list(c(dn1, dn2))
       })
}

## funcao para a simulacao do May
may <- function(S,C,f,nsim=100){
  m <- diag(rep(-1,S))
  ind <- which(m==0,arr.ind=TRUE)
  n <- round((C*S^2)-S,0)
  n.estavel <- 0
  for(i in 1:nsim){
    vals <- rnorm(n,mean=0,sd=f)
    ri <- sample(1:nrow(ind),size=n)
    m2 <- m
    m2[ind[ri,]] <- vals
    autov <- eigen(m2,only.values=TRUE)$values
    n.estavel <- n.estavel+all(Re(autov)<0)
  }
  data.frame(S=S,C=C,f=f,p.estab=n.estavel/nsim)
}

##Simulacao que sorteia matrizes de comunidades entre range de riquezas e conectancias
## e registra o maior valor dos autovalores
yodzis <- function(Smin,Smax,Cmin,Cmax,f=0.2,nsim=100){
  results <- matrix(ncol=4,nrow=nsim,dimnames=list(NULL,c("S","C","f","max.aut")))
  for(i in 1:nsim){
    S <- sample(Smin:Smax,size=1)
    C <- runif(1,Cmin,Cmax)
    m <- diag(rep(-1,S))
    ind <- which(m==0,arr.ind=TRUE)
    n <- round((C*S^2)-S,0)
    vals <- rnorm(n,mean=0,sd=f)
    ri <- sample(1:nrow(ind),size=n)
    m[ind[ri,]] <- vals
    max.autov <- max(Re(eigen(m,only.values=TRUE)$values))
    results[i,] <- c(S,C,f,max.autov)
  }
  as.data.frame(results)
}
