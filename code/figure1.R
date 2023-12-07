##plots for figure 1 of vibe check
##functions 

#Helper functions

logit <- function(theta){log(theta/(1-theta))}
expit <- function(lodds){exp(lodds)/(exp(lodds)+1)}
odds  <- function(theta){theta/(1-theta)}
probs <- function(odds){odds/(odds+1)}

#Posterior distribution of each parameterization
#Theta
post.theta <- function(theta,x=1,n=1){
  dens <- gamma(2+n)/gamma(1+x)/gamma(n-x+1)*
    theta^(x)*(1-theta)^(n-x)
  return(dens)
}

#Log-odds
post.lodds <- function(lodds,x=1,n=1){
  dens <- gamma(2+n)/gamma(1+x)/gamma(n-x+1)*
    exp(lodds*(x+1))/(1+exp(lodds))^(n+2)
  return(dens)
}

#Odds
post.odds <- function(odds,x=1,n=1){
  dens <- gamma(2+n)/gamma(1+x)/gamma(n-x+1)*
    odds^x/(1+odds)^(n+2)
  return(dens)
}

#collect in a list
posteriors <- list("theta"=post.theta, "lodds" = post.lodds, "odds"=post.odds)

#find the MAP for each parameterization
MAP <- function(x=1,n=1){
  map.theta <- x/n
  map.lodds <- log((x+1)/(n-x+1))
  map.odds <- x/(n-x+2)
  return(list("theta"=map.theta,"lodds"=map.lodds,
              "odds"=map.odds))
}

#Find the fisher information for each
#This will help with the HDI algorithm (only used for log-odds currently)
f.info <- function(x,n){
  map <- MAP(x,n)
  f.theta <- n/(map$theta*(1-map$theta))
  f.lodds <- n*exp(map$lodds)/(1+exp(map$lodds))^2
  f.odds  <- n*(map$odds*(1+map$odds)^2)^(-1)
  return(list("theta" = f.theta, "lodds" = f.lodds,
              "odds"=f.odds))
}

#Find HDI for theta
HDI.theta <- function(x=1,n=1,p=.95){
  f <- posteriors$theta
  int <- c(0,1)
  map = MAP(x,n)$theta
  map.dens <- f(map,x,n)
  delta <- .001*map.dens
  eps <- map.dens - delta
  k <- 0
  i <- 1
  roots<-c(0,0)
  while(k < p ){
    if(map==0){ #If MAP = 0, then lower end of HDI is 0
      roots[1] <- 0
    } else roots[1] <- uniroot(function(a)f(a,x,n)-eps,
                               interval = c(int[1],map))$root
    if(map==1){ #If MAP = 1, then upper end of HDI is 1
      roots[2] <- 1
    }else roots[2] <- uniroot(function(a)f(a,x,n)-eps,
                              interval = c(map,int[2]))$root
    k <- integrate(function(b)f(b,x,n), lower=roots[1],upper=roots[2])$value
    eps <- eps - delta
    i <- i+1
  }
  return(HDI=roots)
}

HDI.lodds <- function(x=1,n=1,p=.95){
  f <- posteriors$lodds
  map <- MAP(x,n)$lodds
  info <- f.info(x,n)$lodds
  int <- c(map-2.5/sqrt(info), map+2.5/sqrt(info))
  map.dens <- f(map,x,n)
  delta <- .001*map.dens
  eps <- map.dens - delta
  k <- 0
  i <- 1
  roots<-c(0,0)
  while(k < p ){
    roots[1] <- uniroot(function(a)f(a,x,n)-eps,
                        interval = c(int[1],map))$root
    roots[2] <- uniroot(function(a)f(a,x,n)-eps,
                        interval = c(map,int[2]))$root
    k <- integrate(function(a)f(a,x,n), lower=roots[1],upper=roots[2])$value
    eps <- eps - delta
    i <- i+1
  }
  return("HDI"=roots)
}

HDI.odds <- function(x=1,n=1,p=.95){
  f <- posteriors$odds
  int <- c(0,1e3)
  map <- MAP(x,n)$odds
  map.dens <- f(map,x,n)
  delta <- .001*map.dens
  eps <- map.dens - delta
  k <- 0
  i <- 1
  roots<-c(0,0)
  while(k < p ){
    if(map==0){ #If MAP = 0 (unlikely, but I think possible), then lower end of HDI is 0
      roots[1] <- 0
    }else roots[1] <- uniroot(function(a)f(a,x,n)-eps,
                              interval = c(int[1],map))$root
    roots[2] <- uniroot(function(a)f(a,x,n)-eps,
                        interval = c(map,int[2]))$root
    k <- integrate(function(a)f(a,x,n), lower=roots[1],upper=roots[2])$value
    eps <- eps - delta
    i <- i+1
  }
  return("HDI"=roots)
  
}

##plotting

x <- 32
n <- 47

tlocx = .1
tlocy = .9

map <- MAP(x,n) 
ROPE <- c(.47,.53)
map
rbind(HDI.theta(x,n),HDI.lodds(x,n),HDI.odds(x,n))
test <- rbind(ROPE,HDI.theta(x,n),expit(HDI.lodds(x,n)),probs(HDI.odds(x,n)))
rownames(test) <- c("ROPE","Prob scale", "Logit scale", "Odds scale")
colnames(test) <- c("Lower", "Upper")
round(test,3)

pdf(file='figure1a.pdf',width=6,height=4)
par(mfrow=c(1,1))
plot(1,axes=F,
     type="n",lwd=2, xlab="Probability of success",xlim=c(0,1),
     ylab="Density",ylim=c(0,post.theta(map$theta,x,n)+.1),cex.lab=1.4)
#abline(v=map$theta,lty=2,lwd=2)
rect(xleft = .47, xright = .53, ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("lightgrey"))
lines(seq(0,1,.01),post.theta(seq(0,1,.01),x=x,n=n))
abline(v=HDI.theta(x,n),lty=2)
abline(v=ROPE,lty=1,col="black")
box()
axis(1,cex.axis=1.4)
axis(2,cex.axis=1.4)
#abline(v=HDI.theta(x,n),lty=3,col="red",lwd=2)
mtext("(a)", side = 3, line = -2, adj = .04, cex = 1.4)

dev.off()

pdf(file='figure1b.pdf',width=6,height=4)
plot(1,axes=F,
     type="n",lwd=2, xlab= "Log-odds of success",
     ylab="Density",xlim=c(-2,2),ylim=c(0,post.lodds(map$lodds,x,n)),cex.lab=1.4)
rect(xleft = logit(.47), xright = logit(.53), ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("lightgrey"))
lines(seq(-2,2,.01),post.lodds(seq(-2,2,.01),x=x,n=n))
abline(v=HDI.lodds(x,n),lty=2)
abline(v=logit(ROPE),lty=1,col="black")
box()
axis(1,cex.axis=1.4)
axis(2,cex.axis=1.4)
#abline(v=map$lodds,lty=2,lwd=2)
#abline(v=HDI.lodds(x,n),lty=3,col="red",lwd=2)
mtext("(b)", side = 3, line = -2, adj = .04, cex = 1.4)
dev.off()

pdf(file='figure1c.pdf',width=6,height=4)
plot(1,
     type="n",lwd=2, xlab="Odds of success",
     ylab="Density",xlim=c(0,6),ylim=c(0,post.odds(map$odds,x,n)),axes=F,cex.lab=1.4)  
rect(xleft = odds(.47), xright = odds(.53), ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("lightgrey"))
lines(seq(0,6,.01),post.odds(seq(0,6,.01),x=x,n=n))
abline(v=HDI.odds(x,n),lty=2)
abline(v=odds(ROPE),lty=1,col="black")
box()
axis(1,cex.axis=1.4)
axis(2,cex.axis=1.4)
mtext("(c)", side = 3, line = -2, adj = .04, cex = 1.4)
dev.off()

