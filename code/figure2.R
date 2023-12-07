
library(pBrackets)

#Add alpha to a color
add.alpha <- function(col, alpha=1){#from https://www.r-bloggers.com/how-to-change-the-alpha-value-of-colours-in-r/
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


#Make a version of the tangent lines function for normal cdf
tangentLineNormals <- function(x0=0,xw=.5,mu=0,sd=1,param="normal"){
    x.coords <- c(x0-xw/2, x0+xw/2)
    x.coords2 <- c(x0-xw*3/4, x0+xw*3/4)
    m <- dnorm(x0,mu,sd)
    if(param=="exp"){
      m <- dnorm(log(x0),mu,sd)*1/x0
    }
    x.eval <- x0
    if(param=="exp"){ #for log-normal cdf
      x.eval <- log(x.eval)
    }
    y0 <- pnorm(x.eval,mu,sd)
    y  <- c(y0 + m*(x.coords[1]-x0), y0 + m*(x.coords[2]-x0))  # heights at the ends of the tangent line
    y2  <- c(y0 + m*(x.coords2[1]-x0), y0 + m*(x.coords2[2]-x0))  # heights at the ends of the tangent line
    lines(x.coords, y, col=add.alpha("red",1), lwd=3) #Draw the tangent lines
    lines(x.coords2, y2, col=add.alpha("red",1), lwd=1) #Draw the tangent lines
    points(x0,y0,cex=1,pch=16) #add points
    lines(x.coords,rep(y[1],2),lwd=1,col="black",lty=3) #horizontal bar
    lines(rep(x.coords[2],2),y,lwd=1,col="black",lty=3) #vertical bar
    brackets(x1=x.coords[2]+.1*sd, y1=y[2], x2=x.coords[2]+.1*sd, y2=y[1], ticks = .5, curvature = .2, type = 1,
             col = "black", lwd = 1, lty = 1, xpd = FALSE)
    
}

#Let's do a demo going from normal to log-normal (i.e. from z to exp(z))
mu <-  0
sd <- 1
xw <- sd/2

z <- seq(mu - 10*sd, mu + 10*sd,.001)
t.l <- c(mu-sd,mu,mu+sd) #where to place tangent lines
x.lim <- c(mu - 2.2*sd, mu + 2.2*sd)
max.density <- max(pmax(dnorm(z,mu,sd),dlnorm(exp(z),mu,sd))) #find highest density across both pdfs
y.lim <- c(0, 1.1*max.density) #put the pdfs onto the same y-scale


par(mfrow=c(1,1))

#normal cdf
pdf(file='figure2a.pdf',width=6,height=4)
plot(z,pnorm(z,mu,sd),type="n",ylab="Cumulative probability",main=expression(paste(theta, " CDF")),xlim=x.lim,lwd=3,xlab=expression(theta),cex.lab=1.4,cex.axis=1.4)
lines(z,pnorm(z,mu,sd),lwd=3,col="grey")
tangentLineNormals(t.l[1] ,xw,mu,sd)
tangentLineNormals(t.l[2] ,xw,mu,sd)
tangentLineNormals(t.l[3] ,xw,mu,sd)
mtext("(a)", side = 3, line = -2, adj = .04, cex = 1.4)
dev.off()

#log-normal cdf

pdf(file='figure2b.pdf',width=6,height=4)
plot(exp(z),pnorm(z,mu,sd),type="n",ylab="Cumulative probability",main=expression(paste(exp(theta), " CDF")),xlim=exp(x.lim),xlab=expression(exp(theta)),cex.lab=1.4,cex.axis=1.4)
lines(exp(z),pnorm(z,mu,sd),lwd=3,col="grey")
tangentLineNormals(exp(t.l[1]) ,xw,mu,sd,param="exp")
tangentLineNormals(exp(t.l[2]) ,xw,mu,sd,param="exp")
tangentLineNormals(exp(t.l[3]) ,xw,mu,sd,param="exp")
mtext("(b)", side = 3, line = -2, adj = .04, cex = 1.4)
dev.off()

#pdfs
pdf(file='figure2c.pdf',width=6,height=4)
plot(z,dnorm(z,mu,sd),type="l",xlim=x.lim,ylab="Density", main=expression(paste(theta, " PDF")),ylim=y.lim,lwd=3,col="grey",xlab=expression(theta),cex.lab=1.4,cex.axis=1.4)
points(t.l,dnorm(t.l,mu,sd),pch=16)
mtext("(c)", side = 3, line = -2, adj = .04, cex = 1.4)
dev.off()

pdf(file='figure2d.pdf',width=6,height=4)
plot(exp(z),dlnorm(exp(z),meanlog=mu,sdlog=sd),type="l",xlim=exp(x.lim),ylab="Density",main=expression(paste(exp(theta), " PDF")),ylim=y.lim,col="grey",lwd=3,xlab=expression(exp(theta)),cex.lab=1.4,cex.axis=1.4)
points(exp(t.l), dlnorm(exp(t.l),mu,sd),pch=16)
mtext("(d)", side = 3, line = -2, adj = .95, cex = 1.4)
dev.off()




#######OLD POSTER CODE#########

# 
# ###########################
# #Now we are doing thetas and subsequent transforms
# 
# theta <- seq(.005,.995,.005)
# x <- 1
# n <- 1
# 
# par(mfrow=c(1,2))
# plot(theta,pbeta(theta,x+1,n-x+1),type="l",main="cdf of theta",ylab="Cumulative probability")
# plot(theta,posteriors[["theta"]](theta,x,n),type="l", main="pdf of theta",ylab="Density")
# 
# 
# plot(logit(theta),pbeta(theta,x+1,n-x+1),type="l",main="cdf of logit(theta)",ylab="Cumulative probability")
# 
# plot(logit(theta), posteriors[["lodds"]](logit(theta), x=1,n=1),type="l",main="pdf of logit(theta)",ylab="Density")
# 
# #Tangent line maker for the probability examples
# tangentLine <- function(x0=.5, xw=.2, x=0, n=0, param="theta"){ #make sure to set x0 and xw on param scale
#     x.coords <- c(x0-xw/2, x0+xw/2)  # x-coords for tangent line
#     m <- posteriors[[param]](x0,x,n)  # slope of tangent line
#     x.eval <- x0   #for input to qbeta
#     if(param=="lodds"){
#         x.eval <- expit(x0)        #make sure x0 is on prob. parameterization for qbeta
#     }
#     if(param=="odds"){   
#         x.eval <- probs(x0)
#     }
#     y0 <- pbeta(x.eval, x+1, n-x+1)   # height at midpoint of the tangent line
#     y  <- c(y0 + m*(x.coords[1]-x0), y0 + m*(x.coords[2]-x0))  # heights at the ends of the tangent line
#     
#    lines(x.coords, y, col=add.alpha("red",.6), lwd=4)
#    points(x0,y0,cex=1.5)
# }
# 
# ## Let's test it
# x <- 9
# n <- 10
# t <- c(.40,.60)#,.75,.85, .95) #test values of theta
# 
# par(mfrow=c(2,3))
# 
# #probability scale cdf
# plot(theta,pbeta(theta,x+1,n-x+1),type="n",ylab="probability",main="Theta cdf")
# for(i in 1:length(t)){
#     tangentLine(t[i], .075, x=x, n=n, param = "theta")
# }
# lines(theta,pbeta(theta,x+1,n-x+1))
# 
# #log-odds scale cdf
# plot(logit(theta),pbeta(theta,x+1,n-x+1),type="n",ylab="probability",main="Logit(theta) cdf")
# for(i in 1:length(t)){
#     tangentLine(logit(t[i]), .8, x=x, n=n, param = "lodds")
# }
# lines(logit(theta),pbeta(theta,x+1,n-x+1),type="l")
# 
# #odds scale cdf
# plot(odds(theta),pbeta(theta,x+1,n-x+1),type="n",ylab="probability",main="Odds(theta) cdf")
# for(i in 1:length(t)){
#     tangentLine(odds(t[i]), .1, x=x, n=n, param = "odds")
# }
# lines(odds(theta),pbeta(theta,x+1,n-x+1),type="l")
# #This one is giving me trouble for the tangent lines due to the changing scale!!
# 
# # prob scale pdf
# plot(theta,posteriors[["theta"]](theta,x,n),type="l",ylab="density", main="Theta pdf")
# 
# # log-odds scale pdf
# plot(logit(theta),posteriors[["lodds"]](logit(theta),x,n),type="l",ylab="density", main="Logit(theta) pdf")
# 
# # odds scale pdf
# plot(odds(theta),posteriors[["odds"]](odds(theta),x,n),type="l",ylab="density", main="Odds(theta) pdf")
# 
