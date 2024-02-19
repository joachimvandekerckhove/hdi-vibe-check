rm(list=ls())
HDIofMCMC_flawed <- function(sampleVec,credMass=0.95){
	# Verbatim of Kruschke's 2011 function to
        # compute the HDI.	
	sortedPts <- sort(sampleVec)
	ciIdxInc <- floor(credMass*length(sortedPts))
	nCIs <- length(sortedPts)-ciIdxInc
	ciWidth <- rep(0,nCIs)
	for(i in 1:nCIs){
		ciWidth[i] <- sortedPts[i+ciIdxInc]-sortedPts[i]
	}
	HDImin <- sortedPts[which.min(ciWidth)]
	HDImax <- sortedPts[which.min(ciWidth)+ciIdxInc]
	HDIlim <- c(HDImin,HDImax)
	return(HDIlim)
}

HDIofMCMC_quantile <- function(sampleVec,credMass=0.95){
	# Quantile interval.
	alpha <- (1-credMass)/2
	HDIlim <- quantile(sampleVec,probs=c(alpha,1-alpha))
	return(HDIlim)
}

# Parameters of simulation
successes <- 9
failures <- 1
rope <- c(0.388,0.612)
set.seed(123)
post_samples <- rbeta(9000,1+successes,1+failures)

transformation <- function(samples){
	# Defines the transformation to be applied to the
	# original posterior samples.
	res <- log(samples/(1-samples))
	return(res)
}

scaled_transformation <- function(samples,min=-5,max=5){
	# Rescales the tranformation to 0 and 1, just for
	# plotting purposes.
	original_transformation <- transformation(samples)
	res <- (original_transformation-min)/(max-min)+0
	return(res)
}

# Transform original posterior samples
transformed_samples <- scaled_transformation(post_samples)

scaled_density <- function(density,k=5){
	# Scales densities just for plotting purposes.
	new_density <- 1+(k-1)*density/max(density)
}

linked_segments <- function(original,color='#000000',...){
	# Plots original and transformed values connected by lines.
	segments(x0=original,y0=rep(1,length(original)),
 		x1=scaled_transformation(original),y1=rep(-1,length(original)),
		col=color,...)
}

centered_polygon <- function(center=c(0,0),width=1,height=1,...){
	# Draws a polygon of certain 'width' and 'height' centered
	# at 'center'.
	x <- c(rep(center[1]+width/2,2),rep(center[1]-width/2,2))
	y <- c(center[2]-height/2,rep(center[2]+height/2,2),center[2]-height/2)
	polygon(x,y,...)
}


# Figure begins
pdf(file='figure3.pdf',width=3.5,height=4)
#pdf(file='transformation.pdf',width=3,height=4)
par(mai=rep(0,4),cex.axis=0.7,tck=-0.02)
k <- 5 # ylim in both directions
plot(NULL,xlim=c(-0.05,1.05),ylim=c(-k,k),axes=F,ann=F)
	# 1) ROPE polygon
	rope_st <- scaled_transformation(rope)
	polygon(x=c(rope_st,rope_st[2],rope[2],rope[2],
		    rope[1],rope[1],rope_st[1]),
		y=c(rep(-k,2),-1,1,rep(k,2),1,-1),
		border=F,
		col='#cccccc')
	# 2) Connecting lines
	equally_spaced <- seq(0.05,.95,length.out=21)
	linked_segments(equally_spaced,col='#dddddd')
	# 3) Axes
	# 3a) White boxes around the axes labels
	at <- seq(0,1,0.25)
	for(a in 1:length(at)){
		centered_polygon(c(at[a],.57),.075,.35,col='#ffffffdd',border=F)
		centered_polygon(c(at[a],-.57),.065,.35,col='#ffffffdd',border=F)
	}
	# 4) Panel labels
	text(.2,2,expression(theta))
	text(.2,-2,expression(paste('logit(',theta,')')))
	# 5) Density curves
	ht_original <- hist(post_samples,breaks=seq(0,1,length.out=150),plot=F)
	tr_samples <- transformed_samples[transformed_samples>=0&transformed_samples<=1]
	ht_transformed <- hist(tr_samples,breaks=seq(0,1,length.out=150),plot=F)
	lines(ht_original$mids,scaled_density(ht_original$density),col='#555555')
	lines(ht_transformed$mids,-scaled_density(ht_transformed$density),col='#555555')
	# 5) HDIs
lines(c(HDIofMCMC_quantile(post_samples)[1],HDIofMCMC_quantile(transformed_samples)[1]),c(1,-1),lwd=1.75)
lines(c(HDIofMCMC_quantile(post_samples)[2],HDIofMCMC_quantile(transformed_samples)[2]),c(1,-1),lwd=1.75)
lines(c(HDIofMCMC_flawed(post_samples)[1],HDIofMCMC_flawed(transformed_samples)[1]),c(1,-1),lwd=1.75,lty='21')
lines(c(HDIofMCMC_flawed(post_samples)[2],HDIofMCMC_flawed(transformed_samples)[2]),c(1,-1),lwd=1.75,lty='21')
	# 3b) Axes and axes labels
	axis(1,pos=1,padj=-2,at=seq(0,1,0.25),
		labels=paste(seq(0,1,0.25)))
	axis(1,pos=-1,padj=-6,tck=0.02,
	     at=seq(0,1,0.25),
	     labels=round(transformation(seq(0,1,0.25)),2))
	# 6) Equally-spaced points
	points(equally_spaced,rep(1,length(equally_spaced)),pch=21,cex=0.65,bg='#ffffff',col='#dddddd')
	points(scaled_transformation(equally_spaced),rep(-1,length(equally_spaced)),pch=21,cex=0.65,bg='#ffffff',col='#dddddd')
	# 7) Legend
	legend(0,-4.5,yjust=0.5,xjust=0,lwd=2,cex=.65,
	       seg.len=1.5,x.intersp=0.5,pt.cex=1.25,
	       pch=c(15,NA,NA),col=c('#cccccc','#000000','#000000'),
	       lty=c(NA,'solid','21'),
		legend=c('ROPE','Quantile Interval','HDI'),box.lty='blank')
dev.off()
