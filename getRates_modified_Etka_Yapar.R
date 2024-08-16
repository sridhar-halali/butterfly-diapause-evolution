getrates_ey <- function(tree, bins, maps, type, log=F, uncert="confint") { 
  
  isNumberInRange <- function(number, rng){
    if(number <= rng[2] & number > rng[1]){
      return(TRUE)
    }
    return(FALSE)
  }
  
  
  if (type == "all"){
    changes<-countchanges(tree, maps)
  } else {
    changes<-countchanges(tree, maps)
    changes<-parsetypes(changes, type)
  }
  
	h<-max(nodeHeights(tree))
	b<-bins
	segs<-cbind(seq(0,h-h/b,h/b), seq(1/b*h,h,h/b))
	
	#compute the mean number of sampled changes for each segment
  subset_segments <- apply(segs, 1, function(rng){
    changes_this_range_this_type <- sapply(changes, function(sm){
      if(nrow(sm) == 0){
        return(0)
      }
      inThisRange <- sapply(sm[,1], function(number){
        isNumberInRange(number, rng)
      })
      result = sum(inThisRange)
      return(result)
    })
    N <- length(maps)
    avg <- sum(changes_this_range_this_type) / N
    if(uncert == "confint"){
      se  <- sd(changes_this_range_this_type) / sqrt(N)
      upper <- avg + 2*se
      lower <- avg - 2*se
      return(c(avg, lower, upper))
    }
    if(uncert == "quantile"){
      return(c(avg, quantile(changes_this_range_this_type, c(0.25, 0.75))))
    }
})
  
  subset_segments <- t(subset_segments)
  colnames(subset_segments) <- c("mean", "lwr", "upr")
  nchanges <- subset_segments[,"mean"]
  lchanges <- subset_segments[,"lwr"]
  uchanges <- subset_segments[,"upr"]
	#print("accounting for LTT")
	#control for the total edge length present in a segment using LTT computation
	LTT<-ltt(tree,plot=FALSE)
	
	LTT<-cbind(LTT$ltt[2:(length(LTT$ltt)-1)],
    LTT$times[2:(length(LTT$ltt)-1)],
    LTT$times[3:length(LTT$ltt)])
	
	ii<-1
	edge.length<-rep(0,b)
	
	#print(segs)
	for(i in 1:nrow(segs)){
		done.seg<-FALSE
		while(LTT[ii,2]<=segs[i,2]&&done.seg==FALSE){
			edge.length[i]<-edge.length[i]+
			LTT[ii,1]*(min(segs[i,2],LTT[ii,3])-
			max(segs[i,1],LTT[ii,2]))
		if(LTT[ii,3]>=segs[i,2]) done.seg<-TRUE
		if(LTT[ii,3]<=segs[i,2]) ii<-if(ii<nrow(LTT)) ii+1 else ii
    	}
	}
	
	#print bins thru time
	timesegs<-h-as.vector(t(segs))
	if(log==T) {
	  magnitude<-rbind(nchanges/log(edge.length),nchanges/log(edge.length))
	  lmagnitude<-rbind(lchanges/log(edge.length),lchanges/log(edge.length))
	  umagnitude<-rbind(uchanges/log(edge.length),uchanges/log(edge.length))
	} else {
	  magnitude<-rbind(nchanges/edge.length,nchanges/edge.length)
	  lmagnitude<-rbind(lchanges/edge.length,lchanges/edge.length)
	  umagnitude<-rbind(uchanges/edge.length,uchanges/edge.length)
	}
	
	
	
	out<-list()
	out[[1]]<-segs
  out[[2]]<-timesegs
  out[[3]]<-magnitude
	out[[4]]<-type
	out[[5]]<-lmagnitude
	out[[6]]<-umagnitude
  names(out)<-c("segs", "timesegs", "magnitude", "type", "mag.lwr", "mag.upr")
  
  return(out)

}

rateplot_ey <- function(rates, tree, ylim, spline,color, width, lty, alpha, confint=F){
  if (spline == T) {
    smoothed <- smooth.spline(x=jitter(rates$timesegs, factor=0.1), y=rates$magnitude, cv=T)
    newx<- seq(from = max(rates$timesegs), to=0, length.out = 100)
    pred<-predict(smoothed, x = newx)
    #print(c(max(rates$segs), min(rates$segs)))
    plot(x=newx, y=pred$y, 
         lwd=width, type='l', 
         xlim=c(max(rates$segs), min(rates$segs)), 
         ylab="mean number of changes / total edge length", 
         xlab="time since the present", 
         #main=rates$type, 
         ylim,
         col=color, 
         lty=lty)
    if(confint ==T){
    midtimepoints <- apply(rates$segs, 1, function(x){
      return(x[1] + ((x[2] - x[1])/2) )
    })
    midtimepoints <- rev(midtimepoints)
    midtimepoints <- jitter(midtimepoints, factor=0.2)
    points(midtimepoints, rates$magnitude[1,], col=color)
    segments(midtimepoints, rates$mag.lwr[1,], midtimepoints, rates$mag.upr[1,], col=color)
    }
    
  } else {
    
    plot(rates$timesegs, rates$magnitude, 
         lwd=width, type="l", 
         xlim=c(max(rates$segs), min(rates$segs)),
         lend=0,xlab="time since the present",
         ylab="mean number of changes / total edge length", 
         #main=rates$type, 
         ylim)

  }
  
  plotTree(tree,add=TRUE,ftype="off",lwd=1,color=make.transparent("grey", alpha),
           mar=par()$mar,direction="leftwards",xlim=c(max(rates$segs),min(rates$segs)))
  
  #abline(v=66, col="red", lwd=2)
}

segplot_ey <-function(rates, ylim, color, spline, width, ax=T, lty,
                  y.lab="mean number of changes / total edge length", x.lab="time since the present", confint=F) {
  if (spline == T) {
    smoothed <- smooth.spline(x=jitter(rates$timesegs, factor=0.1), y=rates$magnitude, cv=T)
    newx<- seq(from = max(rates$timesegs), to=0, length.out = 100)
    pred<-predict(smoothed, x = newx)
    plot(x=newx, y=pred$y, 
         lwd=width, type='l', 
         xlim=c(max(rates$segs), min(rates$segs)), 
         ylab=y.lab, 
         xlab=x.lab, 
         #main=rates$type, 
         ylim,
         col=color,
         axes=ax,
         lty=lty)
    if(confint == T){
    midtimepoints <- apply(rates$segs, 1, function(x){
      return(x[1] + ((x[2] - x[1])/2) )
    })
    midtimepoints <- rev(midtimepoints)
    midtimepoints <- jitter(midtimepoints, factor=0.2)
    points(midtimepoints, rates$magnitude[1,], col=color)
    segments(midtimepoints, rates$mag.lwr[1,], midtimepoints, rates$mag.upr[1,], col=color)
    }
    
  } else {
    
    plot(rates$timesegs, rates$magnitude, 
         lwd=width, type="l", 
         xlim=c(max(rates$segs), min(rates$segs)),
         lend=0,xlab="time since the present",
         ylab="mean number of changes / total edge length", 
         #main=rates$type, 
         ylim,
         col=color,
         axes=ax)
    
  }
}
