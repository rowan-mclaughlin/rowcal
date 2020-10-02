#   r o w C a l   

#   Utilities for radiocarbon age calibration in R
#   by T. Rowan McLaughlin email : r.mclaughlin@qub.ac.uk

#   Please cite McLaughlin, T. R. 2019. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory 26: 479-501. https://doi.org/10.1007/s10816-018-9381-3

intcal<-read.csv(url('http://intcal.org/curves/intcal20.14c'),header=FALSE,skip=11)
marine<-read.csv(url('http://intcal.org/curves/intcal20.14c'),header=FALSE,skip=11)
colnames(intcal)<-c("calBP","14Cage","Error","Delta14C","Sigma")
colnames(marine)<-c("calBP","14Cage","Error","Delta14C","Sigma")

rowcal<-function(date, sigma, calcurve=intcal, BC=TRUE, norm=TRUE) {
    curve<-calcurve[(calcurve[,2]<= (date+(sigma*4)) & calcurve[,2]>= (date-(sigma*4)) ),]
	
	#convert from cal. BP to cal. BC if required
    if(BC==TRUE) curve$calBP<-1950-curve$calBP
    
    #build an output table
    calibrated<-as.matrix(curve[,1:2]) 
    if(BC==TRUE) colnames(calibrated)<-c("Cal. BC/AD","Relative p") else colnames(calibrated)<-c("Cal. BP","Relative p")
    
    #apply the normal distribution function to each margin of the calibration curve and average
    H<-dnorm(curve[,2]+curve[,3],mean=date,sd=sigma)
    L<-dnorm(curve[,2]-curve[,3],mean=date,sd=sigma)
    M<-(H+L)/2
    Res<-abs(mean(diff(curve[,1])))
    if(norm==TRUE) calibrated[,2]<-(M/sum(M))/Res else calibrated[,2]<-M/Res
    calibrated
}

#   Function for mixing curves, or applying local delta-R corrections, or both

mixcurves<-function(mix=0.5, uncert=0, c1=intcal, c2=marine, delR=0, error_delR=0) {
	# First interpolate so that the two curves have the same calendar years
	c2y<-approx(c2[,1], c2[,2], c1[,1], rule=2)$y
	c2e<-approx(c2[,1], c2[,3], c1[,1], rule=2)$y
	
	# Apply delta-R if so required
	c2y<-c2y+delR
	c2e<-sqrt(c2e^2 + error_delR^2) 
	
	# Perform curve mixing 
	out<-c1[,1:3] 
    out[,2] <- (1-mix)*c1[,2] + mix * c2y 
    out[,3] <- sqrt(((1-mix)*c1[,3])^2 + (mix*c2e)^2 + (uncert * (c1[,3]-c2e))^2)
    return(out)
}

#____________________________________________________________________________________
#
#function to find modal calibrated date for a C14 date, based on rowcal above
rowcalmode<-function(date,sigma,...) {
 cd<-rowcal(date,sigma,...)
 cd[,1][cd[,2]==max(cd[,2])]
}

#____________________________________________________________________________________
#
#function to find weighted mean date for a C14 date, based on rowcal above
rowcalwm<-function(date,sigma,...) {
   g<-rowcal(date,sigma,...)
   sum(g[,1]*g[,2]/max(g[,2]))/sum(g[,2]/max(g[,2]))
}

#____________________________________________________________________________________
#
#function to find (even just one) age sample for C14 date, based on rowcal above
rowcalsam<-function(date, sigma, N=1, ...) {
  g<-rowcal(date, sigma, ...)
  random.points <- approx(cumsum(g[,2])/sum(g[,2]),g[,1],runif(N))$y
  random.points
}

#____________________________________________________________________________________
#
#function to find median of a C14 date
rowcalmedian<-function(date, sigma, ...) {
   d<-rowcal(date,sigma,...)
   return(d[which.min(abs(cumsum(d[,2])-0.1)),1])
}

#function to make, from the gamma distribution, sapwood estimates or anything else with a gamma distro
sap<-function(date, sigma, shape=sigma, scale=3, N=1, ...) {
	g<-matrix(nrow=201,ncol=2)
	g[,1]<-seq(date-sigma,date-sigma+100,0.5)
	g[,2]<-dgamma(c(1:201),shape=shape, scale=scale, ...)
	g
}


#function to sample from sapwood estimates, like rowcalsam
sapsam<-function(date, sigma, N=1, ...) {
	g<-sap(date, sigma, ...)
  	random.points <- approx(cumsum(g[,2])/sum(g[,2]),g[,1],runif(N))$y
  	random.points
}

#function to make a normal distribution of calendar dates
calen<-function(date, sigma, sigmas=2, res=(if(sigma<5) 0.1 else 1)) {
	# add uncertainty to zero sigmas
  if(sigma==0) sigma<-sigma+0.5
  x<-seq(date-sigmas*sigma,date+sigmas*sigma,res)
	g<-matrix(nrow=length(x),ncol=2)
	g[,1]<-x
	g[,2]<-dnorm(x,mean=date,sd=sigma)
	g
}
	
#function to sample from calendar dates, like rowcalsam
calsam<-function(date, sigma, N=1, ...) {
	g<-calen(date, sigma, ...)
	random.points <- approx(cumsum(g[,2])/sum(g[,2]),g[,1],runif(N))$y
  	random.points
}

#function to turn a three-column list of dates into a list of one sample per date
#third column must be the name of the calibration curve, or 'cal' meaning calendar, 'sap' for sapwood
MCmix<-function(L) {
	colnames(L)<-c('date','error','cc')
	spp<-L[L$cc=='sap',]; cal<-L[L$cc=='cal',]
	cal_points<-c(); sap_points<-c()
	if(nrow(spp)>0) sap_points<-apply(spp[,c('date','error')],1,function(rdate,...) sapsam(rdate['date'],rdate['error'],...) )
	if(nrow(cal)>0) cal_points<-apply(cal[,c('date','error')],1,function(rdate,...) calsam(rdate['date'],rdate['error'],...) )
    C14_points<-c()
    C14<-L[!(L$cc %in% c('sap','cal')),]
    if(nrow(C14)>0) {
    # I know this is as ugly as sin, but it works
	for (N in 1:length(C14[,1])) 
		C14_points[N]<-eval(parse(text=paste("rowcalsam(",C14[N,1],",",C14[N,2],",calcurve=",as.character(C14[N,3]),")",sep='')))
	}
	return(as.vector(c(cal_points,sap_points,C14_points)))
}	
	
#function to turn a two-column list of uncalibrated dates into list of one sample per date
MCsam<-function(L,...) {
   colnames(L)<-c('BP','SD')
   O<-apply(L[,c('BP','SD')],1,function(rdate,...) rowcalsam(rdate['BP'],rdate['SD'],...) )
   O
}
   
#function to make best guess list from same
findwm<-function(L,...) {
   colnames(L)<-c('BP','SD')
   O<-apply(L[,c('BP','SD')],1,function(rdate,...) rowcalwm(rdate['BP'],rdate['SD'],...) )
   O
}
  

# function to make mode-based version of same
findmode<-function(L,...) {
   colnames(L)<-c('BP','SD')
   O<-apply(L[,c('BP','SD')],1,function(rdate,...) rowcalmode(rdate['BP'],rdate['SD'],...) )
   O
}
  

# function for reverse calibrating
revcal<-function(N, calcurve=intcal, BC=TRUE) {
 out<-NA
 if (BC) {
   calcurve$BC<-1950-calcurve$calBP
   out<-calcurve[which(abs(calcurve$BC-N)==min(abs(calcurve$BC-N))),2]	  
   }
   if (!BC) out<-calcurve[which(abs(calcurve$calBP-N)==min(abs(calcurve$calBP-N))),2]
 out

}

# function for inputing data from clipboard independent of platform
CLIP<-function() {
	pasted<-character()
	if (Sys.info()[1]=='Darwin') pasted<-read.table(pipe('pbpaste'),head=FALSE) else pasted<-read.delim('clipboard',head=FALSE)
	pasted
}

# Function for calcualting Gaussian KDE of set of dates using a Monte Carlo method 
# Plots an entertaining build-up depicting each Monte Carlo run
MCdensity2<-function(DATA=CLIP(),N=500,col=rgb(0,0,0,0.01),add=FALSE,bw=30,...) {
	if (!add) plot(density(MCsam(DATA),bw,na.rm=TRUE),col=NA,...)
    for (P in 1:N) lines(density(MCsam(DATA),bw,na.rm=TRUE),col=col)
}

# Function for calculating Gaussian KDE of set of dates using a Monte Carlo method, with output data allowing summary stats 
MCdensity<-function(DATA=CLIP(),N=100,perm_runs=1, perm_size=1, col=rgb(0,0,0,0.01),plot.new=FALSE,add=FALSE,bw=30,...) {
	if (perm_runs==1 & perm_size<1) stop("perm_runs must be specified and >1 if permutations of the data are required")
	if (perm_runs>1 & perm_size<1) {
		oDATA<-DATA
		DATA<-oDATA[sample(1:nrow(oDATA),round(perm_size*nrow(oDATA))),] }
	#by default the densities are calculated at 512 points in time n=512 below
	#this can be changed but needs to be a power of 2 --- 1024 2048 etc etc (must be a Fourier transform at work)
	d<-density(MCsam(DATA),bw,na.rm=TRUE,n=512)
	if (plot.new) {add<-TRUE; plot(d,col=col,...)}
	x1<-min(round(d$x)); x2<-max(round(d$x)); n=d$n
    out<-matrix(nrow = 512, ncol = N*perm_runs+1); out[,1]<-d$x; out[,2]<-d$y; P<-0
	pb <- txtProgressBar(min=2,max=N*perm_runs-1,initial=2)
	for (perm in 1:perm_runs) {
	  for (run in 2:N-1) {
	   	  P<-P+1
		  if (perm_runs>1 & perm_size<1) DATA<-oDATA[sample(1:nrow(oDATA),round(perm_size*nrow(oDATA))),]	
	      d<-density(MCsam(DATA),bw,na.rm=TRUE,from=x1,to=x2,n=512)
		  out[,P+2]<-d$y 
		  if (add & P>2) lines(out[,1],rowMeans(out[,2:P]),col=col)
		  setTxtProgressBar(pb,P+1)
	  }
	}
	close(pb)
	class(out)='MCd'; return(out)
}



# Function for calculating Gaussian KDE of set of dates using a Monte Carlo method, with output data allowing summary stats
# differs from MCdensity in that there should be three input columns, specifying date, error, and curve
# curve can be intcal, or any intcal-like object (marine)
# alternatively, curve can be 'sap' for sapwood models, or 'cal' for calendar years
mixdensity<-function(DATA=CLIP(),N=100,perm_runs=1, perm_size=1, col=rgb(0,0,0,0.01),plot.new=FALSE,add=FALSE,bw=30,...) {
	if (perm_runs==1 & perm_size<1) stop("perm_runs must be specified and >1 if permutations of the data are required")
	if (perm_runs>1 & perm_size<1) {
		oDATA<-DATA
		DATA<-oDATA[sample(1:nrow(oDATA),round(perm_size*nrow(oDATA))),] }
	#by default the densities are calculated at 512 points in time n=512 below
	#this can be changed but needs to be a power of 2 --- 1024 2048 etc etc (must be a Fourier transform at work)
	d<-density(MCmix(DATA),bw,na.rm=TRUE,n=512)
	if (plot.new) {add<-TRUE; plot(d,col=col,...)}
	x1<-min(round(d$x)); x2<-max(round(d$x)); n=d$n
    out<-matrix(nrow = 512, ncol = N*perm_runs+1); out[,1]<-d$x; out[,2]<-d$y; P<-0
	pb <- txtProgressBar(min=2,max=N*perm_runs-1,initial=2)
	for (perm in 1:perm_runs) {
	  for (run in 2:N-1) {
	   	  P<-P+1
		  if (perm_runs>1 & perm_size<1) DATA<-oDATA[sample(1:nrow(oDATA),round(perm_size*nrow(oDATA))),]	
	      d<-density(MCmix(DATA),bw,na.rm=TRUE,from=x1,to=x2,n=512)
		  out[,P+2]<-d$y 
		  if (add & P>2) lines(out[,1],rowMeans(out[,2:P]),col=col)
		  setTxtProgressBar(pb,P+1)
	  }
	}
	close(pb)
	class(out)='MCd'; return(out)
}



# --- Identify phases using hierarchical clustering ---
# (similar to the process implemented by Bevan and Crema 2018)
# returns a random selection of one site per 'bin' whose temporal width is defined by the parameter h (defaults to 30 years)
# if parameter shuffle is TRUE (which it is by default) then a random selection of date-per-phase is selected
# otherwise the first instance of a date in the datelist is selected

phasesam<-function(DATA=CLIP(), h=30, method='median', shuffle=TRUE) {
   if(!(ncol(DATA) %in% c(3,4))) stop('Input should have 3 or 4 columns (siteid, date, error, [calcurve])')
   if(!(method %in% c('median','sample'))) stop('`method` should be `median` to bin using median dates or `sample` for MC draws')
   # make a copy of the input so that the columns can be manipulated and shuffled if required
   datelist<-DATA
   if (shuffle) { 
         datelist<-datelist[sample(nrow(datelist)),]
         rownames(datelist)<-rownames(DATA) }
   sites<-unique(datelist[,1])
   if(ncol(DATA)==3) datelist$calcurve<-'intcal'
   colnames(datelist)<-c('site','bp','sd','calcurve')
   if(method=='median') datelist$pointest<-findmixmedian(datelist[,c(2:4)])
   if(method=='sample') datelist$pointest<-MCmix(datelist[,c(2:4)])
   # The below repeats the MC sampling in case NAs are somehow introduced
   # This is actually a bug in MCmix: note to self to fix 						****
   count<-0
   while(anyNA(datelist$pointest)) {
        count<-count+1
        if(count==10) stop('cannot calibrate dates')
        datelist$pointest<-MCmix(datelist[,c(2:4)])
   }
   datelist$timebin <- 1
   for(N in 1:length(sites)) {
   	  sitedates<-datelist[datelist[,1]==sites[N],'pointest']
   	  if(length(sitedates)>1) {
   	       sitedateindex<-which(datelist[,1]==sites[N])
           dendro<-hclust(dist(sitedates))
   	       sitebins<-cutree(dendro, h=h)
           datelist[sitedateindex,'timebin'] <- sitebins
        }
    }
    output<-datelist[!duplicated(datelist[,c('site','timebin')]),]
    return(output)  
}

phasedensity<-function(DATA=CLIP(),N=100,perm_runs=1, perm_size=1, col=rgb(0,0,0,0.01),plot.new=FALSE,add=FALSE,bw=30,h=30,method='sample',shuffle=TRUE,...) {
  if(!is.data.frame(DATA)) DATA<-as.data.frame(DATA)
  if(ncol(DATA)==3) mix<-FALSE
  if(ncol(DATA)==4) mix<-TRUE
  if(!(ncol(DATA) %in% c(3,4))) stop('Input should have 3 or 4 columns (siteid, date, error, [calcurve])')
  if(perm_runs==1 & perm_size<1) stop("perm_runs must be specified and >1 if permutations of the data are required")
  if (perm_runs>1 & perm_size<1) {
    oDATA<-DATA
    DATA<-oDATA[sample(1:nrow(oDATA),round(perm_size*nrow(oDATA))),] }  
  #by default the densities are calculated at 512 points in time; n=512 below
  #this can be changed but needs to be a power of 2 --- 1024 2048 etc etc
  if(mix)  d<-density(MCmix(phasesam(DATA, h=h)[,2:4]),n=512)
  if(!mix) d<-density(MCsam(phasesam(DATA, h=h)[,2:3]),n=512)
  if (plot.new) {add<-TRUE; plot(d,col=col,...)}
  x1<-min(round(d$x)); x2<-max(round(d$x)); n=d$n
  out<-matrix(nrow = 512, ncol = N*perm_runs+1); out[,1]<-d$x; out[,2]<-d$y; P<-0
  pb <- txtProgressBar(min=2,max=N*perm_runs-1,initial=2)
  for (perm in 1:perm_runs) {
      for (run in 2:N-1) {
      P<-P+1
      if (perm_runs>1 & perm_size<1) DATA<-oDATA[sample(1:nrow(oDATA),round(perm_size*nrow(oDATA))),]	
      # call the MC sampling engine (MCsam is slightly more efficient than MCmix, so is used if no calcurve is specified)
      if(mix)  d<-density(MCmix(phasesam(DATA, method=method, h=h)[,2:4]),bw,na.rm=TRUE,from=x1,to=x2,n=512)
      if(!mix) d<-density(MCsam(phasesam(DATA, method=method, h=h)[,2:3]),bw,na.rm=TRUE,from=x1,to=x2,n=512)
      out[,P+2]<-d$y 
      if (add & P>2) lines(out[,1],rowMeans(out[,2:P]),col=col)
      setTxtProgressBar(pb,P+1)
    }
  }
  close(pb)
  class(out)='MCd'; return(out)
}

# The following function works in the same way as phasesam, but returns the full
# input datelist, with an added column indicating the bin allocated to each date
# (could be easily modified to use MC draws instead of median dates -- see phasesam() )
phaser<-function(DATA=CLIP(), h=30){
   if(!(ncol(DATA) %in% c(3,4))) stop('Input should have 3 or 4 columns (siteid, date, error, [calcurve])')
   # make a copy of the input so that the columns can be manipulated
   datelist<-DATA
   sites<-unique(datelist[,1])
   if(ncol(DATA)==3) datelist$calcurve<-'intcal'
   datelist$median<-findmixmedian(datelist[,c(2:4)])
   datelist$timebin <- 1
   for(N in 1:length(sites)) {
   	  sitedates<-datelist[datelist[,1]==sites[N],'median']
   	  if(length(sitedates)>1) {
   	       sitedateindex<-which(datelist[,1]==sites[N])
           dendro<-hclust(dist(sitedates))
   	       sitebins<-cutree(dendro, h=h)
           datelist[sitedateindex,'timebin'] <- sitebins
        }
    }
    output<-datelist[,c(colnames(DATA),'timebin')]
    if(!quiet) 
    return(output)  
}

# Suggest bandwidth via R's default methods
suggest_bw<-function(DATA=CLIP(),calcurve='intcal') {
   if(!(ncol(DATA) %in% c(2,3))) stop('Input must be two or three column table specifing date, sd, [calcurve]')
   if(ncol(DATA)==2) DATA[,3]<-calcurve
   m<-c()
   for(N in 1:nrow(DATA)) m[N]<-rowcalmedian(DATA[N,1],DATA[N,2],calcurve=eval(parse(text=DATA[N,3])))
   bwselect<-c('nrd0','nrd','ucv','bcv','SJ-ste', 'SJ-dpi')
   out<-c()
   for(N in 1:6) try(out[N]<-density(m,bwselect[N],na.rm=TRUE,n=512)$bw)
   names(out)<-bwselect
   return(round(out,0))  
}

# Method for plotting density models
plot.MCd<-function(x, style=1, add=FALSE,col=rgb(0,0,0.8,0.5),fill=rgb(.5,.5,.5,.5),scalefactor=1,ylab='Density',xlab='Automatic',grid=TRUE,...) { 
  if(style==2) plotMCd2(x,add=FALSE,col=col,fill=fill,scalefactor=scalefactor,xlab=xlab,ylab=ylab,grid=grid,...) else {
  if(style!=1) warning('`style` must be 1 or 2. Defaulting to 1')
  M<-rowMeans(x[,2:ncol(x)],na.rm=TRUE)*scalefactor
  Sd<-apply(x[,2:ncol(x)],1,sd,na.rm=TRUE)*scalefactor
  if (!add) {
  	plot(x[,1],M,col=NA,xlab=NA,ylab=ylab,...)
  	if(grid) grid(lwd=0.5)
    if(!is.na(xlab)) { if(xlab=='Automatic') { 
      xlab<-'Cal. BC/AD'
      pw<-par('xaxp')
      if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
      if(pw[2]<0) xlab<-'Cal. BC'
  	}}
  	title(xlab=xlab)
  }
  lines(x[,1],M+Sd,lwd=1,col=col)
  lines(x[,1],M-Sd,lwd=1,col=col)
  polygon(c(x[,1], rev(x[,1])), c(M+Sd, rev(M-Sd)), col = col, border = NA)
  }
}  

# New style for plotting density models
plotMCd2<-function(x,add=FALSE,col=rgb(0,0,0.8,0.8),fill=rgb(.5,.5,.5,.5),scalefactor=1,ylab='Density',xlab='Automatic',grid=TRUE,...) { 
  M<-rowMeans(x[,2:ncol(x)],na.rm=TRUE)*scalefactor
  Sd<-apply(x[,2:ncol(x)],1,sd,na.rm=TRUE)*scalefactor
  if (!add) {
  	plot(x[,1],M,col=NA,xlab=NA,ylab=ylab,...)
  	if(grid) grid(lwd=0.75)
    if(!is.na(xlab)) { if(xlab=='Automatic') { 
      xlab<-'Cal. BC/AD'
      pw<-par('xaxp')
      if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
      if(pw[2]<0) xlab<-'Cal. BC'
  	}}
  	title(xlab=xlab)
  }
  lines(x[,1],M,lwd=2,col=col)
  lines(x[,1],M+Sd,lwd=1,lty=2,col=col)
  lines(x[,1],M-Sd,lwd=1,lty=2,col=col)
  polygon(c(x[,1], rev(x[,1])), c(M+Sd, rev(M-Sd)), col = fill, border = NA)
} 

# Method for adding an outline of density model to current plot
# e.g. lines(MyDensityModel)
lines.MCd<-function(x,scalefactor=1, ...) {
  M<-rowMeans(x[,2:ncol(x)],na.rm=TRUE)*scalefactor
  Sd<-apply(x[,2:ncol(x)],1,sd,na.rm=TRUE)*scalefactor
  lines(x[,1],M+Sd, ...)
  lines(x[,1],M-Sd, ...)
}  

# Method for calculating median of a density model 
median.MCd<-function(x) {
  out<-matrix(nrow=nrow(x),ncol=2)
  out[,1]<-x[,1]
  out[,2]<-apply(x[,2:ncol(x)],1,median,na.rm=TRUE)
  return(out)
}  


#To test for convergenge:

#test<-MCdensity(SOME_DATA,N=100)
#conv<-c(); nc<-ncol(test)-1
#T<-rowMeans(test[,2:nc+1],na.rm=T)
#for (N in 4:nc-1) conv[N]<-sum(T,na.rm=T)-sum(rowMeans(test[,2:N],na.rm=T))
#plot(conv,pch=1,col=rgb(0.7,0.1,0,0.5),xlab="Run number", ylab="Convergence")
#abline(h=0,col='grey')

MCsimulate<-function(mod,Nd=1000,type='KDE',...) {
   if(!type %in% c('KDE','SPD')) stop('`type` must be `KDE` or `SPD``')
   S<-approx(cumsum(mod$y)/sum(mod$y),mod$x,runif(Nd))$y
   # by default 1000 random simulated dates drawn from 'mod'
   S<-S[!is.na(S)]
   datelist<-matrix(nrow = length(S),ncol=2)
   datelist[,1]<-sapply(S,'revcal')
   datelist[,2]<-sample(25:40,length(S),replace=TRUE)
   if(type=='KDE') {
        out<-MCdensity(datelist,...)
        class(out)='MCd'
   }
   if(type=='SPD') {
        spd<-rowcalsum(datelist)
        out<-SPDboot(spd)
        class(out)<-'SPDboot'
   }
   return(out)
}



#e.g. try the following to add an exponential population model with annual increase of 0.5%
#the plot will start from a point clicked upon with the mouse
#lines(pop_sim(rate=0.5))
	
#function to make a MCdensity-like object from a list of calender dates without 14C uncertainty
#using bootstrapping (sampling with replacement)
bootstrap_density<-function(DATA=CLIP()[,1],N=100, col=rgb(0,0,0,0.01),bw=30,...) {
	#by default the densities are calculated at 512 points in time n=512 below
	#this can be changed but needs to be a power of 2 --- 1024 2048 etc etc
	d<-density(DATA,bw,na.rm=TRUE,n=512)
	x1<-min(round(d$x)); x2<-max(round(d$x)); n=d$n
    out<-matrix(nrow = 512, ncol = N+1); out[,1]<-d$x; out[,2]<-d$y; P<-0
	for (run in 2:N-1) {
		P<-P+1	
	    d<-density(sample(DATA,replace=TRUE),bw,na.rm=TRUE,from=x1,to=x2,n=512)
		out[,P+2]<-d$y 
	}
	class(out)='MCd'; return(out)
}

# Function to align two density models -- the values of 'A' will be interpolated to the timestamps of 'B'
align_densities<-function(A, B) {
    out<-B[,1:ncol(A)]
    for(N in 2:ncol(A)) out[,N]<-approx(x=A[,1],y=A[,N],xout=B[,1])$y
    class(out)<-'MCd'
    return(out)
}


# function to calculate the ratio of two KDE models
# Not yet tested for KDEs of different resolutions

`/.MCd`<-function(A, B) {
    # First extrapolate the shorter (in time) KDE to the longer one
    spanA<-diff(range(A[,1]))
    spanB<-diff(range(B[,1]))
    if( spanA > spanB ) B<-align_densities(B,A) else A<-align_densities(A,B)
    # Use the mean±sd as the envelope to divide by

     AM<-rowMeans(A[,2:ncol(A)],na.rm=TRUE)
     ASd<-apply(A[,2:ncol(A)],1,sd,na.rm=TRUE)
     BM<-rowMeans(B[,2:ncol(B)],na.rm=TRUE)
     BSd<-apply(B[,2:ncol(B)],1,sd,na.rm=TRUE)
     
     out<-matrix(nrow=nrow(A),ncol=3)
     out[,1]<-A[,1]
     out[,2]<-(AM+ASd)/(BM+BSd)
     out[,3]<-(AM-ASd)/(BM-BSd)
     class(out)<-'MCd_summary'
     return(out)
}

# Mean value of a density model, optionally between two points
mean.MCd<-function(A, from=NULL, to=NULL) { 
    if(is.null(from)) from<-min(A[,1])
    if(is.null(to)) to<-max(A[,1])
    return(mean(A[A[,1]>=from & A[,1]<=to, 2:ncol(A)], na.rm=TRUE))
}

# Standard deviation method for density models
# 
# 'sd' is not a generic function so need to hijack its default
# Note to self to remember to add the following when rowcal becomes a proper package:
# formals(sd.default) <- c(formals(sd.default), alist(... = ))
sd <- function(x, ...) UseMethod("sd")
sd.default <- stats::sd

sd.MCd<-function(A, from=NULL, to=NULL) { 
    if(is.null(from)) from<-min(A[,1])
    if(is.null(to)) to<-max(A[,1])
    return(sd(A[A[,1]>=from & A[,1]<=to, 2:ncol(A)], na.rm=TRUE))
}


# Method for plotting summary form of density models
plot.MCd_summary<-function(x,add=FALSE,col=rgb(0,0,0.8,0.5),scalefactor=1,xlab='Automatic',grid=TRUE,...) { 
  x<-x[!is.na(x[,2]),]
  x<-x[!is.na(x[,3]),]
  if (!add) {
  	plot(x[,1],x[,2],col=NA,xlab=NA,ylab='Density',...)
  	if(grid) grid(lwd=0.5)
    if(!is.na(xlab)) { if(xlab=='Automatic') { 
      xlab<-'Cal. BC/AD'
      pw<-par('xaxp')
      if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
      if(pw[2]<0) xlab<-'Cal. BC'
    }}  	
  	title(xlab=xlab)
  }
  lines(x[,1],x[,2],lwd=1,col=col)
  lines(x[,1],x[,3],lwd=1,col=col)
  polygon(c(x[,1], rev(x[,1])), c(x[,2], rev(x[,3])), col = col, border = NA)
}  


#Extract the various probablity densities for a given year from an object of class MCd
year_densities<-function(Y,MCd) MCd[which.min(abs(MCd[,1]-Y)),][-1]
	
#Wee functions for polygon plots, plot.C14 is the method of plot()-ing dates
fill<-function(j,col="grey",border=NA) polygon(c(j[1,1],j[,1],max(j[,1])), c(j[1,2],j[,2],j[1,2]),col=col,border=border)
plot.C14<-function(x,col="grey",...) {plot(x,col=NA,...) ; fill(x,col=col)}


# Functions for summarising radiocarbon dates

# Firstly, a function to return the highest density interval for a date 
# Heavily, er, 'influenced' by the code for same in A Parnell's Bchron !
 
hdr<-function(date, prob = 0.95) {
   date_a<-approx(date[,1],date[,2],c(min(date[,1]):max(date[,1])))
   ag = date_a$x
   de = date_a$y
    # Put the probabilities in order
    o = order(de)
    cu = cumsum(de[o])
    # Find which ones are above the threshold
    good_cu = which(cu>1-prob)
    good_ag = sort(ag[o][good_cu])
    # Pick out the extremes of each range
    breaks = diff(good_ag)>1
    where_breaks = which(diff(good_ag)>1)
    n_breaks = sum(breaks) + 1
    # Store output
    out = vector('list', length = n_breaks)
    low_seq = 1
    high_seq = ifelse(length(where_breaks)==0, length(breaks), where_breaks[1])
    for(i in 1:n_breaks) {
      out[[i]] = c(good_ag[low_seq], good_ag[high_seq])
      curr_dens = round(100*sum(de[o][seq(good_cu[low_seq], good_cu[high_seq])]),1)
      names(out)[[i]] = paste0(as.character(curr_dens),'%')
      low_seq = high_seq + 1
      high_seq = ifelse(i<n_breaks-1, where_breaks[i+1], length(breaks))
    }
    return(out)
}

# Function to return span of date; default is 95% CI

spansc14<-function(rdate,prob=0.95) {
	spans<-unlist(hdr(rdate,prob=prob))
	return(c(min(spans),max(spans)))
}

# Useful summary of C14 date

summary.C14<-function(rdate) {
	wm<-sum(rdate[,1]*rdate[,2]/max(rdate[,2]))/sum(rdate[,2]/max(rdate[,2]))
    spans<-unlist(hdr(rdate))
	out<-c(round(wm), round(min(spans)), round(max(spans)))
	names(out)<-c('wmean','l95','u95')
	return(out)
}

# Function to calibrate a whole datelist, adding a column with text 2-sig range

datelistcal<-function(dl=CLIP()) {
    dl$two_sigma<-NA
    for(N in 1:nrow(dl)) {
    	sm<-summary.C14(rowcal(dl[N,1],dl[N,2]))[2:3]
    	tx<-as.character(abs(sm))
    	if(sm[1]<=0 & sm[2]<=0) out<-paste(tx[1],'to',tx[2],'cal. BC')
    	if(sm[1]<=0 & sm[2]>0) out<-paste(tx[1],'cal. BC to cal. AD',tx[2])
    	if(sm[1]>0 & sm[2]>0) out<-paste('cal. AD',tx[1],'to',tx[2])
    	dl$two_sigma[N]<-out
    }
    return(dl)
}

# Wee functions to print the upper and lower bounds of a calibrated age

U95<-function(d,e) round(max(spansc14(rowcal(d,e))[1,], na.rm=TRUE),-1)
L95<-function(d,e) round(min(spansc14(rowcal(d,e))[1,], na.rm=TRUE),-1)
U68<-function(d,e) round(max(spansc14(rowcal(d,e))[2,], na.rm=TRUE),-1)
L68<-function(d,e) round(min(spansc14(rowcal(d,e))[2,], na.rm=TRUE),-1)

# 'Time matrix' function for summing / animated maps
# Returns a time matrix from input, which is a two-column dataframe
# The argument -- buffer -- describes how much time to wrap round the matrix 
# NOW DEPRECIATED -- use timematrix() instead

time_matrix<-function(datelist,buffer=100,BC=TRUE) {
  
  #set up a big matrix to hold results
  #first we need to "find" the range of the calibrated output
  
  if (ncol(datelist)<3) stop("time_matrix needs at least three columns in the input (ID, BP, error)")
  lower_date<-min(rowcalmode(datelist[(datelist[,2]==max(datelist[,2])),2],40,BC=BC))-buffer
  upper_date<-max(rowcalmode(datelist[(datelist[,2]==min(datelist[,2])),2],40,BC=BC))+buffer
  if (upper_date==-Inf) upper_date<-0
  date_range<-seq(lower_date,upper_date,by=5)
  
  #then make the matrix
    
  l<-length(datelist[,1]) 
  stmap<-as.data.frame(matrix(data=0,nrow=l,ncol=ncol(datelist)+(upper_date-lower_date)/5+1))
  colnames(stmap)<-c(colnames(datelist),as.character(date_range))  
  stmap[,1:ncol(datelist)]<-datelist                                               
  
  #perform the calculations, by default if a fourth column exists, it should contain the name of the calcurve to be used
  pb <- txtProgressBar(min=1,max=l,initial=1)
  for (d in 1:l) {
      cd<-rowcal(stmap[d,2],stmap[d,3],BC=BC)
      if (ncol(datelist)>3) {if (stmap[d,4]=='M') cd<-rowcal(stmap[d,2],stmap[d,3],BC=BC,calcurve=marine)}
	  tcd<-t(cd[(cd[,1] >=lower_date & cd[,1]<=upper_date),1:2])
      colnames(tcd)<-tcd[1,]
      tcd<-tcd[-1,]
      stmap[d,]<-replace(stmap[d,],values=tcd,list=names(tcd))
      setTxtProgressBar(pb,d)
    }
    close(pb)
	return(stmap)
}

#Version of the time_matrix function for samples with different calibration curves
#Custom calibration curves can be made with the mixcurves() function
#This is a good way of summing marine dates with custom deltaRs for example
# NOW DEPRECIATED -- use timematrix() instead

time_matrix_mix<-function(datelist,buffer=100,BC=TRUE) {
  
  #set up a big matrix to hold results
  #first we need to "find" the range of the calibrated output
  
  if (ncol(datelist)<3) stop("time_matrix needs at least four columns in the input (ID, BP, error, calcurve)")
  lower_date<-min(rowcalmode(datelist[(datelist[,2]==max(datelist[,2])),2],40,BC=BC))-buffer
  upper_date<-max(rowcalmode(datelist[(datelist[,2]==min(datelist[,2])),2],40,BC=BC))+buffer
  if (upper_date==-Inf) upper_date<-0
  date_range<-seq(lower_date,upper_date,by=5)
  
  #then make the matrix
    
  l<-length(datelist[,1]) 
  stmap<-as.data.frame(matrix(data=0,nrow=l,ncol=ncol(datelist)+(upper_date-lower_date)/5+1))
  colnames(stmap)<-c(colnames(datelist),as.character(date_range))  
  stmap[,1:ncol(datelist)]<-datelist                                               
  
  #perform the calculations, by default if a fourth column exists, it should contain the name of the calcurve to be used
  pb <- txtProgressBar(min=1,max=l,initial=1)
  for (d in 1:l) {
      cd<-eval(parse(text=paste("rowcal(",stmap[d,2],",",stmap[d,3],",calcurve=",as.character(stmap[d,4]),")",sep='')))
      tcd<-t(cd[(cd[,1] >=lower_date & cd[,1]<=upper_date),1:2])
      colnames(tcd)<-tcd[1,]
      tcd<-tcd[-1,]
      stmap[d,]<-replace(stmap[d,],values=tcd,list=names(tcd))
      setTxtProgressBar(pb,d)
    }
    close(pb)
	return(stmap)
}

timematrix<-function(datelist,buffer=100,BC=TRUE,res=10,default_curve='intcal',quiet=FALSE) {
  
  # make datalist a data frame is not already
  
  datelist<-as.data.frame(datelist)
  
  #first we need to "find" the range of the calibrated output
  
  if (ncol(datelist)<2) stop("timematrix needs at two or three columns in the input (BP, error, [calcurve])")
  if (ncol(datelist)==2) datelist[,3]<-default_curve 
  lower_date<-min(rowcalmode(datelist[(datelist[,1]==max(datelist[,1])),1],40,BC=BC))-buffer
  upper_date<-max(rowcalmode(datelist[(datelist[,1]==min(datelist[,1])),1],40,BC=BC))+buffer
  if (upper_date==-Inf) upper_date<-1950*BC
  # build a range of round numbers based on the required resolution 'res' 
  date_range<-seq(round(lower_date,-round(log10(res))),round(upper_date,-round(log10(res))),by=res)
  
  #then make the matrix
    
  l<-length(datelist[,1]) 
  timegrid<-matrix(data=0,ncol=l,nrow=length(date_range))
  
  #perform the calculations, by default if a third column exists, it should contain the name of the calcurve to be used
  if (!quiet) pb <- txtProgressBar(min=1,max=l,initial=1)
  for (d in 1:l) {
      cd<-eval(parse(text=paste("rowcal(",datelist[d,1],",",datelist[d,2],",calcurve=",as.character(datelist[d,3]),")",sep='')))
      timegrid[,d]<-approx(x=cd[,1], y=cd[,2], xout=date_range)$y
      if (!quiet) setTxtProgressBar(pb,d)
    }
    if (!quiet) close(pb)
	out<-list(datelist=datelist,calyears=date_range,timematrix=timegrid,res=res)
	class(out)<-'timematrix'
	return(out)
}


# code to sum and plot sums of C14 dates
# startcol can be specified to point to the first column of the input matrix that contains the radiocarbond data,
# by default this is the first numeric column name

# 'offset' can be specified to 'float' the sum into the yaxis some distance

# NOW DEPRECIATED -- use rowcalsum() instead

psum<-function(stt,add=FALSE, lines=TRUE, fill=TRUE,xlab="Cal. BC/AD",ylab="Summed prob. dens.",
				col.fill=rgb(0,0,0.8,0.5),col.lines=rgb(0,0,0.8,0.5), inv=FALSE, lty=1, lwd=1, 
				startcol=suppressWarnings(length(as.numeric(colnames(stt))[is.na(as.numeric(colnames(stt)))]))+1, offset=0,...) 
{  w<-length(stt[1,])
   h<-stt[,startcol:w]
   s<-colSums(h)
   x<-names(s)
   if (inv==TRUE) s<-s*-1
   if (add==FALSE) plot(x,s+offset,xlab=xlab,ylab=ylab,col=NA,...)
   if (lines==TRUE) lines(x,s+offset,col=col.lines,lty=lty,lwd=lwd)
   if (fill==TRUE) polygon(c(x[1],x,x[w-5]),c(offset,s+offset,offset),col=col.fill,border=NA)
   if (inv==TRUE) s<-s*-1 
   if (lines==FALSE && add==FALSE && fill==FALSE) return(s+offset)
}  

# User-friendly function to quickly sum any list of dates
# NOW DEPRECIATED -- use rowcalsum() instead
rowcalsum_old<-function(DATA=CLIP(),...) {mat<-time_matrix(cbind(NA,DATA)); psum(mat,...)}

#New way to sum dates that is much more flexible
rowcalsum<-function(DATA=CLIP(),norm=TRUE,quiet=FALSE) {
	nr<-ncol(DATA)
	mat<-timematrix(DATA, quiet=quiet)
	# Unpack the old-style time matrix (will optimise in future)
	S<-rowSums(mat$timematrix,na.rm=TRUE)
	if(norm==TRUE) S<-(S/nrow(mat$datelist))
	out<-list(
		dates=mat$datelist[,1],
		stds=mat$datelist[,2],
		curves=mat$datelist[,3],
		yrs=mat$calyears,		
		sum=S
	)
	class(out)='c14sum'
	return(out)
}

plot.c14sum<-function(S,add=FALSE,offset=0,xlab='Automatic',ylab='Summed prob.',type='f',col=rgb(0.8,0,0,0.8),lty=1,lwd=1,...) {
   if(!(type %in% c('f','l'))) stop("'type' must be 'f' for filled polygons or 'l' for lines")
   # bookend sum with 0s in case the plot is left 'hanging'
   x<-c( S$yrs[1], S$yrs, S$yrs[length(S$yrs)] )
   s<-c( 0, S$sum, 0)
   # replace NAs with 0s to avoid hanging plots
   s[is.na(s)]<-0
   if (add==FALSE) {
   	plot(x,s,xlab=NA,ylab=ylab,col=NA,...)
   	if(!is.na(xlab)) { if(xlab=='Automatic') { 
       xlab<-'Cal. BC/AD'
       pw<-par('xaxp')
       if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
       if(pw[2]<0) xlab<-'Cal. BC'
  	 }}
  	 title(xlab=xlab)
  	}
   if (type=='f') polygon(c(x[1],x,x[length(x)]),c(offset,s+offset,offset),col=col,border=NA)
   if (type=='l') lines(x,s+offset,col=col,lty=lty,lwd=lwd)
}

`+.c14sum`<-function(A,B) {
 	L=A; S=B
 	if(length(A$yrs)<length(B$yrs)) { L=B; S=A }
 	# First interpolate so that the two curves have the same calendar years
	add<-approx(x=S$yrs, y=S$sum, xout=L$yrs)
	add$y[is.na(add$y)]<-0
	out<-list(
		dates=c(S$dates, L$dates),
		stds=c(S$stds, L$stds),
		curves=c(S$stds, L$stds),
		yrs=add$x,		
		sum=add$y+L$sum
	)
	class(out)<-'c14sum'
	return(out)		
}

`-.c14sum`<-function(A,B) {
 	L=A; S=B; reverse=FALSE
 	if(length(A$yrs)<length(B$yrs)) { L=B; S=A; reverse=TRUE }
 	# First interpolate so that the two curves have the same calendar years
	add<-approx(x=S$yrs, y=S$sum, xout=L$yrs)
	add$y[is.na(add$y)]<-0
	if(reverse==TRUE) result<-add$y-L$sum else result<-L$sum-add$y
	out<-list(
		dates=c(S$dates, L$dates),
		stds=c(S$stds, L$stds),
		curves=c(S$stds, L$stds),
		yrs=add$x,		
		sum=result
	)
	class(out)<-'c14sum'
	return(out)		
}

`/.c14sum`<-function(A,B) {
 	L=A; S=B; reverse=FALSE
 	if(length(A$yrs)<length(B$yrs)) { L=B; S=A; reverse=TRUE }
 	# First interpolate so that the two curves have the same calendar years
	add<-approx(x=S$yrs, y=S$sum, xout=L$yrs)
	add$y[is.na(add$y)]<-0
	if(reverse==TRUE) result<-add$y/L$sum else result<-L$sum/add$y
	out<-list(
		dates=c(S$dates, L$dates),
		stds=c(S$stds, L$stds),
		curves=c(S$stds, L$stds),
		yrs=add$x,		
		sum=result
	)
	class(out)<-'c14sum'
	return(out)		
}


# Implementation of bootstrapping method for calculating SPD confidence intervals (Fernandez-Lopez de Pablo et al)
	
SPDboot<-function(SPD,Nboot=100,plot.new=FALSE,addlines=FALSE, col=rgb(0.8,0,0,0.1),...) { 
	g<-matrix(c(SPD$yrs,y=SPD$sum),ncol=2)
    stds<-c(quantile(SPD$stds)[2]:quantile(SPD$stds)[3])
	out<-matrix(0,nrow=length(SPD$yrs),ncol=Nboot)
	rownames(out)<-SPD$yrs	
	pb <- txtProgressBar(min=1,max=Nboot,initial=1)

	for(N in 1:Nboot) {
	   # pick X number of years from SPD where X is number of 14C dates 
	   sam<-approx(cumsum(g[,2])/sum(g[,2]),g[,1],runif(length(SPD$dates)),rule=2)$y
	   # reverse calibrate these years
		datelist<-matrix(nrow = length(sam),ncol=2)
        datelist[,1]<-sapply(sam,'revcal')
    	datelist[,2]<-sample(stds,length(sam),replace=TRUE)
       # make new SPD
        bootstrap<-rowcalsum(as.data.frame(datelist),quiet=TRUE)
        if(plot.new==TRUE & N==1) plot(bootstrap,col=NA,...)
        if(addlines==TRUE) plot(bootstrap,type='l',add=T,col=col)
        setTxtProgressBar(pb,N)
        out[,N]<-approx(x=bootstrap$yrs,y=bootstrap$sum,xout=SPD$yrs)$y
    }
    close(pb)
    class(out)<-'SPDboot'
    return(out)
}

lines.SPDboot<-function(SPDb,col=rgb(0,.8,0,.1),...) for(N in 1:ncol(SPDb)) lines(as.numeric(rownames(SPDb)),SPDb[,N],col=col,...)

# return CI envelope and mean
summary.SPDboot<-function(SPDb,probs=c(.05,.9)) {
	qu<-t(apply(SPDb,1,'quantile',probs=probs,na.rm=TRUE))
	class(qu)<-'SPDboot'
	return(qu)
}	


plot.SPDboot<-function(SPDb, probs=c(.16,.84), col=rgb(.1,.7,.1,.8), add=FALSE, fill=TRUE, xlab='Automatic',ylab='Summed prob.', ...) {
  sy<-summary(SPDb,probs=probs)
  x<-as.numeric(rownames(SPDb)); U<-sy[,2]; L<-sy[,1]
  if (add==FALSE) {
   	 plot(x,U,xlab=NA,ylab=ylab,col=NA,...)
   	  if(!is.na(xlab)) { if(xlab=='Automatic') { 
        xlab<-'Cal. BC/AD'
        pw<-par('xaxp')
        if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
        if(pw[2]<0) xlab<-'Cal. BC'
  	   }}
  	   title(xlab=xlab)
  }
  lines(x,U,lwd=1,col=col)
  lines(x,L,lwd=1,col=col)
  if(fill==TRUE) polygon(c(x, rev(x)), c(U, rev(L)), col = col, border = NA)
}

# Find the median of a set of SPD bootstraps (i.e. 'population proxy')

median.SPDboot<-function(SPDb) {
  y<-apply(SPDb,1,FUN='median')
  x<-as.numeric(rownames(SPDb))
  return(list(x=x,y=y)) }

SPDsigniftest<-function(SPDb,MOD) {
  e<-summary(SPDsignif(MOD / MOD[,ncol(MOD):1]))
  o<-summary(SPDsignif(MOD / SPDb))
  expected<-sum(e[3:4])
  observed<-sum(o[3:4])
  totalexp<-(e[1]/round(e[2]))*ncol(MOD)
  totalobs<-(e[1]/round(e[2]))*ncol(SPDb)
  ctab<-as.table(rbind(c(expected, totalexp-expected),c(observed, totalobs-observed)))
  dimnames(ctab)<-list(SPD=c('Model','Data'),c('Different','NotDiff'))
  return(chisq.test(ctab))
}

SPDsignif<-function(SPDb,probs=c(.05,.9), threshold=1) {
  SPDb<-summary(SPDb, probs=probs)
  sig_high<-SPDb[,1]>threshold
  sig_low<-SPDb[,2]<threshold
  out<-data.frame('year'=as.numeric(rownames(SPDb)), 'sig_high'=sig_high, 'sig_low'=sig_low)
  class(out)<-'SPDb_sig'
  return(out)
}

polygon <- function(x, ...) UseMethod("polygon")
polygon.default <- graphics::polygon

polygon.SPDb_sig<-function(S, add=TRUE, colhigh=rgb(1,0,0,.4), collow=rgb(0,0,1,.4)) {
  ylims<-par("usr")[3:4]
  xmax<-par("usr")[2]
  xlims_high<-c() ; if(S$sig_high[1]) xlims_high<-S$year[1]
  xlims_low<-c() ; if(S$sig_low[1]) xlims_low<-S$year[1]
  for(N in 2:length(S$year)) { 
    if( S$sig_high[N] & !S$sig_high[N-1] ) xlims_high<-c(xlims_high, S$year[N])
    if( S$sig_low[N] & !S$sig_low[N-1] ) xlims_low<-c(xlims_low, S$year[N])
    if( !S$sig_high[N] & S$sig_high[N-1] ) xlims_high<-c(xlims_high, S$year[N])
    if( !S$sig_low[N] & S$sig_low[N-1] ) xlims_low<-c(xlims_low, S$year[N])
  }
  if(length(xlims_high) %% 2 == 1) xlims_high<-c(xlims_high, xmax)
  if(length(xlims_low) %% 2 == 1) xlims_low<-c(xlims_low, xmax)
  polygon(
    x=rep(xlims_high, each=2), 
    y= rep(c(ylims,rev(ylims)),length(xlims_high)/2), 
    col=colhigh, border=NA)
  polygon(
    x=rep(xlims_low, each=2), 
    y= rep(c(ylims,rev(ylims)),length(xlims_low)/2), 
    col=collow, border=NA)
}

summary.SPDb_sig<-function(S) {
  out<-c(NA,NA,NA,NA)
  names(out)<-c('span','resolution','Cases sig. high', 'Cases sig. low')
  out[1]<-diff(range(S$year))
  out[2]<-mean(diff(S$year))
  out[3]<-sum(S$sig_high)
  out[4]<-sum(S$sig_low)
  return(out)
}


# Function to compute the ratio of two SPDboot classes

`/.SPDboot`<-function(A,B){
  # first establish what the output will look like
  # i.e. minimum of the number of input bootstaps, max of years
  Nboot<-min(dim(A)[2],dim(B)[2])
  yA<-as.numeric(rownames(A)); yB<-as.numeric(rownames(B))
  minyear<-min(c(yA,yB))
  maxyear<-max(c(yA,yB))
  #to do: introduce way of implementing other resolutions
  res<-10
  output_years<-seq(minyear,maxyear,res)
  #Make output
  out<-matrix(0,nrow=length(output_years),ncol=Nboot)
  rownames(out)<-output_years
  for(N in 1:Nboot) {
    bootA<-approx(x=yA, y=A[,N], xout=output_years)
    bootB<-approx(x=yB, y=B[,N], xout=output_years)
    out[,N]<-bootA$y/bootB$y
  }
  # Remove empty years and return
  out[is.na(out)]<-0
  empty<-rowSums(out)
  out<-out[empty>0,]
  class(out)<-'SPDboot'
  return(out)
}

 
# Function for converting the median of a 'SPDboot' SPDboot-like object back into rowcalsum 'c14sum' type
as.c14sum<-function(SPDb,FUN='median',...) {
	M<-t(apply(SPDb,1,FUN=FUN,na.rm=TRUE,...))
	out<-list(
		dates=NULL,
		stds=NULL,
		curves=NULL,
		yrs=as.numeric(rownames(SPDb)),		
		sum=M
	)
	class(out)<-'c14sum'
    return(out)
}

# User-friendly function to plot both sum and KDE a la CBR's 2017 paper
# Example use:
#  Step (1) Highlight and copy list of dates in spreadsheet
#  Step (2) rowcaldensum() 
# So.... simple to use, but will take some time for big lists with thousands of dates because all the calculations are performed each time the graphs are plotted

rowcaldensum<-function(DATA=CLIP(),N=100,perm_runs=1, perm_size=1, bw=30,col.fill=rgb(0,0,0.8,0.5),col.lines=rgb(0,0,0.8,1),col.den=rgb(0.8,0,0,0.4),...) {
	den<-MCdensity(DATA,N=N,perm_runs=perm_runs, perm_size=perm_size, bw=bw)
	Height<-max(den[,-1])*nrow(DATA)
	plot(den,col=col.den,scalefactor=nrow(DATA),ylim=c(0,Height*2),yaxt='n',...)
	tm<-pretty(c(0,Height)); tm<-tm[-length(tm)]
	axis(2,at=tm,label=tm)
        axis(4,at=Height+tm,label=tm)
	mat<-time_matrix(cbind(NA,DATA))
	psum(mat,add=TRUE,col.fill=col.fill,col.lines=col.lines,offset=Height)
}


# Generic radiocarbon animated map

# DATA must contain 5 columns 1: ID 2: BP 3: SD 4: X 5: Y
# Alternatively, you can specify via stmap= an object created by time_matrix, but the first 5 columns must be as above
# threshold refers to a minium p-value for plotting date
# Example:
# ire<-readShapeSpatial("../GISMisc/Britain/Ireland_OSGB.shp")
# data<-BIRE[BIRE$Where=='Ireland' & BIRE$mcode=='Sc' & BIRE$ccode=='I',c(1,2,3,6,7)]
# make_animated_map(data,ire,out='Cereals_from_Industry.pdf')

make_animated_map<-function(DATA,coast,stmap=NA,out='output.pdf',pt.cex=1,pt.col=rgb(1,0,0,0.4),threshold=0.002,dropempty=FALSE,BP=FALSE,width=6,height=6,...) {
  pdf(out,width=width,height=height)
  print('Calibrating dates...')
  if(is.na(stmap)) stmap<-time_matrix(DATA)
  print('Drawing map')
  pb <- txtProgressBar(min=6,max=ncol(stmap),initial=6)
  for(COL in 6:ncol(stmap)){
    if(dropempty==FALSE | sum(stmap[,COL]>0)) {
       plot(coast,...)
       year<-as.numeric(colnames(stmap)[COL])
       year_str<-paste(-year,"cal. BC")
       if(year > 0 ) year_str<-paste(year,"cal. AD")
       if(BP==TRUE) year_str<-paste(1950-year,"cal. BP")
       submap<-na.omit(stmap[stmap[,COL]>threshold,c(4:5,COL)])
       if(nrow(submap)>0) { for (n in 1:nrow(submap)) {
            # now plot the point, size scaled from the prob
            psize<-pt.cex*25*submap[n,3]^0.5     # scaling factor
            if (psize > 3.5) psize<-3.5   # sets a maximum size 
            points(submap[n,1],submap[n,2],pch=15,col=pt.col,cex=psize)
    }}}
    title(main=year_str)
    setTxtProgressBar(pb,COL)
	}
	close(pb)
	dev.off()
}


DMlm<-function(DM, from=0, to=0, add=FALSE, ...) {
	# Select relevant part of KDE model 
	if(from==0 & to==0) {
		loc<-locator()
		from<-loc$x[1]; to<-loc$x[2] }
	fromi<-which(abs(DM[,1]-from)==min(abs(DM[,1]-from)))
	toi<-which(abs(DM[,1]-to)==min(abs(DM[,1]-to)))
	SUB<-DM[fromi:toi,]
	
	# Compute linear regression for each KDE simulation
	
	rate<-c(); Pr<-c(); r.squared<-c()
	for(N in 2:ncol(SUB)){
		mod<-lm(SUB[,N] ~ SUB[,1])
		slm<-summary(mod)	
		rate[N]<-slm$coefficients[2,1]
		 Pr[N]<-slm$coefficients[2,4]
		r.squared[N]<-slm$r.squared
		if(add) abline(mod,...)
		}
	return(data.frame(rate=rate, Pr=Pr, r.squared=r.squared))
}

DMgrowth<-function(DM, from=0, to=0, add=FALSE, ...) {
	# Select relevant part of KDE model 
	if(from==0 & to==0) {
		loc<-locator()
		from<-loc$x[1]; to<-loc$x[2] }
	fromi<-which(abs(DM[,1]-from)==min(abs(DM[,1]-from)))
	toi<-which(abs(DM[,1]-to)==min(abs(DM[,1]-to)))
	SUB<-DM[fromi:toi,]
	
	# Compute linear regression for each KDE simulation
	# And calculate an annual growth rate for each
	growth<-c(); Pr<-c(); r.squared<-c(); int<-c()
	for(N in 2:ncol(SUB)){
		mod<-lm(SUB[,N] ~ SUB[,1])
		slm<-summary(mod)	
		 k<-slm$coefficients[2,1]
   		 C<-slm$coefficients[1,1]
    	 y2<-to*k+C
    	 y1<-from*k+C
    	if(y2<y1) growth[N]<-(-y1/y2)/(to-from) else growth[N]<-(y2/y1)/(to-from)
    	#growth[N]<-(y2-y1)/(y1*(to-from))
		#growth[N]<-(((y2-y1)/y1)*100)/(to-from)
		r.squared[N]<-slm$r.squared
		Pr[N]<-slm$coefficients[2,4] 
		if(add) abline(mod,...)
		}
	out<-data.frame(from=from,to=to,growth=growth, Pr=Pr, r.squared=r.squared)
	out<-out[-1,]
	class(out)='DMfit'; return(out)
}

summary.DMfit<-function(DMf) {
   stopifnot(inherits(DMf, "DMfit"))
   from<-min(DMf$from)
   to<-max(DMf$to)
   g<-sort(DMf$growth)
   # Exclude the 30 strangest regressions
   l<-length(g); g<-g[15:l-15]
   cat("\n", 
        sprintf("\tSpan: %s to %s",round(from,0),round(to,0)),"\n",
        sprintf("\tGrowth rate (R-squared): %s±%s (%s)",round(median(g)*100,2),round(mad(g)*100,2),round(mean(DMf$r.squared),2)),"\n",
        sprintf("\tPr|t| : %s±%s",signif(mean(DMf$Pr),2),signif(sd(DMf$Pr),2)),"\n")
}

totalden<-function(DM, from, to, density=FALSE) {
   tot<-sum(DM[,-1])/(ncol(DM)-1)
   sub<-DM[which(DM[,1]>=from & DM[,1]<to),-1]
   r<-(sum(sub)/ncol(sub))/tot
   if(density) r<-r/(to-from)	
   return(r)
}

# Function to build a timeseries model by drawing on the current plot
       
drawmodel<-function(n=100,add=TRUE,...) {
	model_raw<-locator(n)
	mod<-approx(x=model_raw$x, y=model_raw$y, xout=round(min(model_raw$x)):round(max(model_raw$x)), rule=2)
	if (add) lines(mod,...)
	mod
}

# Population simulation

pop_sim<-function(rate=1,startpop=NA,from=NA,to=1700){
	if(is.na(startpop) | is.na(from)) {
		l<-locator()
		from=l$x[1]
		startpop=l$y[1] }
	pop<-c()
	pop[1]<-startpop; Y=0
	for(Y in 2:(to-from+1)) pop[Y]<-pop[Y-1]*(rate/100)+pop[Y-1]
	return(data.frame(x=c(from:to),y=pop))
}


# Function that adds the calibration curve across the current plot

lines.calcurve<-function(calcurve=intcal, BP=FALSE,...) {
	xr<-par('usr')[1:2]; yr<-par('usr')[3:4]
	if(BP==FALSE) mod=1950 else mod=0
	calx<-calcurve[calcurve[,1]<(mod-xr[1]) & calcurve[,1]>mod-xr[2],]
	calx$y<-(calx[,2]-min(calx[,2]))/diff(range(calx[,2]))*yr[2]
	calx$err<-calx[,3]/diff(range(calx[,2]))*yr[2]
	lines(mod-calx[,1],calx$y+calx$err,...)
	lines(mod-calx[,1],calx$y-calx$err,...)
	calx$err<-calx[,3]/diff(range(calx[,2]))*yr[2]
}


# Function to transform an object of class 'MCd' into geometric growth rate (in % per year)
#     via numeric differentiation

ggr<-function(MCD, bw=100) {
   out<-MCD[-nrow(MCD),]
   yrs<-MCD[,1]
   for(N in 2:ncol(MCD)) out[,N]<- ( (diff(MCD[,N]) / MCD[-nrow(MCD),N]) / diff(MCD[,1])) *100
   class(out)<-'diffMCd'
   return(out)
}

# Method for plotting same

plot.diffMCd<-function(MCD, style=2, ..., ylab='Growth rate / %', ylim=c(-2,2))  plot.MCd(MCD, style=style, ylab=ylab, ylim=ylim, ...)


# Function to summerise a MCd or diffMCd object
 
summary.MCd<-function(MCD,probs=c(.05,.9)) {
  out<-MCD[,1:(length(probs)+1)]
  qu<-t(apply(MCD[,2:ncol(MCD)],1,'quantile',probs=probs,na.rm=TRUE))
  out[,2:(length(probs)+1)]<-qu
  return(out)
}	

# function to offset the timebase of a MCd or diffMCd object
timeshift<-function(MCD, shift=0) {
	original_class<-class(MCD)
	out<-MCD
	out[,1]<-out[,1]+shift
	class(out)<-original_class
	return(out)
}

# function to find all the years of significant +ve or -ve growth in a differentiated Monte Carlo density model (diffMCd)

ggrsignif<-function(MCD,probs=c(.05,.9), threshold=0) {
  MCD<-summary.MCd(MCD, probs=probs)
  sig_high<-MCD[,2]>threshold
  sig_low<-MCD[,3]<threshold
  out<-data.frame('year'=MCD[,1], 'sig_high'=sig_high, 'sig_low'=sig_low)
  class(out)<-'SPDb_sig'
  return(out)
}

# return list of weighted means of samples, given three column input (BP, SD, Calcurve)
findmixmedian<-function(L,...) {
   O<-c()
   for(N in 1:nrow(L)) eval(parse(text=paste('O[N]<-rowcalmedian(L[N,1],L[N,2],calcurve=',L[N,3],')',sep='')))
   return(O)
}

# Function to find all the bootstrapped years of stationary growth
stationary<-function(MCD, from=NULL, to=NULL, tolerance=0.005) {
   out<-c()
   if(!is.null(from)) MCD<-MCD[MCD[,1]>from,]
   if(!is.null(to)) MCD<-MCD[MCD[,1]<to,]
   yrs<-MCD[-nrow(MCD),1]
   for(N in 2:ncol(MCD)) {
       bootstrap <- ( (diff(MCD[,N]) / MCD[-nrow(MCD),N]) / diff(MCD[,1])) *100
       out<-c(out,yrs[bootstrap>-tolerance & bootstrap<tolerance])
   }
   return(out)
}


# Function to bootstrap a histogram distribution, much like MCdensity but in the 
#      format of a historgram where an error bar representes calibration uncertainty
# Bin width is specified by the parameter 'bw'

# Plot a histogram (if plot=TRUE) or outputs the matrix used to calculate it (if plot=FALSE)

MChist<-function(datelist, bw=100, buffer=300, Nboot=100, default_curve='intcal', plot=TRUE, xlab='Cal. BC/AD', col.error='black',lwd.error=1,lty.error=1, ... ) {
	
 # find range of input 

  if (ncol(datelist)<2) stop("MChist needs at two or three columns in the input (BP, error, [calcurve])")
  if (ncol(datelist)==2) datelist[,3]<-default_curve 
  lower_date<-min(rowcalmode(datelist[(datelist[,1]==max(datelist[,1])),1],40))-buffer
  upper_date<-max(rowcalmode(datelist[(datelist[,1]==min(datelist[,1])),1],40))+buffer
  if (upper_date==-Inf) upper_date<-2000

  # build a range of round numbers based on the required resolution 'res' 
  bins<-seq(round(lower_date,-round(log10(bw))),round(upper_date,-round(log10(bw))),by=bw)
  
  # make temporary histogram structure and matrix for the bootstraps
  A<-hist(findmixmedian(datelist), breaks=bins, plot=FALSE)
  midpoints<-A$mids
  out<-matrix(NA, nrow=length(midpoints), ncol=Nboot)
  rownames(out)<-midpoints  

  # Do the bootstrap resampling via MCmix
  # bootstrap wrapped in try because sometimes an outling value is sampled, which breaks hist()
  # if a lot of these error messages are recieved, try the function again with a bigger buffer
  # e.g, MChist(mylist, buffer=500)
  pb <- txtProgressBar(min=1,max=Nboot,initial=1)
  for(N in 1:Nboot) {
       try(out[,N]<-hist(MCmix(datelist),breaks=bins,plot=FALSE)$counts)
       setTxtProgressBar(pb,N)
  }

  # Compute summary stats and store this in histogram structure
  M<-rowMeans(out,na.rm=TRUE)
  Sd<-apply(out,1,sd,na.rm=TRUE)
  A$counts<-M  
  A$bootstrapped_counts<-out
  class(A)<-'MChistogram'

  #plot the output
  if(plot) plot(A, xlab=xlab, ...)      
  else return(A)
}

# Method for plotting the output of MChist

plot.MChistogram<-function(A, xlab='Cal. BC/AD', col.error='black',lwd.error=1,lty.error=1, ... ) {
      class(A)<-'histogram'
      M<-rowMeans(A$bootstrapped_counts,na.rm=TRUE)
      Sd<-apply(A$bootstrapped_counts,1,sd,na.rm=TRUE)
      plot(A, main='', xlab=xlab, ...)      
      for(N in 1:length(M)) {
          x=rep(A$mids[N],2)
          y=c(M[N]-Sd[N],M[N]+Sd[N])
          if(y[1]<0) y[1]<-0
          lines(x,y,col=col.error,lwd=lwd.error,lty=lty.error)
      }
}

# Method for plotting density models akin to the 'battleship' seriation graphs of yore
# see 'battleships() below'

battleship.MCd<-function(x, os=0, add=FALSE,col=rgb(0,0,0.8,0.5),lcol=rgb(0,0,0.8,0.8),lwd=1, lty=1, scalefactor=1,ylab='Density',xlab='Automatic',grid=TRUE,...) { 
  M<-rowMeans(x[,2:ncol(x)],na.rm=TRUE)*scalefactor
  Sd<-apply(x[,2:ncol(x)],1,sd,na.rm=TRUE)*scalefactor
  if (!add) {
    ylm<-max(pretty(rowMeans(x[,2:ncol(x)])))
  	plot(x[,1],M,col=NA,xlab=NA,ylab=ylab,ylim=c(-ylm,ylm), ...)
  	if(grid) grid(lwd=0.5)
    if(!is.na(xlab)) { if(xlab=='Automatic') { 
      xlab<-'Cal. BC/AD'
      pw<-par('xaxp')
      if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
      if(pw[2]<0) xlab<-'Cal. BC'
  	}}
  	title(xlab=xlab)
  }
  mn<-(M-Sd)
  mx<-(M+Sd)
  polygon(c(x[,1], rev(x[,1])), c(mn+os, rev(-mn)+os), col = col, border = NA)
  lines(x[,1],mx+os, lwd=lwd, lty=lty, col=lcol)
  lines(x[,1],-mx+os,lwd=lwd, lty=lty, col=lcol)
}  

# function to plot multiple 'battleship' curves for one or a list of MCd objects
# e.g.   >battleships(mydenmod1) 
# or     >battleships(list('North'=mydenmod1, 'South & East'=mydenmod2, 'West'=mydenmod3))

battleships<-function(L, cols=rep(c(1:9),10), space=1.5, ylas=1, xlab='Automatic', ylab='', ...) {
   if(class(L)=='MCd') L<-list(x=L)
   if(! class(L) %in% c('MCd', 'list')) stop('input should be a MCd object or a list of MCd objects')
   Nships<-length(L)

   # find age range
   x<-c()
   for(N in 1:Nships) x<-c(x, min(L[[N]][,1]), max(L[[N]][,1]))
   
   # find probability range
   ylm<-c()
   for(N in 1:Nships) ylm<-c(ylm,max(pretty(rowMeans(L[[N]][,2:ncol(L[[N]])]))))

   # Set up plot
   plot(c(min(x),max(x)), c(-ylm[1],max(ylm)*Nships*space), col=NA, xlab=NA,ylab=ylab, yaxt='n', ...)
   if(!is.na(xlab)) { if(xlab=='Automatic') { 
     xlab<-'Cal. BC/AD'
     pw<-par('xaxp')
     if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
     if(pw[2]<0) xlab<-'Cal. BC'
  	}}
  	title(xlab=xlab)
   
   # Plot the ships  
   for(N in 1:Nships) battleship.MCd(L[[N]], col=cols[N], lcol=cols[N], os=(N-1)*space*max(ylm), add=TRUE)
   
   # Label y-axis
   axis(2,at=(c(1:Nships)-1)*space*max(ylm), lab=names(L), las=ylas)
}
	    
