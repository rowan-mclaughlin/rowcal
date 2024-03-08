# Calibrate radiocarbon date with output in rows of two columns
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

# Function for mixing curves, or applying local delta-R corrections, or both

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


#function to find modal calibrated date for a C14 date, based on rowcal above
rowcalmode<-function(date,sigma,...) {
 cd<-rowcal(date,sigma,...)
 cd[,1][cd[,2]==max(cd[,2])]
}


#function to find weighted mean date for a C14 date, based on rowcal above
rowcalwm<-function(date,sigma,...) {
   g<-rowcal(date,sigma,...)
   sum(g[,1]*g[,2]/max(g[,2]))/sum(g[,2]/max(g[,2]))
}


#function to find (even just one) age sample for C14 date, based on rowcal above
rowcalsam<-function(date, sigma, calcurve=intcal, N=1, ...) {
  g<-rowcal(date, sigma, calcurve, ...)
  random.points <- approx(cumsum(g[,2])/sum(g[,2]),g[,1],runif(N))$y
  random.points
}


#function to find median of a C14 date
rowcalmedian<-function(date, sigma, ...) {
   d<-rowcal(date,sigma,...)
   return(d[which.min(abs(cumsum(d[,2])-0.1)),1])
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
#third column must be the name of the calibration curve, or 'cal' meaning calendar
MCmix<-function(L) {
	colnames(L)<-c('date','error','cc')
	cal<-L[L$cc=='cal',]
	cal_points<-c()
	if(nrow(cal)>0) cal_points<-apply(cal[,c('date','error')],1,function(rdate,...) calsam(rdate['date'],rdate['error'],...) )
    C14_points<-c()
    C14<-L[!(L$cc %in% c('cal')),]
    if(nrow(C14)>0) {
    # I know, this is as ugly as sin
	for (N in 1:length(C14[,1])) 
		C14_points[N]<-eval(parse(text=paste("rowcalsam(",C14[N,1],",",C14[N,2],",calcurve=",as.character(C14[N,3]),")",sep='')))
	}
	return(as.vector(c(cal_points,C14_points)))
}	


#function to turn a two-column list of uncalibrated dates into list of one sample per date
MCsam<-function(L,...) {
   colnames(L)<-c('BP','SD')
   O<-apply(L[,c('BP','SD')],1,function(rdate,...) rowcalsam(rdate['BP'],rdate['SD'],...) )
   O
}

#version for already-calibrated lists
MCsam.list<-function(L){
  out<-c()
  for(N in 1:length(L)) {
    g<-L[[N]]
    out[N]<-approx(cumsum(g[,2])/sum(g[,2]),g[,1],runif(1))$y
  }
  out
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

# return list of weighted means of samples, given three column input (BP, SD, Calcurve)
findmixmedian<-function(L,...) {
   	colnames(L)<-c('date','error','cc')
	spp<-L[L$cc=='sap',]; cal<-L[L$cc=='cal',]
	cal_points<-c(); sap_points<-c()
	if(nrow(spp)>0) sap_points<-apply(spp[,c('date','error')],1,row,...) 
	if(nrow(cal)>0) cal_points<-apply(cal[,c('date','error')],1,function(rdate,...) calsam(rdate['date'],rdate['error'],...) )
    C14_points<-c()
    C14<-L[!(L$cc %in% c('sap','cal')),]
    if(nrow(C14)>0) {
    # I know this is as ugly as sin, but it works
	for (N in 1:length(C14[,1])) 
		C14_points[N]<-eval(parse(text=paste("rowcalmedian(",C14[N,1],",",C14[N,2],",calcurve=",as.character(C14[N,3]),")",sep='')))
	}
	return(as.vector(c(cal_points,sap_points,C14_points)))
}	
 

# function for inputing data from clipboard independent of platform
CLIP<-function() {
	pasted<-character()
	if (Sys.info()[1]=='Darwin') pasted<-read.table(pipe('pbpaste'),head=FALSE) else pasted<-read.delim('clipboard',head=FALSE)
	pasted
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
# alternatively, curve 'cal' for a mean of calendar years which will be sampled from the gaussian
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
        #if(count==10) stop('cannot calibrate dates')
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

# Compute density estimate like MCdensity via the above phase sampler
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
  if(mix)  d<-density(MCmix(phasesam(DATA, h=h)[,2:4]),n=512, na.rm=TRUE)
  if(!mix) d<-density(MCsam(phasesam(DATA, h=h)[,2:3]),n=512, na.rm=TRUE)
  if (plot.new) {add<-TRUE; plot(d,col=col,...)}
  x1<-min(round(d$x))-100; x2<-max(round(d$x))+100; n=d$n
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


# Method for plotting all of the above density models
plot.MCd<-function(x, add=FALSE,col=1,fill='#00000022',scalefactor=1,smax=FALSE,ylab='Density',xlab='Automatic',Toffset=0,grid=FALSE,lwd=1, lty=1, ...) { 
  x[,1]<-x[,1]+Toffset
  if(smax) {
    M<-rowMeans(x[,-1],na.rm=TRUE)
    scalefactor=1/max(M)
  }
  M<-rowMeans(x[,-1],na.rm=TRUE)*scalefactor
  Sd<-apply(x[,-1],1,sd,na.rm=TRUE)*scalefactor
  if(smax) scalefactor=1/max(M)
  if (!add) {
  	plot(x[,1],M+Sd,col=NA,xlab=NA,ylab=ylab,...)
  	if(grid) grid(lwd=0.5)
    if(!is.na(xlab)) { if(xlab=='Automatic') { 
      xlab<-'Cal. BC/AD'
      pw<-par('xaxp')
      if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
      if(pw[2]<0) xlab<-'Cal. BC'
      if(abs(mean(x[,1]))<100) xlab<-'kyr BP'
  	}}
  	title(xlab=xlab)
  }
  polygon(c(x[,1], rev(x[,1])), c(M+Sd, rev(M-Sd)), col = fill, border = NA)
  lines(x[,1],M,lwd=lwd,lty=lty,col=col)
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

plot_MCd_offsets<-function(x, norm=TRUE, add=FALSE,col=1,fill='#e5e5e5e5',scalefactor=1,
                           ylab='Density',xlab='Automatic',
                           yaxt='Automatic',yaxt_side=2, Toffset=0,yoffset=0,yaxes=TRUE,
                           grid=FALSE,lwd=1, lty=1, ...) { 
  x[,1]<-x[,1]+Toffset
  if(norm) x[,2:ncol(x)]<-x[,2:ncol(x)]/max(x[,2:ncol(x)])
  x[,2:ncol(x)]<-x[,2:ncol(x)]+yoffset
  M<-rowMeans(x[,2:ncol(x)],na.rm=TRUE)*scalefactor
  Sd<-apply(x[,2:ncol(x)],1,sd,na.rm=TRUE)*scalefactor
  if (!add) {
    plot(x[,1],M,col=NA,xlab=NA,ylab=ylab,yaxt='n', ...)
    if(grid) grid(lwd=0.5)
    if(!is.na(xlab)) { if(xlab=='Automatic') { 
      xlab<-'Cal. BC/AD'
      pw<-par('xaxp')
      if(pw[1]>0 & pw[2]>0) xlab<-'Cal. AD'
      if(pw[2]<0) xlab<-'Cal. BC'
    }}
    title(xlab=xlab)
  }
  polygon(c(x[,1], rev(x[,1])), c(M+Sd, rev(M-Sd)), col = fill, border = NA)
  lines(x[,1],M,lwd=lwd,lty=lty,col=col)
  signifs<-round(abs(log10(max(x[,-1]))))+1
  locs<-round(pretty(rowMeans(x[,-1])-yoffset),signifs)
  if(yaxes) axis(side=yaxt_side, at=locs+yoffset, lab=locs)
}


# Standard deviation method for density models
# 
# 'sd' is not a generic function so need to hijack its default
sd <- function(x, ...) UseMethod("sd")
sd.default <- stats::sd

sd.MCd<-function(A, from=NULL, to=NULL) { 
    if(is.null(from)) from<-min(A[,1])
    if(is.null(to)) to<-max(A[,1])
    return(sd(A[A[,1]>=from & A[,1]<=to, 2:ncol(A)], na.rm=TRUE))
}

# Log transform a density model
log.MCd<-function(x){
  out<-x
  out[,2:ncol(out)]<-log(out[,2:ncol(out)])
  class(out)='MCd'
  return(out)
}


#Wee functions for polygon plots
fill<-function(j,col="grey",border=NA,...) polygon(c(j[1,1],j[,1],max(j[,1])), c(j[1,2],j[,2],j[1,2]),col=col,border=border,...)

# 'Time matrix' function for summing etc

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



#Calculate SPD for radiocarbon dates in a flexible format
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

# Function to transform an object of class 'MCd' into geometric growth rate (in % per year)
#     via numeric differentiation
ggr<-function(MCD, threshold=6e-6) {
   out<-MCD[-nrow(MCD),]
   yrs<-MCD[,1]
   M<-rowMeans(MCD[,-1],na.rm = TRUE)
   # filter: change tiny density estimates to zero
   MCD[which(M<threshold),2:ncol(MCD)] <- 0
   for(N in 2:ncol(MCD)) out[,N]<- ( (diff(MCD[,N]) / MCD[-nrow(MCD),N]) / diff(MCD[,1])) *100
   class(out)<-'diffMCd'
   return(out)
}

# Method for plotting same
plot.diffMCd<-function(MCD, ..., endzero=FALSE, ylab='Growth rate / %', ylim=c(-2,2)) {
  #add 0 at beginning and end of density model (maybe incoporate this to ggr() instead?)
  if(endzero) MCD<-rbind(c(MCD[1,1]-diff(MCD[1:2,1]),rep(0,ncol(MCD)-1)),MCD,c(MCD[nrow(MCD),1]+diff(MCD[1:2,1]),rep(0,ncol(MCD)-1)))
  M<-rowMeans(MCD)
  zerorows<-which(M<1e-7)
  #send only finite values to plot.MCd
  plot.MCd(MCD[which(is.finite(M)),], ylab=ylab, ylim=ylim, ...)
}


# Function to summerise a MCd or diffMCd object
summary.MCd<-function(MCD,probs=c(.05,.9)) {
  out<-MCD[,1:(length(probs)+1)]
  qu<-t(apply(MCD[,2:ncol(MCD)],1,'quantile',probs=probs,na.rm=TRUE))
  out[,2:(length(probs)+1)]<-qu
  return(out)
}	


# function to find all the years of significant +ve or -ve growth in a differentiated Monte Carlo density model (diffMCd)
ggrsignif<-function(MCD,probs=c(.05,.9), threshold=0) {
  MCD<-summary.MCd(MCD, probs=probs)
  sig_high<-MCD[,2]>threshold
  sig_low<-MCD[,3]<threshold
  sig_high[which(is.na(sig_high))]<-FALSE
  sig_low[which(is.na(sig_low))]<-FALSE
  out<-data.frame('year'=MCD[,1], 'sig_high'=sig_high, 'sig_low'=sig_low)
  class(out)<-'SPDb_sig'
  return(out)
}

# Method for plotting same. First hijack default method
polygon <- function(x, ...) UseMethod("polygon")
polygon.default <- graphics::polygon

polygon.SPDb_sig<-function(S, add=TRUE, colhigh='#FF000020', collow='#0000FF20') {
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

# Function for drawing pretty axes with tickmarks
ax<-function(side=1, tick=100, ticksize=-0.01, labs=NULL) {
  ats<-pretty(par('usr')[1:2])
  if(is.null(labs)) labs<-c(abs(ats[which(ats< -1)]),ats[which(ats>-1)]+1)
  axis(side, at=ats,lab=labs)
  rug(seq(ats[1]-500,ats[length(ats)]+500,tick),ticksize = ticksize, side=side,quiet = TRUE)
}
