#   r o w C a l   

#   Utilities for radiocarbon age calibration in R
#   by T. Rowan McLaughlin email : r.mclaughlin@qub.ac.uk

#   Please cite McLaughlin, T. R. 2018. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory

intcal<-read.csv("rowcal/intcal13.14c",header=FALSE,skip=11)
marine<-read.csv("rowcal/marine13.14c",header=FALSE,skip=11)
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

plot.MCd<-function(x,add=FALSE,col=rgb(0,0,0.8,0.5),scalefactor=1,xlab='Automatic',grid=TRUE,...) { 
  M<-rowMeans(x[,2:ncol(x)],na.rm=TRUE)*scalefactor
  Sd<-apply(x[,2:ncol(x)],1,sd,na.rm=TRUE)*scalefactor
  if (!add) {
  	plot(x[,1],M,col=NA,xlab=NA,ylab='Density',...)
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

#To test for convergenge:

#test<-MCdensity(SOME_DATA,N=100)
#conv<-c(); nc<-ncol(test)-1
#T<-rowMeans(test[,2:nc+1],na.rm=T)
#for (N in 4:nc-1) conv[N]<-sum(T,na.rm=T)-sum(rowMeans(test[,2:N],na.rm=T))
#plot(conv,pch=1,col=rgb(0.7,0.1,0,0.5),xlab="Run number", ylab="Convergence")
#abline(h=0,col='grey')

MCsimulate<-function(mod,Nd=1000,...) {
   S<-approx(cumsum(mod$y)/sum(mod$y),mod$x,runif(Nd))$y
   # by default 1000 random simulated dates drawn from 'mod'
   S<-S[!is.na(S)]
   datelist<-matrix(nrow = length(S),ncol=2)
   datelist[,1]<-sapply(S,'revcal')
   datelist[,2]<-sample(25:40,length(S),replace=TRUE)
   out<-MCdensity(datelist,...)
   class(out)='MCd'; return(out)
}

drawmodel<-function(n=100,add=TRUE,...) {
	model_raw<-locator(n)
	mod<-approx(x=model_raw$x, y=model_raw$y, xout=round(min(model_raw$x)):round(max(model_raw$x)), rule=2)
	if (add) lines(mod,...)
	mod
}

pop_sim<-function(rate=1,startpop=10000,from=200,to=1700){
	pop<-c()
	pop[1]<-startpop; Y=0
	for(Y in 2:(to-from+1)) pop[Y]<-pop[Y-1]*(rate/100)+pop[Y-1]
	return(data.frame(x=c(from:to),y=pop))
}

#e.g. try the following to add an exponential population model with annual increase of 0.5%
#the plot will start from a point clicked upon with the mouse
#l<-locator(); lines(pop_sim(rate=0.5,from=l$x[1],to=1700,start=l$y[1]))
	
#function to make a MCdensity-like object from a list of calender dates without 14C uncertainty
#using bootstrapping (sampling with replacement)
bootstrap_density<-function(DATA=CLIP(),N=100, col=rgb(0,0,0,0.01),bw=30,...) {
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

compare_densities<-function(D1,D2,...) {
    sub1<-D1[which(round(D2[,1],-1) %in% round(D1[,1],-1)),]
    sub2<-D2[which(round(D1[,1],-1) %in% round(D2[,1],-1)),]
    if(nrow(sub1)!=nrow(sub2)) stop("Sorry can't align inputs")
    plot(sub1[,1],sub1[,2]/sub2[,2],col=NA,...)
    for (N in 2:ncol(sub1)-1) lines(sub1[,1],sub1[,N]/sub2[,N],col=rgb(0,0,0,0.05)) 
    abline(h=1,col='red')
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

# Wee functions to print the upper and lower bounds of a calibrated age
## NB for some reason these sometimes fail when mapply-ing, presumably due to a bug in hdrcode (?)
## Best to apply these to a dataframe in a loop wrapped with try(), e.g.: for (n in 1:length(b[,1])) try(b$L68[n]<-L68(b[n,2],b[n,3]))

U95<-function(d,e) round(max(spansc14(rowcal(d,e))[1,], na.rm=TRUE),-1)
L95<-function(d,e) round(min(spansc14(rowcal(d,e))[1,], na.rm=TRUE),-1)
U68<-function(d,e) round(max(spansc14(rowcal(d,e))[2,], na.rm=TRUE),-1)
L68<-function(d,e) round(min(spansc14(rowcal(d,e))[2,], na.rm=TRUE),-1)

# 'Time matrix' function for summing / animated maps
# Returns a time matrix from input, which is a two-column dataframe
# The argument -- buffer -- describes how much time to wrap round the matrix 

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



# code to sum and plot sums of C14 dates
# startcol can be specified to point to the first column of the input matrix that contains the radiocarbond data,
# by default this is the first numeric column name

# 'offset' can be specified to 'float' the sum into the yaxis some distance

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
rowcalsum_old<-function(DATA=CLIP(),...) {mat<-time_matrix(cbind(NA,DATA)); psum(mat,...)}

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

