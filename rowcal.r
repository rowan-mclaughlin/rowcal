#   r o w C a l   

#   Utilities for radiocarbon age calibration in R
#   by T. Rowan McLaughlin email : r.mclaughlin@qub.ac.uk

#   Please cite McLaughlin, T. R. Submitted. On applications of space-time modelling with open-source 14C age calibration. Journal of Archaeological Method and Theory

intcal<-read.csv("rowcal/intcal13.14c",header=FALSE,skip=11)
marine<-read.csv("rowcal/marine13.14c",header=FALSE,skip=11)
colnames(intcal)<-c("calBP","14Cage","Error","Delta14C","Sigma")
colnames(marine)<-c("calBP","14Cage","Error","Delta14C","Sigma")

rowcal<-function(date, sigma, calcurve=intcal, BC=TRUE) {
    curve<-calcurve[(calcurve[,2]<= (date+(sigma*4)) & calcurve[,2]>= (date-(sigma*4)) ),]
	
	#convert from cal. BP to cal. BC if required
    if (BC==TRUE) curve$calBP<-1950-curve$calBP
    
    #build an output table
    calibrated<-as.matrix(curve[,1:2]) 
    if (BC==TRUE) colnames(calibrated)<-c("Cal. BC/AD","Relative p") else colnames(calibrated)<-c("Cal. BP","Relative p")
    
    #apply the normal distribution function to each margin of the calibration curve and average
    H<-dnorm(curve[,2]+curve[,3],mean=date,sd=sigma)
    L<-dnorm(curve[,2]-curve[,3],mean=date,sd=sigma)
    ol<-length(curve[,1])
    Res<-((curve[,1]-curve[1,1])/ol-1)[ol]
    calibrated[,2]<-(H+L)/2/Res/sum((H+L)/2)
    calibrated
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
	pb <- txtProgressBar(min=2,max=N-1,initial=2)
	for (run in 2:N-1) {
		P<-P+1
		if (perm_runs>1 & perm_size<1) DATA<-oDATA[sample(1:nrow(oDATA),round(perm_size*nrow(oDATA))),]	
	    d<-density(MCsam(DATA),bw,na.rm=TRUE,from=x1,to=x2,n=512)
		out[,P+2]<-d$y 
		if (add & P>2) lines(out[,1],rowMeans(out[,2:P]),col=col)
		setTxtProgressBar(pb,run)
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
    sub1<-D1[which(round(D1[,1],-1) %in% round(D2[,1],-1)),]
    sub2<-D2[which(round(D2[,1],-1) %in% round(D1[,1],-1)),]
    if(nrow(sub1)!=nrow(sub2)) stop("Sorry can't align inputs")
    plot(sub1[,1],sub1[,2]/sub2[,2],col=NA,...)
    for (N in 2:ncol(sub1)-1) lines(sub1[,1],sub1[,N]/sub2[,N],col=rgb(0,0,0,0.05)) 
    abline(h=1,col='red')
}

#Extract the various probablity densities for a given year from an object of class MCd
year_densities<-function(Y,MCd) MCd[which.min(abs(txd[,1]-Y)),][-1]
	
#Wee functions for polygon plots, plot.C14 is the method of plot()-ing dates
fill<-function(j,col="grey",border=NA) polygon(c(j[1,1],j[,1],max(j[,1])), c(j[1,2],j[,2],j[1,2]),col=col,border=border)
plot.C14<-function(x,col="grey",...) {plot(x,col=NA,...) ; fill(x,col=col)}

# Obtain confidence intervals using Highest Density Regions (requires hrdcde package)
try(library(hdrcde))

spansc14<-function(rdate) {
	d<-list(x=rdate[,1],y=rdate[,2])
	spans<-hdr(den=d,prob=c(68,95))$hdr
	spans
}

summary.C14<-function(rdate) {
	d<-list(x=rdate[,1],y=rdate[,2])
	spans<-hdr(den=d,prob=c(68.27,95.45))
	spans
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
  
  #perform the calculations, by default if the fourth column contains M for marine samples the marine curve will be used
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

#Optimised time_matrix function
time_matrix2<-function(datelist,buffer=100,BC=TRUE) {
  
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
  
  #perform the calculations, by default if the fourth column contains M for marine samples the marine curve will be used
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
   if (lines==FALSE && plot==FALSE && fill==FALSE) return(s+offset)
}  

# User-friendly function to quickly sum any list of dates
rowcalsum<-function(DATA=CLIP(),...) {mat<-time_matrix(cbind(NA,DATA)); psum(mat,...)}

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

make_animated_map<-function(DATA,coast,stmap=NA,out='output.pdf',pt.cex=1,pt.col=rgb(1,0,0,0.4),threshold=0.002,width=6,height=6,...) {
  pdf(out,width=width,height=height)
  print('Calibrating dates...')
  if(is.na(stmap)) stmap<-time_matrix(DATA)
  print('Drawing map')
  pb <- txtProgressBar(min=6,max=ncol(stmap),initial=6)
  for(COL in 6:ncol(stmap)){
    plot(coast,...)
    year<-as.numeric(colnames(stmap)[COL])
    year_str<-paste(-year,"cal. BC")
    if(year > 0 ) year_str<-paste(year,"cal. AD")
    submap<-na.omit(stmap[stmap[,COL]>threshold,c(4:5,COL)])
    if(nrow(submap)>0) { for (n in 1:nrow(submap)) {
         # now plot the point, size scaled from the prob
         psize<-pt.cex*25*submap[n,3]^0.5     # scaling factor
         if (psize > 3.5) psize<-3.5   # sets a maximum size 
         points(submap[n,1],submap[n,2],pch=15,col=pt.col,cex=psize)
    }}
    title(main=year_str)
    setTxtProgressBar(pb,COL)
	}
	close(pb)
	dev.off()
}
