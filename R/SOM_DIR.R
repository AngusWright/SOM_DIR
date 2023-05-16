#!/usr/bin/Rscript
#
# Script to generate DIR weights from an input refr-z catalogue
# Written by: A.H. Wright (16-01-2019)
#

#Function prints the help document /*fold*/{{{
.help.print<-function() { 
  cat(paste("\nSOM_DIR.R: script for creating DIR weights for an arbitrary survey\n",
            "Script calling Syntax:\n",
            "Rscript SOM_DIR.R [options] -i <InputReferenceCat> <InputTrainingCat1> <InputTrainingCat2>...\n",
            "Available options:\n",
            "                -p : switch which create plots of the data and SOMs\n",
            "                -q : switch which causes the code to execute quietly (no prompts)\n",
            "                -i : list of input catalogues\n",
            "                -k : list of catalogue keywords. Can be length 2 or more:\n",
            "                     i.e. -k RAkey DECkey FACTORkey1 FACTORkey2 ...\n",
            "                     The first two factors must be those that describe the RA & DEC in the input catalogue.\n",
            "                     The final FACTOR keys, if present, determine the variables that are used for grouping\n",
            "                     sources into systematically variable groups. E.g. it can be just the pointing ID,\n",
            "                     so that randoms are allowed to vary in density between adjacent observations,\n",
            "                     thus better matching the observed source number density variations on sky.\n",
            "                     Or it can be a combination of SEEING, SKYBRIGHTNESS, etc, which are then grouped in a SOM\n",
            "                     and modelled simultaneously by their number density.\n",
            "                     ---> If not set, the default keys are ALPHA_J2000 DELTA_J2000 THELI_NAME\n",
            "                -b : Angular binsize for the randoms grid, in arcmin. Default: 1/3\n",
            "                -l : Path to the LDAC binary file ldactoasc. Default: `which ldactoasc`\n",
            "                -n : number of Randoms to generate. Default: 13E6\n",
            "                -m : multiple factor of input catalogue length, used for number of Randoms to generate.\n",
            "                     i.e. Number of Randoms = M*length(input catalogue). \n",
            "                     ---> NB: OVERRULES THE -n FLAG WHEN PROVIDED!\n",
            "            --seed : seed number to use for Randoms generation. Default: 666\n",
            "--factor.nbins -fn : define the number of hierarchical clusters generated from the FACTORs. Default: 100\n",
            "     --som.dim -sd : define the dimension of the SOM. Default: 55 55 \n",
            "    --som.rate -sr : define the rate of convergence of the SOM. Default: 0.1 0.05 \n",
            "  --som.method -sm : define the method for SOM iteration. Default: 'qbatch' \n",
            "   --som.cores -sc : define the number of cores used during SOM computation. Default: -1 (All) \n",
            "    --som.iter -si : define the number of SOM iterations. Default: 100 \n",
            "    --non.toro -nt : switch which makes the SOM non-toroidal in behaviour\n",
            "            --test : switch which run the pipeline in testing mode; uses a low-res defaults and a thinned input catalogue\n",
            "                -h : Print this help\n",
            "               -hd : Print the default parameter values\n\n\n"))
  quit()
}
.default.print<-function() { 
  cat(paste("\ngenerate_randoms_2D.R: script for creating 2D randoms for an arbitrary survey geometry\n",
            "Script calling Syntax:\n",
            "Rscript SOM_DIR.R [options] -i <InputReferenceCat> <InputTrainingCat1> <InputTrainingCat2>...\n",
            "Available option default values:\n",
            "                  -k : ALPHA_J2000 DELTA_J2000 THELI_NAME \n",
            "                  -l : `which ldactoasc`\n",
            "              --seed : 666\n",
            "  --factor.nbins -fn : 100\n",
            "       --som.dim -sd : 55 55 \n",
            "      --som.rate -sr : 0.05 0.01 \n",
            "    --som.method -sm : 'qbatch' \n",
            "     --som.cores -sc : -1 \n",
            "      --som.iter -si : 100 \n",
            "              --test : modifies these default parameters to:\n",
            "                       --> -sd : 12x12 \n",
            "                       --> -si : 20 \n",
            "                       -->  -b : 4 \n",
            "                       -->  -m : 1 \n",
            "                       --> -sa : 4 \n",
            "NB: items not shown here have no defaults because they behave as switches or are have superseding effects (-m).\n\n"))
  quit()
}
#/*fend*/}}}

#Read the command line options /*fold*/{{{
inputs<-commandArgs(TRUE)
if (length(inputs)==0) { .help.print() } 
#/*fend*/}}}

#Start the timer and load the required Packages /*fold*/ {{{
start.timer<-proc.time()[3]
suppressWarnings(suppressPackageStartupMessages(require(kohonen)))
suppressWarnings(suppressPackageStartupMessages(require(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(require(helpRfuncs)))
#suppressWarnings(suppressPackageStartupMessages(require(LAMBDAR)))
suppressWarnings(suppressPackageStartupMessages(require(data.table)))
suppressWarnings(suppressPackageStartupMessages(require(doParallel)))
suppressWarnings(suppressPackageStartupMessages(require(KernSmooth)))
suppressWarnings(suppressPackageStartupMessages(require(itertools)))
#/*fend*/}}}

#Functions /*fold*/ {{{
#Function to scale a colour palette to desired quantiles /*fold*/ {{{
scale.palette<-function(X,probs=c(0.1,0.9),palette=terrain.colors,n=1e3) { 
  cells<-seq(min(X,na.rm=T),max(X,na.rm=T),length=n+1)
  lims<-quantile(X,probs=probs,na.rm=TRUE)
  palette.cols<-palette(n)
  new.col<-rep(NA,n+1)
  new.col[cells<=min(lims)]<-palette.cols[1]
  new.col[cells>=max(lims)]<-palette.cols[n]
  new.col[cells>min(lims)&cells<max(lims)]<-palette(length(which(cells>min(lims)&cells<max(lims))))
  scale.palette<-colorRampPalette(new.col)
  return=scale.palette
}
#/*fend*/}}}

#Create the function to plot the Density Bar /*fold*/ {{{
.density.bar<-function(palette,z,zrange=range(z,na.rm=na.rm),nBarLines=1e3,bw=diff(zrange)/25/sqrt(12),kernel='rect',col='black',lwd=2,na.rm=TRUE) { 
  #Get the plot limits
  #usr<-par('usr')
  #Create the bar-lines
  zbar<-seq(zrange[1],zrange[2],length=nBarLines)
  #Get the colours
  barcol<-palette(nBarLines)
  #Create the density Plot
  dens<-density(z,bw=bw,kernel=kernel,from=zrange[1],to=zrange[2],na.rm=na.rm)
  #Normalise to unity
  dens$y<-dens$y/max(dens$y)
  #Setup the plot
  plot(NA,axes=F,xlim=zrange,xlab='',ylab='',type='n',ylim=c(0,1.2),xaxs='i')
  #Draw the bar
  arrows(zbar,rep(-0.3,nBarLines),zbar,rep(1.5,nBarLines),length=0,col=barcol) 
  #Draw the density
  lines(dens,lwd=lwd,col=col)
  #Draw the margins
  magaxis(side=c(1,3),labels=c(F,F),majorn=5,minorn=5,mtline=0.2,frame.plot=TRUE)
  mtext(side=2,line=0.2,text='Dens')
  #Draw the axis labels
  sci.tick = maglab(zrange, n = 5, log = F, exptext = F)
  mgp<-par("mgp")
  mgp[2]<-0
  axis(side=1,at=sci.tick$labat,label=sci.tick$exp,tick=FALSE,mgp=mgp)
}                                                                
#/*fend*/}}}

#Function sets a new subplot within the current plot area /*fold*/ {{{
#at the fractional x and y positions frac.fig=c(xlo,xhi,ylo,yhi)
.set.frac.to.fig<-function(frac.fig,cur.fig=par('fig'),new=TRUE,bg=col2alpha('white',0.7),bg.expand=c(0.06,0.01,0.05,0.01)) { 
  if (length(frac.fig)!=4) { stop("frac.fig must be length 4") }
  if (length(bg.expand)!=4) { bg.expand<-rep(bg.expand[1],4) }
  #Define the new fig values
  new.fig<-c(cur.fig[1]+frac.fig[1]*diff(cur.fig[1:2]),
             cur.fig[1]+frac.fig[2]*diff(cur.fig[1:2]),
             cur.fig[3]+frac.fig[3]*diff(cur.fig[3:4]),
             cur.fig[3]+frac.fig[4]*diff(cur.fig[3:4]))
  #Add the background
  if (!is.na(bg)&!is.null(bg)) { 
    usr<-par('usr')
    rect(xleft  =usr[1]+(frac.fig[1]-bg.expand[1])*diff(usr[1:2]),
         xright =usr[1]+(frac.fig[2]+bg.expand[2])*diff(usr[1:2]),
         ybottom=usr[3]+(frac.fig[3]-bg.expand[3])*diff(usr[3:4]),
         ytop   =usr[3]+(frac.fig[4]+bg.expand[4])*diff(usr[3:4]),
         col=bg,border=NA)
  }
  #Set the new fig value
  par(fig=new.fig,new=new)
  return=par('fig')
}
#/*fend*/}}}

#Function sets the mainplot panel for a new figure assuming a /*fold*/ {{{
#matrix array of positions (layout), and the plot panel to start 
#(panelnum). 
.mainplot<-function(layout,panelnum,new=TRUE,flip.yaxis=TRUE) { 
  #Generate the layout
  if (flip.yaxis) { 
    layout<-layout[nrow(layout):1,,drop=FALSE]
  }
  xy=which(layout==panelnum,arr.ind=TRUE)
  if (length(xy)==0) { 
    stop('Panel number does not exist in layout!')
  } else if (length(xy)==2) { 
    new.fig<-c((1/ncol(layout))*(xy[2]-1),
               (1/ncol(layout))*(xy[2]),
               (1/nrow(layout))*(xy[1]-1),
               (1/nrow(layout))*(xy[1]))
  } else {
    new.fig<-c((1/ncol(layout))*min(xy[,2]-1),
               (1/ncol(layout))*max(xy[,2]),
               (1/nrow(layout))*min(xy[,1]-1),
               (1/nrow(layout))*max(xy[,1]))
  }
  par(fig=new.fig,new=new)
  return=par('fig')
}
#/*fend*/}}}

#Function plots the segments for each hierarchical cluster /*fold*/ {{{
plot.segments<-function(codes,col.segments=terrain.colors(ncol(codes)+1),bgcol='transparent',add=FALSE,...) { 
  mar<-par(mar=c(0,0,0,0),oma=c(4,4,4,4))
  ret<-stars(codes,draw.segments=TRUE,label=NULL,len=0,add=add,col.segments=NA,lty=0)
  current.plot <- par("mfg")
  plot.width <- par("usr")[1:2]
  plot.height <- par("usr")[3:4]
  cex <- 1 # First see if the legend fits
  leg.result <- legend('bottom',legend = colnames(codes),cex=cex,
                       plot=FALSE,ncol = 2,fill = terrain.colors(ncol(codes)))
  while (leg.result$rect$w > diff(plot.width)) {
    cex <- cex*0.9 # if too large, decrease text size
    leg.result <- legend('bottom',legend = colnames(codes),cex=cex,
                         plot=FALSE,ncol = 2,fill = terrain.colors(ncol(codes)))
  } # until it fits!
  leg.result <- legend('bottom',legend = colnames(codes),cex=cex,
                       plot=FALSE,ncol = 2,fill = terrain.colors(ncol(codes)),...)
  if (!add) { 
    par(mfg = current.plot)
    plot(ret,type='n',axes=FALSE,ylim=rev(c(max(plot.height) + min(plot.height)+1.5*abs(diff(plot.height))/(max(ret[,2]-0.5)),-1.2*leg.result$rect$h)),xlim=plot.width,xaxs='i',yaxs='i',add=add)
  }
  ret<-stars(codes,draw.segments=TRUE,label=1:nrow(codes),col.segments=col.segments,len=0.8,add=TRUE)
  if (any(bgcol!='transparent')) { 
    symbols(ret[, 1], ret[, 2],squares = rep(2, nrow(ret)), inches = FALSE, add = TRUE, fg = NULL, bg = bgcol)
    ret<-stars(codes,draw.segments=TRUE,label=1:nrow(codes),col.segments=col.segments,len=0.8,add=TRUE)
  }
  #if (length(bgcol)==1) { 
    legend('bottom',legend = colnames(codes),cex=cex,
         plot=TRUE,ncol = 2,fill = terrain.colors(ncol(codes)),...)
  #} else { 
  #  legend('bottomright',legend = colnames(codes),cex=cex,
  #       plot=TRUE,ncol = 2,fill = terrain.colors(ncol(codes)),...)
  #} 
  par(mar=mar)
}
#/*fend*/}}}

#Print Rounding functions /*fold*/ {{{
fsignif<-function(x,digits) gsub('-','\U2212',sapply(signif(x,digits), sprintf, fmt=paste0("%#.",digits,"g")))
fround<-function(x,digits) gsub('-','\U2212',gsub('-0.00','0.00',(sapply(round(x,digits), sprintf, fmt=paste0("%#.",digits,"f")))))
#/*fend*/}}}
#Show bandwidth function /*fold*/ {{{
showbw<-function(dens,kernel=NULL,loc='topleft',scale=0.2,inset=c(0.1,0.1),cex=1,col='black',type='s') {
  #Get location parameters
  usercoord = par()$usr
  xlogcheck = FALSE
  ylogcheck = FALSE
  if (par()$xlog) {
    par(xlog = FALSE)
    par(usr = c(log10(par()$usr[1:2]), par()$usr[3:4]))
    xlogcheck = TRUE
  }
  if (par()$ylog) {
    par(ylog = FALSE)
    par(usr = c(par()$usr[1:2], log10(par()$usr[3:4])))
    ylogcheck = TRUE
  }
  if (length(inset)==1) {
    inset<-c(inset,inset)
  }
  xlo = usercoord[1]
  xhi = usercoord[2]
  ylo = usercoord[3]
  yhi = usercoord[4]
  xdiff = xhi - xlo
  ydiff = yhi - ylo
  xl = xlo + xdiff/2 - xdiff * scale/2
  yb = ylo + ydiff/2 - ydiff * scale/2
  xr = xlo + xdiff/2 + xdiff * scale/2
  yt = ylo + ydiff/2 + ydiff * scale/2
  if (grepl("bottom",loc)) {
    yb = ylo + ydiff * inset[2]
    yt = ylo + ydiff * inset[2] + ydiff * scale
  }
  if (grepl("top",loc)) {
    yb = yhi - ydiff * inset[2] - ydiff * scale
    yt = yhi - ydiff * inset[2]
  }
  if (grepl("left",loc)) {
    xl = xlo + xdiff * inset[1]
    xr = xlo + xdiff * inset[1] + xdiff * scale
  }
  if (grepl("right",loc)) {
    xl = xhi - xdiff * inset[1] - xdiff * scale
    xr = xhi - xdiff * inset[1]
  }
  dy<-yt-yb
  dx<-xr-xl
  #Get density information from call
  if (is.null(kernel)) {
    kernel<-dens$call$kernel
    if (is.null(kernel)) { kernel<-'gaussian' }
  }
  #Get number information from call
  dens.n<-dens$n
  if (is.null(dens.n)) { dens.n<-512 }
  bw<-dens$bw
  #Create Kernel
  kern<-density(rep(0,10),from=-xdiff/2,to=xdiff/2,bw=bw,kernel=kernel,na.rm=TRUE,n=dens.n)
  kern.sum<-sum(kern$y)
  kern$y<-kern$y/max(kern$y)
  kern$y[rev(which(zapsmall(kern$y[1:which.max(kern$y)])==0))[-3:-1]]<-NA
  kern$y[which.max(kern$y):length(kern$y)][which(zapsmall(kern$y[which.max(kern$y):length(kern$y)])==0)[-3:-1]]<-NA
  #ind<-range(which(!is.na(kern$y)))
  #kern$y<-kern$y[ind[1]:ind[2]]
  #kern$x<-kern$x[ind[1]:ind[2]]
  kern$y<-dy*(kern$y*0.8 + 0.1) + yb
  kern$x<-kern$x+xl+dx/2
  text(xl+dx/2,max(kern$y,na.rm=T),lab="Density Kernel",cex=cex,pos=3,col=col)
  text(xl+dx/2,min(kern$y,na.rm=T),lab=bquote(paste("log"[10],"(bw) =",.(fsignif(log10(kern$bw),digits=2)))),cex=cex,pos=1,col=col)
  lines(kern,col=col,type=type)
  par(xlog = xlogcheck)
  par(ylog = ylogcheck)
  par(usr = usercoord)
  return=kern.sum
}
#/*fend*/}}}

#The colour palette /*fold*/ {{{
BlRd<-function(n,alpha=1){rainbow(n,end=2/3,alpha=alpha)[n:1]}
#/*fend*/}}}
#/*fend*/}}}

#Define the colour palette and figure placement/*fold*/ {{{
frac.fig<-c(0.1,0.9,0.1,0.2)
BlBuRd<-colorRampPalette(c('black',rev(brewer.pal(10,"RdBu"))))
WtBuRd<-colorRampPalette(c('white',rev(brewer.pal(10,"RdBu")[-(5)])))
#/*fend*/}}}

#Read the options /*fold*/ {{{
#Default parameter values /*fold*/ {{{
refr.truth<-only.som<-force<-sparse.som<-reuse<-useMult<-quiet<-FALSE
optimise.HCs<-do.zcalib<-short.write<-refr.flag<-train.flag<-FALSE
optimize.z.threshold<-0.01
loop.start<-1
plot<-0
sparse.min.density<-50
som.data.file<-sparse.var<-NULL
seed<-666
res<-200
min.gal.per.core<-1000
maxNAfrac=1
data.threshold<-c(0,80)
data.missing<--99
count.variable.r<-count.variable.t<-''
ldac.options.1<-ldac.options.2<-''
addstr<-''
output.path<-'./'
keys<-"-k MAG_GAAP_u MAG_GAAP_g MAG_GAAP_r MAG_GAAP_i MAG_GAAP_Z MAG_GAAP_Y MAG_GAAP_J MAG_GAAP_H MAG_GAAP_Ks"
zr.label<-"z_spec"
zt.label<-"Z_B"
factor.label<-paste0('MAG_GAAP_',c('u','g','r','i','Z','Y','J','H','Ks'))
factor.nbins<-100
som.dim<-c(55,55)
som.topo<-'hexagonal'
som.toroidal<-TRUE
som.iter<-100
som.method<-"pbatch"
som.cores<-12
som.rate<-c(0.05,0.01)
ldactoasc<-""# system('which ldactoasc',intern=TRUE)
train.catalogues<-NULL
testing<-FALSE
do.QC<-TRUE
#/*fend*/}}}
#Loop through the command arguments /*fold*/ {{{
while (length(inputs)!=0) {
  #Check the options syntax /*fold*/ {{{
  while (length(inputs)!=0 && inputs[1]=='') { inputs<-inputs[-1] }  
  if (!grepl('^-',inputs[1])) {
    print(inputs)
    stop(paste("Incorrect options provided!",
               "Check the lengths for each option!\n",
               "Only -i and -k parameters can have more than 1 item"))
  }
  #Check for test variable /*fold*/ {{{
  if (any(inputs=='--test')) { 
    inputs<-c(inputs[which(inputs=='--test')],inputs[which(inputs!='--test')])
  }
  #/*fend*/}}}
  #/*fend*/}}}
  if (inputs[1]=='-p') { 
    #Create the plots /*fold*/ {{{
    plot<-1
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-pp') {
    #Do not create the plots /*fold*/ {{{
    plot<-2
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-np') {
    #Do not create the plots /*fold*/ {{{
    plot<-0
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-as') {
    #read the addstring  /*fold*/ {{{
    inputs<-inputs[-1]
    addstr<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-r') {
    #Read the input reference catalogue(s) /*fold*/ {{{
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) { 
      refr.catalogues<-inputs[1:(which(grepl('^-',inputs))[1]-1)]
      inputs<-inputs[-(1:(which(grepl('^-',inputs))[1]-1))]
    } else { 
      refr.catalogues<-inputs
      inputs<-NULL
    } 
    #/*fold*/}}}
  } else if (inputs[1]=='-t') {
    #Read the input training catalogue(s) /*fold*/ {{{
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) { 
      train.catalogues<-inputs[1:(which(grepl('^-',inputs))[1]-1)]
      inputs<-inputs[-(1:(which(grepl('^-',inputs))[1]-1))]
    } else { 
      train.catalogues<-inputs
      inputs<-NULL
    } 
    #/*fend*/}}}
  } else if (inputs[1]=='-cr') {
    #Define the count variable for the reference cat /*fold*/ {{{
    inputs<-inputs[-1]
    count.variable.r<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-ct') {
    #Define the count variable for the training cat(s) /*fold*/ {{{
    inputs<-inputs[-1]
    count.variable.t<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--noqc') {
    #Read the quality control expressions /*fold*/ {{{
    do.QC<-FALSE
    inputs<-inputs[-1]
    warning("Not performing QC of clusters!") 
    #/*fend*/}}}
  } else if (inputs[1]=='--zt.calib') {
    #Read the quality control expressions /*fold*/ {{{
    inputs<-inputs[-1]
    zcalib.expr<-inputs[1]
    do.zcalib<-TRUE
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-qc') {
    #Read the quality control expressions /*fold*/ {{{
    inputs<-inputs[-1]
    qc.expr<-inputs[1]
    inputs<-inputs[-1]
    warning("Modifying QC expression away from Fiducial!") 
    #/*fend*/}}}
  } else if (inputs[1]=='-k') {
    #Read the catalogue keywords /*fold*/ {{{
    if (any(grepl('^-',inputs[-1]))) {
      key.ids<-1:(which(grepl('^-',inputs[-1]))[1])
    } else { 
      key.ids<-1:length(inputs)
    }
    keys<-inputs[key.ids[c(-1)]]
    inputs<-inputs[-(key.ids)]
    if (length(key.ids)<4) { 
      stop("There are no columns to group the data on!")
    } else if (length(key.ids)>=4) { 
      suppressWarnings(suppressPackageStartupMessages(require(kohonen)))
      factor.label<-keys
    } else { 
      stop(paste0("Catalogue Keywords must have length 3 or greater. It is ",length(key.ids),":-> ",keys," <-\n",
                  "i.e.:\n",
                  "Rscript SOM_DIR.R -k MAG_GAAP_u MAG_GAAP_g MAG_GAAP_r MAG_GAAP_i -r <InputReferenceCat> -t <InputTrainingCat1> ...\n",
                  "Rscript SOM_DIR.R -k MAG_GAAP_u-MAG_GAAP_g MAG_GAAP_r-MAG_GAAP_i MAG_GAAP_Z-MAG_GAAP_Y MAG_GAAP_J-MAG_GAAP_H MAG_GAAP_r -r <InputReferenceCat> -t <InputTrainingCat1> ...\n"))
    }
    #Create the full keys string/*fold*/ {{{
    #keys<-paste(keys,collapse=' ')
    keys<-paste('-k',paste(unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",keys),' '))),collapse=' '))
    #/*fend*/}}}
    #/*fend*/}}}
  } else if (inputs[1]=='-of') {
    #Read the output file name(s) /*fold*/ {{{
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) { 
      output.file<-inputs[1:(which(grepl('^-',inputs))[1]-1)]
      inputs<-inputs[-(1:(which(grepl('^-',inputs))[1]-1))]
    } else { 
      output.file<-inputs
      inputs<-NULL
    } 
    #/*fold*/}}}
  } else if (inputs[1]=='-o') {
    #Define the output directory  /*fold*/ {{{
    inputs<-inputs[-1]
    output.path<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-ls') {
    #Define the loop.start counter /*fold*/ {{{
    inputs<-inputs[-1]
    loop.start<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-l') {
    #Define the ldactoasc binary  /*fold*/ {{{
    inputs<-inputs[-1]
    ldactoasc<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-lor') {
    #Define the ldacoptions for reference catalogue /*fold*/ {{{
    inputs<-inputs[-1]
    ldac.options.1<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-lot') {
    #Define the ldacoptions for training catalogue(s)  /*fold*/ {{{
    inputs<-inputs[-1]
    ldac.options.2<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-nt'|inputs[1]=='--nontoroidal') {
    #Define whether the SOM is toroidal /*fold*/ {{{
    inputs<-inputs[-1]
    som.toroidal<-FALSE
    #/*fend*/}}}
  } else if (inputs[1]=='--zt.label') {
    #Define the z_spec label in the training catalogue /*fold*/ {{{
    inputs<-inputs[-1]
    zt.label<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--zr.label') {
    #Define the z_phot label in the training & reference catalogues /*fold*/ {{{
    inputs<-inputs[-1]
    zr.label<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--topo') {
    #Define whether the SOM is toroidal /*fold*/ {{{
    inputs<-inputs[-1]
    som.topo<-inputs[1]
    inputs<-inputs[-1]
    if (!som.topo%in%c("rectangular","hexagonal")) { 
      stop(paste0('SOM Topology must be "rectangular" or "hexagonal", not: ',som.topo))
    }
    #/*fend*/}}}
  } else if (inputs[1]=='--train.flag') {
    #Define whether the training weights should be output as a 0/1 flag /*fold*/ {{{
    inputs<-inputs[-1]
    train.flag<-TRUE
    #/*fend*/}}}
  } else if (inputs[1]=='--refr.flag') {
    #Define whether the reference weights should be output as a 0/1 flag /*fold*/ {{{
    inputs<-inputs[-1]
    refr.flag<-TRUE
    #/*fend*/}}}
  } else if (inputs[1]=='--short.write') {
    #Define whether the to write short catalogues /*fold*/ {{{
    inputs<-inputs[-1]
    short.write<-TRUE
    #/*fend*/}}}
  } else if (inputs[1]=='--refr.truth') {
    #Define whether we want to compute the true z of the reference sample /*fold*/ {{{
    inputs<-inputs[-1]
    refr.truth<-TRUE
    #/*fend*/}}}
  } else if (inputs[1]=='--optimise') {
    #Define whether we want to optimise the number of HCs for representation /*fold*/ {{{
    inputs<-inputs[-1]
    optimise.HCs<-TRUE
    #/*fend*/}}}
  } else if (inputs[1]=='--toroidal') {
    #Define whether the SOM is toroidal /*fold*/ {{{
    inputs<-inputs[-1]
    som.toroidal<-TRUE
    #/*fend*/}}}
  } else if (inputs[1]=='--factor.nbins'|inputs[1]=='-fn') {
    #Define the number of factor bins /*fold*/ {{{
    inputs<-inputs[-1]
    factor.nbins<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--sparse.var') {
    #/*fold*/ {{{
    inputs<-inputs[-1]
    sparse.var<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--sparse.som') {
    #Load an already calculated som from file /*fold*/ {{{
    inputs<-inputs[-1]
    sparse.frac<-as.numeric(inputs[1])
    if (sparse.frac>=1) { 
      stop("Sparse sampling fraction must be less than 1!") 
    } else if (sparse.frac <=0) { 
      stop("Sparse sampling fraction must be greater than 0!") 
    }
    inputs<-inputs[-1]
    sparse.som<-TRUE
    #/*fend*/}}}
  } else if (inputs[1]=='--old.som') {
    #Load an already calculated som from file /*fold*/ {{{
    #inputs<-inputs[-1]
    #som.data.file<-inputs[1]
    #inputs<-inputs[-1]
    #reuse<-TRUE
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) { 
      som.data.file<-inputs[1:(which(grepl('^-',inputs))[1]-1)]
      inputs<-inputs[-(1:(which(grepl('^-',inputs))[1]-1))]
    } else { 
      som.data.file<-inputs
      inputs<-NULL
    } 
    reuse<-TRUE
    #/*fold*/}}}
  } else if (inputs[1]=='--som.dim'|inputs[1]=='-sd') {
    #Define the SOM dimension /*fold*/ {{{
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) {
      dim.ids<-1:(which(grepl('^-',inputs))[1]-1)
    } else { 
      dim.ids<-1:length(inputs)
    }
    if (length(dim.ids)!=2) { stop("som.dim must be of length 2; i.e. --som.dim 55 55") } 
    som.dim<-as.numeric(inputs[dim.ids])
    inputs<-inputs[-dim.ids]
    #/*fend*/}}}
  } else if (inputs[1]=='--som.rate'|inputs[1]=='-sr') {
    #Define the SOM rate /*fold*/ {{{
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) {
      rate.ids<-1:(which(grepl('^-',inputs))[1]-1)
    } else { 
      rate.ids<-1:length(inputs)
    }
    if (length(rate.ids)!=2) { stop("som.rate must be of length 2; i.e. --som.rate 0.05 0.01") } 
    som.rate<-as.numeric(inputs[rate.ids])
    inputs<-inputs[-rate.ids]
    #/*fend*/}}}
  } else if (inputs[1]=='--som.method'|inputs[1]=='-sm') {
    #Define the SOM method /*fold*/ {{{
    inputs<-inputs[-1]
    som.method<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--som.cores'|inputs[1]=='-sc') {
    #Define the number of SOM cores /*fold*/ {{{
    inputs<-inputs[-1]
    som.cores<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--som.na.max'|inputs[1]=='-na') {
    #Define the number of SOM iterations /*fold*/ {{{
    inputs<-inputs[-1]
    maxNAfrac<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--som.iter'|inputs[1]=='-si') {
    #Define the number of SOM iterations /*fold*/ {{{
    inputs<-inputs[-1]
    som.iter<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--data.missing') {
    #Define the value for missing data /*fold*/ {{{
    inputs<-inputs[-1]
    data.missing<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--data.threshold') {
    #Define the SOM rate /*fold*/ {{{
    inputs<-inputs[-1]
    data.threshold<-as.numeric(inputs[1:2])
    inputs<-inputs[-2:-1]
    if (any(is.na(data.threshold))) { stop("data.threshold parameters are NA. For no thresholding, set --data.threshold -Inf Inf") }
    #/*fend*/}}}
  } else if (inputs[1]=='--seed') {
    #Define the seed for SOM generation /*fold*/ {{{
    inputs<-inputs[-1]
    seed<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--force') { 
    #Force catalogue creation /*fold*/ {{{
    force<-TRUE
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--only.som') { 
    #Run in testing mode /*fold*/ {{{
    only.som<-TRUE
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--test') { 
    #Run in testing mode /*fold*/ {{{
    testing<-TRUE
    inputs<-inputs[-1]
    #Set the testing defaults /*fold*/ {{{
    som.iter<-20
    som.dim<-c(12,12)
    #/*fend*/}}}
    #/*fend*/}}}
  } else if (inputs[1]=='-q') {
    #Run Quietly /*fold*/ {{{
    quiet<-TRUE
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-hd') {
    #Print the default values help function/*fold*/ {{{
    .default.print()
    return(NULL)
    #/*fend*/}}}
  } else if (inputs[1]=='-h') {
    #Print the help function/*fold*/ {{{
    .help.print()
    return(NULL)
    #/*fend*/}}}
  } else {
    stop(paste("Unknown option",inputs[1]))
  }
}
#/*fend*/}}}
#/*fend*/}}}

#Setup the QC expression {{{ 
if (!exists("qc.expr") & do.QC) {
  if (count.variable.t=="" & count.variable.r=="") { 
    #Straight mean for train.cat and refr.cat
    qc.expr<-"abs(mean(train.cat[[zt.label]])-mean(refr.cat[[zr.label]]))<=matrixStats::rowMaxs(cbind(5*mad(full.train.cat[[zt.label]]-full.train.cat[[zr.label]]),0.4),na.rm=T)" 
  } else if (count.variable.r=="") { 
    #Weighted mean for train.cat & Straight mean for refr.cat 
    qc.expr<-"abs(weighted.mean(train.cat[[zt.label]],train.cat[[count.variable.t]])-mean(refr.cat[[zr.label]]))<=matrixStats::rowMaxs(cbind(5*mad(full.train.cat[[zt.label]]-full.train.cat[[zr.label]]),0.4),na.rm=T)" 
  } else if (count.variable.t=="") { 
    #Straight mean for train.cat & Weighted mean for refr.cat 
    qc.expr<-"abs(mean(train.cat[[zt.label]])-weighted.mean(refr.cat[[zr.label]],refr.cat[[count.variable.r]]))<=matrixStats::rowMaxs(cbind(5*mad(full.train.cat[[zt.label]]-full.train.cat[[zr.label]]),0.4),na.rm=T)" 
  } else {
    #Weighted mean for train.cat & Weighted mean for refr.cat 
    qc.expr<-"abs(weighted.mean(train.cat[[zt.label]],train.cat[[count.variable.t]])-weighted.mean(refr.cat[[zr.label]],refr.cat[[count.variable.r]]))<=matrixStats::rowMaxs(cbind(5*mad(full.train.cat[[zt.label]]-full.train.cat[[zr.label]]),0.4),na.rm=T)" 
  }
}
#}}}

#Check for the output path /*fold*/ {{{
if (!dir.exists(output.path)) { 
  dir.create(output.path,recursive=TRUE)
}
#/*fend*/}}}

#Check for the input catalogues /*fold*/ {{{
n.catalogues<-max(length(train.catalogues),length(refr.catalogues),length(som.data.file))
if (n.catalogues==0) { 
  stop("No catalogues provided! Call Syntax:\nRscript SOM_DIR.R [options] -i <InputReferenceCat> <InputTrainingCat1> ...") 
}
#Make sure catalogue lists have matching length {{{
if (length(train.catalogues)!=n.catalogues){ 
  if (length(train.catalogues)==1) { 
    train.catalogues<-rep(train.catalogues,n.catalogues)
  } else if (n.catalogues%%length(train.catalogues)==0) { 
    train.catalogues<-rep(train.catalogues,n.catalogues/length(train.catalogues))
  } else { 
    stop("Input training catalogue list is neither length 1, length(reference catalogue list), nor a factor of length(reference catalogue list)")
  }
}
if (length(refr.catalogues)!=n.catalogues){ 
  if (length(refr.catalogues)==1) { 
    refr.catalogues<-rep(refr.catalogues,n.catalogues)
  } else if (n.catalogues%%length(refr.catalogues)==0) { 
    refr.catalogues<-rep(refr.catalogues,n.catalogues/length(refr.catalogues))
  } else { 
    stop("Input reference catalogue list is neither length 1, length(reference catalogue list), nor a factor of length(reference catalogue list)")
  }
}
#}}}
#Make sure input SOM lists have matching length {{{
if (reuse) { 
  if (length(som.data.file)!=n.catalogues){ 
    if (length(som.data.file)==1) { 
      som.data.file<-rep(som.data.file,n.catalogues)
    } else if (n.catalogues%%length(som.data.file)==0) { 
      som.data.file<-rep(som.data.file,n.catalogues/length(som.data.file))
    } else { 
      stop("Input reference catalogue list is neither length 1, length(reference catalogue list), nor a factor of length(reference catalogue list)")
    }
  }
}
#}}}
#/*fend*/}}}

#Prompt and start /*fold*/ {{{
if (!quiet) { 
  cat("Creating DIR weights with input training catalogue(s):\n ") 
  if (n.catalogues > 5) { 
    cat(paste0('  --> ',train.catalogues[1:3],'\n'))
    cat(paste0('  ....\n'))
    cat(paste0('  --> ',rev(train.catalogues)[3:1],'\n'))
  } else {
    cat(paste0('  --> ',train.catalogues,'\n'))
  }
  cat("And corresponding input reference catalogue(s):\n ") 
  if (n.catalogues > 5) { 
    cat(paste0('  --> ',refr.catalogues[1:3],'\n'))
    cat(paste0('  ....\n'))
    cat(paste0('  --> ',rev(refr.catalogues)[3:1],'\n'))
  } else {
    cat(paste0('  --> ',refr.catalogues,'\n'))
  }
}
#/*fend*/}}}

#Set the seed for randomisation /*fold*/ {{{
set.seed(seed)
#/*fend*/}}}

#Initialise loop counters /*fold*/ {{{
io.clock<-0
loop.num<-loop.start-1 #this is just here for name convenience when the pipe errors
loop.length<-n.catalogues+loop.start-1
#/*fend*/}}}

#Check for the output catalogues /*fold*/ {{{
if (!exists('output.file')) { 
  #Training filenames 
  train.catnam=vecsplit(train.catalogues,'/',-1)
  output.file<-train.catnam
  #Define the filename addition /*fold*/ {{{
  addstr<-paste0('_SOM_',paste(factor.label,collapse='_'))
  
  if (nchar(addstr>25)) { 
    print.length<-length(which(cumsum(nchar(gsub("[aeiou_]","",factor.label,perl=T,ignore.case=T)))<25))
    addstr<-paste0('_SOM_',paste(gsub("[aeiou_]","",factor.label[1:print.length],
                                      perl=T,ignore.case=T),collapse='_'),
                   '_pl_',length(factor.label)-print.length)
    addstr<-rep(addstr,n.catalogues)
  }
  #/*fend*/}}}
} else if (length(output.file)!=n.catalogues) { 
  addstr<-paste0("_",1:n.catalogues)
  output.file<-rep(output.file,n.catalogues)
}
if (length(addstr)!=n.catalogues) { addstr<-rep(addstr,n.catalogues) } 
output.ending<-vecsplit(output.file,'.',-1,fixed=T)
#/*fend*/}}}

for (catpath.count in 1:n.catalogues) { 
#Define the training and reference catalogue paths {{{
train.catpath<-train.catalogues[catpath.count]
refr.catpath<-refr.catalogues[catpath.count]
if (reuse) {
  som.datfile<-som.data.file[catpath.count]
}
#}}}
#Loop through catalogues /*fold*/ {{{
loop.num<-loop.num+1
#Get the catalogue name and ending from the catalogue path given #/*fold*/ {{{
if (!quiet) { cat(paste("Working on Catalogues:\n    ",train.catpath,"\n    ",refr.catpath,"\n")) }
train.catnam<-rev(strsplit(train.catpath,'/')[[1]])[1]
train.ending<-rev(strsplit(train.catnam,'.',fixed=TRUE)[[1]])[1]
##Check for the output catalogues /*fold*/ {{{
#if (!exists('output.file')) { 
#  output.file<-train.catnam
#  #Define the filename addition /*fold*/ {{{
#  addstr<-paste0('_SOM_',paste(factor.label,collapse='_'))
#  
#  if (nchar(addstr>25)) { 
#    print.length<-length(which(cumsum(nchar(gsub("[aeiou_]","",factor.label,perl=T,ignore.case=T)))<25))
#    addstr<-paste0('_SOM_',paste(gsub("[aeiou_]","",factor.label[1:print.length],
#                                      perl=T,ignore.case=T),collapse='_'),
#                   '_pl_',length(factor.label)-print.length)
#  }
#  #/*fend*/}}}
#} else if (length(output.file)!=n.catalogues) { 
#  addstr<-paste0("_",loop.num)
#}
#output.ending<-rev(strsplit(output.file,'.',fixed=TRUE)[[1]])[1]
##/*fend*/}}}
#/*fend*/}}}

#Check if the output catalogue already exists /*fold*/ {{{
if (grepl('.cat',output.file[loop.num],fixed=T)) {
  #Output an FITS catalogue /*fold*/ {{{
  file =paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_DIRsom',addstr[loop.num],'.fits'),output.file[loop.num],fixed=TRUE))
  file2=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_refr_DIRsom',addstr[loop.num],'.fits'),output.file[loop.num],fixed=TRUE))
  #/*fend*/}}}
} else if (grepl('.fits',output.file[loop.num],fixed=T)) {
  #Output an fits catalogue /*fold*/ {{{
  file =paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_DIRsom',addstr[loop.num],'.fits'),output.file[loop.num],fixed=TRUE))
  file2=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_refr_DIRsom',addstr[loop.num],'.fits'),output.file[loop.num],fixed=TRUE))
  #/*fend*/}}}
} else if (grepl('.Rdata',output.file[loop.num],fixed=T)) {
  #Output an Rdata catalogue /*fold*/ {{{
  file =paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_DIRsom',addstr[loop.num],'.Rdata'),output.file[loop.num],fixed=TRUE))
  file2=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_refr_DIRsom',addstr[loop.num],'.Rdata'),output.file[loop.num],fixed=TRUE))
  #/*fend*/}}}
} else { 
  #Output a CSV catalogue /*fold*/ {{{
  file =paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_DIRsom',addstr[loop.num],'.csv'),output.file[loop.num],fixed=TRUE))
  file2=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_refr_DIRsom',addstr[loop.num],'.csv'),output.file[loop.num],fixed=TRUE))
  #/*fend*/}}}
}
somfile<-paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0(addstr[loop.num],'_SOMdata.Rdata'),output.file[loop.num],fixed=TRUE))
if (file.exists(file) & file.exists(file2) | (only.som & file.exists(somfile))) { 
  if (!force) { 
    cat("  ### Output Catalogues already exist! Not using [--force], so skipping! ###\n")
    next
  } else { 
    cat("  ### Output Catalogues already exist, but continuing because of --force ###\n")
  } 
}
#/*fend*/}}}

#Read in the input training catalogue /*fold*/ {{{
if (catpath.count!=1 && train.catalogues[catpath.count-1]==train.catpath) { 
  if (!quiet) { 
    cat(paste("  > Skipping Training Catalogue Read (same as previous loop!)")) 
    timer<-proc.time()[3]
  }
} else { 
  if (!quiet) { 
    cat(paste("  > Reading Training Catalogue")) 
    timer<-proc.time()[3]
  }
  if (short.write | only.som) { 
    cols<-unique(vecsplit(factor.label,"[+/-]|\\*",fixed=FALSE))
    if (!only.som) { cols<-c(cols,zt.label,zr.label) }
    if (count.variable.t!="") { cols<-c(cols,count.variable.t) }
    print(cols)
    train.cat<-read.file(train.catpath,cols=cols)
  } else { 
    train.cat<-read.file(train.catpath)
  }
}
if (nrow(train.cat)==0) { 
  #The training catalogue is empty?!
  stop(paste0("Training catalogue was read successfully but has no rows!\n"))
}
#Make sure that the table is a data.table/*fold*/ {{{
if (!is.data.table(train.cat)){ 
  train.cat<-as.data.table(train.cat)
}
train.catnam<-rev(strsplit(train.catpath,'/')[[1]])[1]
train.ending<-rev(strsplit(train.catnam,'.',fixed=TRUE)[[1]])[1]
if (testing) {
  train.cat<-train.cat[runif(nrow(train.cat))<0.2,]
}
#Catalogue length/*fold*/ {{{
train.cat.len<-nrow(train.cat)
if (!quiet) { 
  io.clock<-io.clock+proc.time()[3]-timer
  cat(paste(" - Done",as.time(proc.time()[3]-timer,digits=0),"\n"))
  if (testing) { 
    cat("  # Running in Testing mode! Many parameters degraded to improve speed (--som.dim) #\n")
  }
}
#/*fend*/}}}
#/*fend*/}}}
#/*fend*/}}}

#Read in the input reference catalogue /*fold*/ {{{
if (catpath.count!=1 && refr.catalogues[catpath.count-1]==refr.catpath) { 
  if (!quiet) { 
    cat(paste("  > Skipping Reference Catalogue Read (same as previous loop!)")) 
    timer<-proc.time()[3]
  }
} else { 
  if (only.som) { 
    if (!quiet) { 
      cat(paste("  > Skipping Reference Catalogue Read (only training a SOM!)")) 
      timer<-proc.time()[3]
      #Assign the training cat to the reference catalogue (won't be used meaningfully) 
      refr.cat<-train.cat
    }
  } else { 
    if (!quiet) { 
      cat(paste("  > Reading Reference Catalogue")) 
      timer<-proc.time()[3]
    }
    if (short.write) { 
      cols<-c(unique(vecsplit(factor.label,"[+/-]|\\*",fixed=FALSE)),zr.label)
      if (refr.truth) { cols<-c(cols,zt.label) }
      if (count.variable.r!="") { cols<-c(cols,count.variable.r) }
      print(cols)
      refr.cat<-read.file(refr.catpath,cols=cols)
    } else { 
      refr.cat<-read.file(refr.catpath)
    } 
  }
}
#Make sure that the table is a data.table
if (!is.data.table(refr.cat)){ 
  refr.cat<-as.data.table(refr.cat)
}
if (testing) {
  refr.cat<-refr.cat[runif(nrow(refr.cat))<0.2,]
}
#Catalogue length
refr.cat.len<-nrow(refr.cat)
if (!quiet) { 
  io.clock<-io.clock+proc.time()[3]-timer
  cat(paste(" - Done",as.time(proc.time()[3]-timer,digits=0),"\n"))
  if (testing) { 
    cat("  # Running in Testing mode! Many parameters degraded to improve speed (--som.dim) #\n")
  }
}

#/*fend*/}}}

#Generate the DIR wieights /*fold*/ {{{
#Notify /*fold*/ {{{
if (!quiet) { 
  cat(paste("  > Generating the SOM from Training catalogue "))
  short.timer<-timer<-proc.time()[3]
}
#/*fend*/}}}
#Measure the SOM /*fold*/ {{{
if (length(factor.label)<2) { 
  #There is a insufficient data for SOM /*fold*/ {{{
  stop("SOM Construction requires 2 or more data dimensions!")
  #/*fend*/}}}
} else { 
  #There are multiple factors: Cluster the data with self organising map /*fold*/ {{{
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat(" (multiple numeric factors) {\n")
    cat("    -> constructing scaled data vector")
  }#/*fend*/}}}
  #Check for character columns /*fold*/ {{{
  seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",factor.label),' ')))
  if (any(sapply(train.cat[,seperated.labels,with=F],class)=='character')) { 
    if (!quiet) { cat(" (converting character cols to numeric!)") }
    for (i in sperated.labels[which(sapply(train.cat[1,seperated.labels,with=F],class)=='character')]) { 
      train.cat[,i,with=F]<-as.numeric(train.cat[,i,with=F])
    }
  }
  #/*fend*/}}}
  #Load or Generate the SOM /*fold*/ {{{
  if (reuse && file.exists(som.datfile)) {
    #Use a previously constructed SOM /*fold*/ {{{
    #Notify /*fold*/ {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> loading previously calculated SOM (!!)")
    }#/*fend*/}}}
    #Load the SOM /*fold*/ {{{
    name<-load(som.datfile)
    if (name!="train.som") { 
      train.som<-get(name)
    }
    som.dim<-c(train.som$grid$xdim,train.som$grid$ydim)
    #/*fend*/}}}
    #Notify /*fold*/ {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> passing training set into previous SOM (!!)")
    }#/*fend*/}}}
    #Check for previous SOM versions {{{
    if (length(train.som$whiten.param)==0 & length(train.som$rescale.param)!=0) { 
      warning("Updating SOM file to current style: rescale.param->whiten.param")
      train.som$whiten.param<-train.som$rescale.param
    } else if (length(train.som$whiten.param)==0) { 
      stop("SOM has no whitening parameters")
    }
    #}}}
    #Get the data positions from the trained SOM /*fold*/ {{{
    train.som<-kohparse(som=train.som,data=train.cat,train.expr=factor.label,data.missing=data.missing,data.threshold=data.threshold,
                        n.cores=som.cores,max.na.frac=1,quiet=TRUE)
    #/*fend*/}}}
    #/*fend*/}}}
  } else { 
    #Construct a new SOM /*fold*/ {{{
    #If we wanted to reuse, stop with error /*fold*/ {{{
    if (reuse) { 
      stop("SOM data file does not exist! Cannot use previous SOM!\n",som.datfile)
    }
    #/*fend*/}}}
    #Run the training /*fold*/ {{{
    train.som<-kohtrain(data=train.cat,train.expr=factor.label,som.dim=som.dim,som.topo=som.topo,
                        som.toroidal=som.toroidal,som.iter=som.iter,som.rate=som.rate,n.cores=som.cores,
                        train.sparse=sparse.som,sparse.min.density=sparse.min.density,sparse.var=sparse.var,
                        data.missing=data.missing,data.threshold=data.threshold,quiet=quiet,seed=seed)
    #/*fend*/}}}
    #If we removed training rows because of NAs and/or sparsity, update the prediction /*fold*/ {{{
    if (length(train.som$unit.classif)!=nrow(train.cat)) {
      #Check the maxNA values /*fold*/ {{{
      if (maxNAfrac==1 & !sparse.som) { 
        stop(paste0("We have lost training rows while maxNAfrac==1 and sparse.som==FALSE!\n",
             "There must be a parallelisation error! Try reducing the number of requested threads."))
      } else if (maxNAfrac!=train.som$maxNA.fraction) { 
        stop("BUG: The SOM training used the wrong maxNAfrac?!")
      }
      #/*fend*/}}}
      #Get the data positions from the trained SOM /*fold*/ {{{
      train.som<-kohparse(som=train.som,data=train.cat,train.expr=factor.label,data.missing=data.missing,data.threshold=data.threshold,
                          n.cores=som.cores,max.na.frac=1,quiet=TRUE)
      #/*fend*/}}}
    }
    #/*fend*/}}}
    #Notify /*fold*/ {{{
    if (!quiet) { 
      short.timer<-proc.time()[3]
      cat("\n    -> Saving SOM ")
    }#/*fend*/}}}
    #Save the SOM to file /*fold*/ {{{
    save(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0(addstr[loop.num],'_SOMdata.Rdata'),output.file[loop.num],fixed=TRUE)),train.som)
    #/*fend*/}}}
    #Skip the rest of the measurements?! /*fold*/ {{{
    if (only.som) { 
      #Notify /*fold*/ {{{
      if (!quiet) { 
        cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
        short.timer<-proc.time()[3]
        cat("\n    -> Skipping remaining loop because '--only.som' flag was set!")
        cat("\n       (In the future, you could just use the kohtrain function directly)\n")
        next
      }#/*fend*/}}}
    }
    #/*fend*/}}}
    #/*fend*/}}}
  } 
  #/*fend*/}}}
  #Plot the SOM and additional information /*fold*/ {{{
  if (plot>0) { 
    #Notify /*fold*/ {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> Plotting SOM figures")
    }#/*fend*/}}}
    #Plot the som codes /*fold*/ {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_SOM_codes.png',output.file[loop.num],fixed=TRUE)),height=10*res,width=10*res,res=res)
    plot(train.som, type="codes",shape='straight',codeRendering='lines',border=NA)
    dev.off()
    #/*fend*/}}}
    #Plot the cluster makeup /*fold*/ {{{
    #Plot the changes at each iteration /*fold*/ {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_SOM_changes.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=5*res,res=res)
    plot(train.som, type="changes")
    dev.off()
    #/*fend*/}}}
    if (plot>1) {
      #Plot the various SOM components /*fold*/ {{{
      png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_SOM_props%02d.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=5*res,res=res)
      #Set the margins /*fold*/ {{{
      par(mar=c(0,0,0,0))
      #/*fend*/}}}
      #Get the colour palette /*fold*/ {{{
      BlRd<-function(n,alpha=1){rainbow(n,end=2/3,alpha=alpha)[n:1]}
      #/*fend*/}}}
      #Plot the counts per SOM cell /*fold*/ {{{
      cell.counts<-plot(train.som, type="count",shape='straight',border=NA,heatkeywidth=1,zlog=TRUE)
      #/*fend*/}}}
      #Plot the neighbour distances per SOM cell /*fold*/ {{{
      plot(train.som, type="dist.neighbours",shape='straight',border=NA,heatkeywidth=1)
      #/*fend*/}}}
      #Plot the quality per SOM cell /*fold*/ {{{
      plot(train.som, type="quality",shape='straight',border=NA,heatkeywidth=1)
      #/*fend*/}}}
      #Plot the number density per SOM cell /*fold*/ {{{
      plot(train.som, type="count",shape='straight',border=NA,heatkeywidth=1)
      #/*fend*/}}}
      for (i in 1:length(train.som$codes[[1]][1,])) { 
        #Plot the i^th data property per SOM cell /*fold*/ {{{
        plot(train.som, type = "property", property = i,
             main=paste0(colnames(train.som$data[[1]])[i],'(scaled)'), palette.name=BlRd,shape='straight',border=NA,heatkeywidth=1)
        #/*fend*/}}}
      }
      #Close the plot device /*fold*/ {{{
      dev.off()
      #/*fend*/}}}
      #/*fend*/}}}
    }
    #/*fend*/}}}
  }
  #/*fend*/}}}
  #/*fend*/}}}
} 
#/*fend*/}}}
#Notify /*fold*/ {{{
if (!quiet) { 
  cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
  short.timer<-proc.time()[3]
  cat("\n    -> Clustering the SOM cells")
}#/*fend*/}}}
#Check the number of HCs /*fold*/ {{{
if (factor.nbins>=prod(som.dim)) { 
  if (!quiet) { 
    cat("\n  # WARNING: Number of hierarchical clusters is >= number of SOM cells!           #\n")
    cat("  #          We will not use the hierarchical clustering, just split by SOM cell! #\n")
    cat("\n    -> Contining with clustering")
  }
  factor.nbins<-prod(som.dim)
}
#/*fend*/}}}
#If optimising HCs, save the original factor.nbins {{{
if (optimise.HCs & !exists("factor.nbins.original")) { 
  factor.nbins.original<-factor.nbins
} else if (optimise.HCs) { 
  factor.nbins<-factor.nbins.original
}
#}}}
#Get the training positions from the trained SOM /*fold*/ {{{
train.som<-generate.kohgroups(som=train.som,n.cluster.bins=factor.nbins,
                             train.expr=factor.label,data.missing=data.missing,
                             data.threshold=data.threshold,n.cores=som.cores,max.na.frac=1,quiet=TRUE)
#/*fend*/}}}
#Notify /*fold*/ {{{
if (!quiet) { 
  cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
  short.timer<-proc.time()[3]
  cat("\n    -> Passing reference data vector into trained & clustered SOM")
}#/*fend*/}}}
#Get the data positions from the trained SOM /*fold*/ {{{
refr.som<-generate.kohgroups(som=train.som,n.cluster.bins=factor.nbins,new.data=refr.cat,
                             train.expr=factor.label,data.missing=data.missing,
                             data.threshold=data.threshold,n.cores=som.cores,max.na.frac=1,quiet=TRUE)
#/*fend*/}}}
#Optimise the number of HCs? /*fold*/ {{{
if (optimise.HCs) { 
  #Define the HC.steps {{{
  HC.steps<-ceiling(seq(factor.nbins*0.01,factor.nbins,length=100))
  HC.steps<-HC.steps[which(!duplicated(HC.steps))]
  #}}}
  #Notify {{{
  if (!quiet) { 
    cat(paste(" - Pausing",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat(paste("\n    -> Optimising number of heirarchical clusters"))
    cat(paste("\n    -> constructing",length(HC.steps),"sets of hierarchical clusters"))
  }#}}}
  #Cut the cluster dendrogram at the required cluster numbers {{{
  if (is.null(train.som$hclust)) {
    train.hc.mat = cutree(hclust(dist(train.som$codes[[1]])), as.numeric(HC.steps))
  } else {
    train.hc.mat = cutree(train.som$hclust, as.numeric(HC.steps))
  }
  #}}}
  #Notify {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Dispersing training data into cluster bins")
  }#}}}
  #Optimise by datum or by cell {{{
  optimise.cell<-TRUE
  if (optimise.cell) { 
    #By Cell {{{
    #If we're doing QC, prepare the QC expressions {{{
    if (do.QC) { 
      #Check the qc expression for errors {{{ 
      if (grepl("count.variable.t",qc.expr)&count.variable.t=='') { 
        stop("count.variable.t is not defined but is used in the cluster QC")
      }
      if (grepl("count.variable.r",qc.expr)&count.variable.r=='') { 
        stop("count.variable.r is not defined but is used in the cluster QC")
      }
      if (count.variable.t!=''&!count.variable.t%in%colnames(train.cat)) { 
        stop("count.variable.t is defined but is not present in the training catalogue!")
      }
      if (count.variable.r!=''&!count.variable.r%in%colnames(refr.cat)) { 
        stop("count.variable.r is defined but is not present in the reference catalogue!")
      }
      if (zt.label!=''&!zt.label%in%colnames(train.cat)) { 
        stop("zt.label is defined but is not present in the training catalogue!")
      }
      if (zr.label!=''&!zr.label%in%colnames(train.cat)) { 
        stop("zr.label is defined but is not present in the training catalogue!")
      }
      if (zr.label!=''&!zr.label%in%colnames(refr.cat)) { 
        stop("zr.label is defined but is not present in the reference catalogue!")
      }
      if (refr.truth & zt.label!='' & !zt.label%in%colnames(refr.cat)) { 
        stop("zt.label is defined but is not present in the reference catalogue!")
      }
      #}}}
      #Seperate the QC expression into individual terms {{{
      split.expr<-split.expr(qc.expr,ignore=c('abs','na.rm','cbind','rbind','matrixStats::rowMaxs','colMaxs'))
      split.names<-names(split.expr$components)
      keep<-rep(TRUE,length(split.names))
      for (ind in 1:length(split.names)) { 
        keep[ind]<-grepl(split.names[ind],split.expr$replace.expr)
      }
      split.expr$components<-split.expr$components[keep]
      #Check for expression errors {{{
      if (any(grepl("train.cat",split.expr$components)&grepl("refr.cat",split.expr$components))) { 
        stop("QC expression attempts to combine train.cat and refr.cat catalogues (of different lengths!)")
      }
      #}}}
      train.expression<-split.expr$components[which(grepl("train.cat",split.expr$components)&
                                                   !grepl("full.train",split.expr$components))]
      refr.expression<-split.expr$components[which(grepl("refr.cat",split.expr$components)&
                                                  !grepl("full.refr",split.expr$components))]
      full.train.expression<-split.expr$components[which(grepl("full.train",split.expr$components))]
      full.refr.expression<-split.expr$components[which(grepl("full.refr",split.expr$components))]
      #}}}
      #Convert the train & refr cat names to 'data' {{{
      train.expression<-gsub('train.cat','data',train.expression)
      refr.expression<-gsub('refr.cat','data',refr.expression)
      full.train.expression<-gsub('full.train.cat','train.cat',full.train.expression)
      full.refr.expression<-gsub('full.refr.cat','refr.cat',full.refr.expression)
      #}}}
      #If needed, run any full catalogue QC components {{{
      if (any(grepl("full.train",split.expr$components))) {
        #Run the training cat QC components {{{
        full.train.qc.vals<-NULL
        for (expr in full.train.expression) {
          full.train.qc.vals<-cbind(full.train.qc.vals,eval(parse(text=expr)))
        }
        colnames(full.train.qc.vals)<-names(split.expr$components)[which(grepl("full.train",split.expr$components))]
        #}}}
      }
      if (any(grepl("full.refr",split.expr$components))) {
        #Run the reference cat QC components {{{
        full.refr.qc.vals<-NULL
        for (expr in full.refr.expression) {
          full.refr.qc.vals<-cbind(full.refr.qc.vals,eval(parse(text=expr)))
        }
        colnames(full.refr.qc.vals)<-names(split.expr$components)[which(grepl("full.refr",split.expr$components))]
        #}}}
      }
      #}}}
    }
    #}}}
    #Notify {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> Computing per-cell training stats")
    }#}}}
    #Measure the per-cell training properties {{{ 
    if (count.variable.t=='') { 
      stat.expression<-c(count="nrow(data)")
    } else { 
      stat.expression<-c(count="sum(data[[count.variable.t]],na.rm=T)")
    } 
    stat.expression<-c(stat.expression,meanz="mean(data[[zt.label]],na.rm=T)")
    if (do.QC) { 
      stat.expression<-c(stat.expression,train.expression)
    }
    train.cell.props<-generate.kohgroup.property(som=train.som,data=train.cat,
                                              expression=stat.expression,
                                              expr.label=names(stat.expression),
                                              n.cores=som.cores,n.cluster.bins=prod(som.dim),quiet=TRUE)
    if (do.QC) { 
      #Add to the QC frame {{{
      qc.frame<-as.data.table(train.cell.props$property)
      #}}}
    }
    #}}}
    #Notify {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> Computing per-cell reference stats")
    }#}}}
    #Measure the per-cell reference properties {{{ 
    if (count.variable.r=='') { 
      stat.expression<-c(count="nrow(data)",countsq="nrow(data)")
      stat.expression<-c(stat.expression,
                         meanz="mean(data[[zr.label]],na.rm=T)")
      if (refr.truth) { 
        stat.expression<-c(stat.expression,
                           meanz_true="mean(data[[zt.label]],na.rm=T)")
      }
    } else { 
      stat.expression<-c(count="sum(data[[count.variable.r]],na.rm=T)",
                         countsq="sum(data[[count.variable.r]]^2,na.rm=T)")
      stat.expression<-c(stat.expression,
                         meanz="weighted.mean(data[[zr.label]],data[[count.variable.r]],na.rm=T)")
      if (refr.truth) { 
        stat.expression<-c(stat.expression,
                         meanz_true="weighted.mean(data[[zt.label]],data[[count.variable.r]],na.rm=T)")
      }
    } 
    if (do.QC) { 
      stat.expression<-c(stat.expression,refr.expression)
    }
    refr.cell.props<-generate.kohgroup.property(som=refr.som,data=refr.cat,
                                              expression=stat.expression,
                                              expr.label=names(stat.expression),
                                              n.cores=som.cores,n.cluster.bins=prod(som.dim),quiet=TRUE)
    #Add to the QC frame {{{
    qc.frame<-as.data.table(cbind(qc.frame,refr.cell.props$property[,-which(colnames(refr.cell.props$property)=='group.id'),drop=F]))
    #}}}
    #}}}
    #If needed, add full catalogue properties to QC frame {{{
    if (do.QC) { 
      if (any(grepl("full.train",split.expr$components))) {
        #Add to the QC frame {{{
        qc.frame<-as.data.table(cbind(qc.frame,full.train.qc.vals))
        #}}}
      } 
      if (any(grepl("full.refr",split.expr$components))) {
        #Add to the QC frame {{{
        qc.frame<-as.data.table(cbind(qc.frame,full.refr.qc.vals))
        #}}}
      } 
    }
    #}}}
    #Notify {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> Computing Statistics for each Clustering level")
    }#}}}
    #Conglomerate values for each HC {{{ 
    #Training results {{{
    train.res<-foreach(step=1:length(HC.steps))%dopar% {
      #Get the cell-to-clust IDs for this setup
      group.fact<-factor(train.hc.mat[,step],levels=seq(max(HC.steps)))
      #Compute counts per cluster
      count<-sapply(split(train.cell.props$property$count,group.fact),sum)
      #Compute (muz_i*N_i / sum(N_i)) counts per cluster
      meanz<-sapply(split(train.cell.props$property$meanz*train.cell.props$property$count,group.fact),sum,na.rm=T)
      meanz<-meanz/count
      #Return the results
      return=data.frame(count=count,meanz=meanz)
    }
    #}}}
    # Reference results {{{
    refr.res<-foreach(step=1:length(HC.steps))%dopar% {
      group.fact<-factor(train.hc.mat[,step],levels=seq(max(HC.steps)))
      count<-sapply(split(refr.cell.props$property$count,group.fact),sum)
      countsq<-sapply(split(refr.cell.props$property$countsq,group.fact),sum)
      meanz<-sapply(split(refr.cell.props$property$meanz*refr.cell.props$property$count,group.fact),sum,na.rm=T)
      meanz<-meanz/count
      if (refr.truth) { 
        meanz_true<-sapply(split(refr.cell.props$property$meanz_true*refr.cell.props$property$count,group.fact),sum,na.rm=T)
        meanz_true<-meanz_true/count
        return=data.frame(count=count,countsq=countsq,meanz=meanz,meanz_true=meanz_true)
      } else { 
        return=data.frame(count=count,countsq=countsq,meanz=meanz)
      }
    }
    #}}}
    if (do.QC) { 
      #Get the QC results per cluster {{{
      qc.res<-foreach(step=1:length(HC.steps))%dopar% {
        group.fact<-factor(train.hc.mat[,step],levels=seq(max(HC.steps)))
        #Combine the QC components 
        comb.comps<-data.table(group.id=seq(max(HC.steps)))
        #Compute counts per cluster
        for (comp in names(split.expr$components)) { 
          comb.comps[[comp]]<-sapply(split(qc.frame[[comp]],group.fact),mean,na.rm=T)
        }
        #Evaluate the QC
        qc.res.expr<-paste0("data.table(group.id=group.id,QCeval=",split.expr$replace.expr,")")
        qc.result<-comb.comps[,eval(parse(text=qc.res.expr))]
        #Return the QC result 
        return=data.frame(QC=ifelse(qc.result$QCeval,1,0))
      }
      #}}}
      #Reformat the results {{{
      qc.vals<-ctsq.refr<-ct.refr<-ct.train<-muz.refr<-muz.train<-matrix(NA,nrow=length(HC.steps),ncol=max(HC.steps))
      for (step in 1:length(HC.steps)) { 
        #Training Counts per setup 
        ct.train[step,]<-train.res[[step]]$count
        #Training meanz per setup 
        muz.train[step,]<-train.res[[step]]$meanz
        if (refr.truth) { 
          #Reference meanz per setup 
          muz.refr[step,]<-refr.res[[step]]$meanz_true
        }
        #Reference Counts per setup 
        ct.refr[step,]<-refr.res[[step]]$count
        #Reference Counts per setup 
        ctsq.refr[step,]<-refr.res[[step]]$countsq
        #QC per setup
        qc.vals[step,]<-qc.res[[step]]$QC
      }
      #}}}
    } else { 
      #Reformat the results {{{
      ctsq.refr<-ct.refr<-ct.train<-muz.refr<-muz.train<-matrix(NA,nrow=length(HC.steps),ncol=max(HC.steps))
      for (step in 1:length(HC.steps)) { 
        #Training Counts per setup 
        ct.train[step,]<-train.res[[step]]$count
        #Training meanz per setup 
        muz.train[step,]<-train.res[[step]]$meanz
        if (refr.truth) { 
          #Reference meanz per setup 
          muz.refr[step,]<-refr.res[[step]]$meanz_true
        }
        #Reference Counts per setup 
        ct.refr[step,]<-refr.res[[step]]$count
        #Reference Counts squared per setup 
        ctsq.refr[step,]<-refr.res[[step]]$countsq
      }
      #}}}
    }
    #}}}
    #}}}
  } else { 
    #By Datum {{{
    #Sort the training data into the hierarchical clusters {{{
    #Get the individual data-to-SOMcell IDs {{{
    classif<-replicate(length(HC.steps),train.som$unit.classif)
    ij<-as.matrix(expand.grid(seq(nrow(classif)),seq(ncol(classif)))); 
    #}}}
    #Convert the IDs from data-to-SOMcell into data-to-cluster {{{
    train.somclust<-array(train.hc.mat[cbind(classif[ij],col(classif)[ij])],dim=dim(classif))
    #}}}
    #}}}
    #Notify {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> Dispersing reference data into cluster bins")
    }#}}}
    #Sort the refr data into the hierarchical clusters {{{
    #Get the individual data-to-SOMcell IDs {{{
    classif<-replicate(length(HC.steps),refr.som$unit.classif)
    ij<-as.matrix(expand.grid(seq(nrow(classif)),seq(ncol(classif)))); 
    #}}}
    #Convert the IDs from data-to-SOMcell into data-to-cluster {{{
    refr.somclust<-array(train.hc.mat[cbind(classif[ij],col(classif)[ij])],dim=dim(classif))
    #}}}
    #}}}
    #Notify {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> Generating the SOM DIR weights: ")
    }
    #}}}
    #Reference counts {{{
    if (!quiet) { 
      cat("\n       -Reference counts per cluster")
    }
    if (count.variable.r=='') { 
      ct.refr<-matrixStats::colTabulates(refr.somclust,values=1:max(HC.steps))
    } else { 
      if (nrow(refr.somclust)!=nrow(refr.cat)) { 
        stop(paste("somclust and catalogue for reference are of different lengths?!\n",
             nrow(refr.somclust),"!=",nrow(refr.cat),'\n'))
      }
      ct.refr<-colWeightedTabulates(refr.somclust,w=refr.cat[[count.variable.r]],values=1:max(HC.steps),cores=som.cores)
    } 
    if (any(dim(ct.refr)!=c(length(HC.steps),max(HC.steps)))) { 
      stop(paste0("\nReference Counts is not of expected length; likely failure in parallelisation!\n",
           "{",paste(dim(ct.refr),collapse='x'),"} != {",length(HC.steps),"x",max(HC.steps)))
    }
    #}}}
    #Training counts {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n       -Training counts per cluster")
    }
    if (count.variable.t=='') { 
      ct.train<-matrixStats::colTabulates(train.somclust,values=1:max(HC.steps))
    } else {
      if (nrow(train.somclust)!=nrow(train.cat)) { 
        stop(paste("somclust and catalogue for training are of different lengths?!\n",
             nrow(train.somclust),"!=",nrow(train.cat),'\n'))
      }
      ct.train<-colWeightedTabulates(train.somclust,w=train.cat[[count.variable.t]],values=1:max(HC.steps))
    }
    #}}}
    #Training mean zs {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n       -Training <z> per cluster")
    }
    muz.train<-colWeightedTabulates(train.somclust,w=train.cat[[zt.label]],values=1:max(HC.steps))/ct.train
    #}}}
    #}}}
  }
  #}}}
  #Cluster Weights {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n       -Training weights per cluster")
  }
  wt.final<-ct.refr/ct.train
  wt.final[which(!is.finite(wt.final))]<-0
  #}}}
  #Calculate the mean-z for each HC step {{{
  if (do.QC) { 
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n       -Change in mean-z a.f.o. cluster with QC")
    }
    muz<-rowSums(ct.train*qc.vals*wt.final*muz.train,na.rm=TRUE)/rowSums(ct.train*qc.vals*wt.final,na.rm=TRUE)
    mask<-ifelse(wt.final>0,1,0)
    dneff<-(rowSums(ct.refr*qc.vals*mask,na.rm=TRUE)^2/rowSums(ctsq.refr*qc.vals*mask,na.rm=TRUE))/
           (rowSums(ct.refr             ,na.rm=TRUE)^2/rowSums(ctsq.refr             ,na.rm=TRUE))
    if (refr.truth) { 
      muz.true<-rowSums(ct.refr*qc.vals*mask*muz.refr,na.rm=TRUE)/rowSums(ct.refr*qc.vals*mask,na.rm=TRUE)
    }
  } else { 
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n       -Change in mean-z a.f.o. cluster")
    }
    muz<-rowSums(ct.train*wt.final*muz.train,na.rm=TRUE)/rowSums(ct.train*wt.final,na.rm=TRUE)
    mask<-ifelse(wt.final>0,1,0)
    dneff<-(rowSums(ct.refr*mask,na.rm=TRUE)^2/rowSums(ctsq.refr*mask,na.rm=TRUE))/
           (rowSums(ct.refr     ,na.rm=TRUE)^2/rowSums(ctsq.refr     ,na.rm=TRUE))
    if (refr.truth) { 
      muz.true<-rowSums(ct.refr*mask*muz.refr,na.rm=TRUE)/rowSums(ct.refr*mask,na.rm=TRUE)
    }
  }
  #}}}
  #Notify  {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n       - Saving statistics")
  }
  #}}}
  #Save the statistics {{{
  HCoptim<-list(ct.train=ct.train,ct.refr=ct.refr,ctsq.refr=ctsq.refr,muz.train=muz.train,wt.final=wt.final,HC.steps=HC.steps,muz=muz,dneff=dneff)
  if (do.QC) HCoptim$qc.vals<-qc.vals
  if (refr.truth) { HCoptim$muz.true=muz.true }
  save(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0(addstr[loop.num],'_HCoptim.Rdata'),output.file[loop.num],fixed=TRUE)),HCoptim)
  #}}}
  #Notify  {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n       - Determine Optimal Nclust")
  }
  #}}}
  #Select the optimal HC step {{{
  muz.fiducial<-muz[which(HC.steps==factor.nbins)]
  dneff.fiducial<-dneff[which(HC.steps==factor.nbins)]
  HCs.possible<-HC.steps[which(abs(muz-muz.fiducial)<=optimize.z.threshold)]
  HC.optimal<-min(HCs.possible)
  muz.optimal<-muz[which(HC.steps==HC.optimal)]
  if (refr.truth) {  
    muz.true.fiducial<-muz.true[which(HC.steps==factor.nbins)]
    muz.true.optimal<-muz.true[which(HC.steps==HC.optimal)] 
  }
  dneff.optimal<-dneff[which(HC.steps==HC.optimal)]
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat(paste0("\n       -mean-z range: <z> in [",paste(round(digits=3,range(muz,na.rm=T)),collapse=','),"]"))
    cat(paste0("\n       -Fiducial mu_z: <z> = ",round(digits=3,muz.fiducial)," @ factor.nbins = ",factor.nbins))
    if (refr.truth) { 
      cat(paste0("\n       -Fiducial Bias: <z>-<z>_true = ",round(digits=3,muz.fiducial-muz.true.fiducial)," @ factor.nbins = ",factor.nbins))
    }
    cat(paste0("\n       -Fiducial dneff: Delta neff = ",round(digits=3,dneff.fiducial)," @ factor.nbins = ",factor.nbins))
    cat(paste0("\n       -Optimal mu_z:  <z> = ",round(digits=3,muz.optimal)," @ factor.nbins = ",HC.optimal))
    if (refr.truth) { 
      cat(paste0("\n       -Optimal Bias: <z>-<z>_true = ",round(digits=3,muz.optimal-muz.true.optimal)," @ factor.nbins = ",HC.optimal))
    }
    cat(paste0("\n       -Optimal dneff: Delta neff = ",round(digits=3,dneff.optimal)," @ factor.nbins = ",HC.optimal))
    cat(paste0("\n       -Updating training and reference cluster assignments with factor.nbins = ",HC.optimal))
  }
  #}}}
  #Update the factor.nbins, and cluster assignments {{{
  factor.nbins<-HC.optimal
  #Notify /*fold*/ {{{
  if (!quiet) { 
    short.timer<-proc.time()[3]
    cat("\n    -> Re-assigning the training data to clusters")
  }#/*fend*/}}}
  #Get the training positions from the trained SOM /*fold*/ {{{
  train.som<-generate.kohgroups(som=train.som,n.cluster.bins=factor.nbins,
                               train.expr=factor.label,data.missing=data.missing,
                               data.threshold=data.threshold,n.cores=som.cores,max.na.frac=1,quiet=TRUE)
  #/*fend*/}}}
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Re-assigning the reference data to clusters")
  }#/*fend*/}}}
  #Get the data positions from the trained SOM /*fold*/ {{{
  refr.som<-generate.kohgroups(som=refr.som,n.cluster.bins=factor.nbins,
                               train.expr=factor.label,data.missing=data.missing,
                               data.threshold=data.threshold,n.cores=som.cores,max.na.frac=1,quiet=TRUE)
  #/*fend*/}}}
  #}}}
}
#/*fend*/}}}
#Record the cluster assignments and write the SOMs to file {{{
#Training Catalogue {{{
train.cat$GroupFactor<-train.som$clust.classif
save(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0(addstr[loop.num],'_SOMdata.Rdata'),output.file[loop.num],fixed=TRUE)),train.som)
#}}}
#Reference Catalogue {{{
refr.cat$GroupFactor<-refr.som$clust.classif
save(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0(addstr[loop.num],'_refr_SOMdata.Rdata'),output.file[loop.num],fixed=TRUE)),refr.som)
#}}}
#}}}
#Plot the SOM and additional information /*fold*/ {{{
if (plot>0) { 
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Plotting SOM figures")
  }#/*fend*/}}}
  #Plot the som codes /*fold*/ {{{
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_refr_SOM_codes.png',output.file[loop.num],fixed=TRUE)),height=10*res,width=10*res,res=res)
  plot(refr.som, type="codes", bgcol=rainbow(factor.nbins)[train.som$cell.clust],shape='straight',codeRendering='lines',border=NA)
  add.cluster.boundaries(refr.som,train.som$cell.clust,lwd=3,col=hsv(1,0,0,0.8))
  dev.off()
  #/*fend*/}}}
  #Plot the changes at each iteration /*fold*/ {{{
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_refr_SOM_changes.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=5*res,res=res)
  plot(refr.som, type="changes")
  dev.off()
  #/*fend*/}}}
  if (plot>1) { 
    #Plot the various SOM components /*fold*/ {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_refr_SOM_props%02d.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=12*res,res=res)
    #Set the margins /*fold*/ {{{
    par(mar=c(0,0,0,0))
    #/*fend*/}}}
    ##Define the plot layout /*fold*/ {{{
    #suppressWarnings(lay.mat<-matrix(c(1,2,3,4,c(4+1:length(refr.som$codes[[1]][1,]),rep(0,100))),
    #                          ncol=4,nrow=floor(length(refr.som$codes[[1]][1,])/4)+2,byrow=T))
    ##Remove any useless rows 
    #if (any(rowSums(lay.mat)==0)) { 
    #  lay.mat<-lay.mat[-which(rowSums(lay.mat)==0),]
    #}
    #layout(lay.mat)
    ##/*fend*/}}}
    #Plot the counts per SOM cell /*fold*/ {{{
    cell.counts<-plot(refr.som, type="count",shape='straight',border=NA,heatkeywidth=1,zlog=TRUE)
    add.cluster.boundaries(refr.som,train.som$cell.clust,lwd=1,col=hsv(1,0,0,0.3))
    #/*fend*/}}}
    #Plot the neighbour distances per SOM cell /*fold*/ {{{
    plot(refr.som, type="dist.neighbours",shape='straight',border=NA,heatkeywidth=1)
    add.cluster.boundaries(refr.som,train.som$cell.clust,lwd=1,col=hsv(1,0,0,0.3))
    #/*fend*/}}}
    #Plot the quality per SOM cell /*fold*/ {{{
    plot(refr.som, type="quality",shape='straight',border=NA,heatkeywidth=1)
    add.cluster.boundaries(refr.som,train.som$cell.clust,lwd=1,col=hsv(1,0,0,0.3))
    #/*fend*/}}}
    #Plot the number density per SOM cell /*fold*/ {{{
    plot(refr.som, type="count",shape='straight',border=NA,heatkeywidth=1)
    add.cluster.boundaries(refr.som,train.som$cell.clust,lwd=1,col=hsv(1,0,0,0.3))
    #/*fend*/}}}
    for (i in 1:length(refr.som$codes[[1]][1,])) { 
      #Plot the i^th data property per SOM cell /*fold*/ {{{
      plot(refr.som, type = "property", property = refr.som$codes[[1]][,i],
           main=paste0(colnames(refr.som$data[[1]])[i],'(scaled)'), palette.name=BlRd,shape='straight',border=NA,heatkeywidth=1)
      add.cluster.boundaries(refr.som,train.som$cell.clust,lwd=1,col=hsv(1,0,0,0.3))
      #/*fend*/}}}
    }
    #Close the plot device /*fold*/ {{{
    dev.off()
    #/*fend*/}}}
    #/*fend*/}}}
  }
}
#/*fend*/}}}
#Notify /*fold*/ {{{
if (!quiet) { 
  cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
  short.timer<-proc.time()[3]
  cat("\n    -> Generating the SOM DIR weights ")
}
#/*fend*/}}}
#Generate the weights /*fold*/ {{{
if (count.variable.r=='' & count.variable.t=='') { 
  wtind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=TRUE)%dopar%{
    train.index<-which(train.cat$GroupFactor==i)
    refr.index<-which(refr.cat$GroupFactor==i)
    wt<-length(refr.index)/length(train.index)
    return=cbind(i,wt)
  }
} else if (count.variable.r=='') { 
  wtind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=TRUE)%dopar%{
    train.index<-which(train.cat$GroupFactor==i)
    refr.index<-which(refr.cat$GroupFactor==i)
    wt<-length(refr.index)/sum(train.cat[[count.variable.t]][train.index])
    return=cbind(i,wt)
  }
} else if (count.variable.t=='') { 
  wtind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=TRUE)%dopar%{
    train.index<-which(train.cat$GroupFactor==i)
    refr.index<-which(refr.cat$GroupFactor==i)
    wt<-sum(refr.cat[[count.variable.r]][refr.index],na.rm=T)/length(train.index)
    return=cbind(i,wt)
  }
} else { 
  wtind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=TRUE)%dopar%{
    train.index<-which(train.cat$GroupFactor==i)
    refr.index<-which(refr.cat$GroupFactor==i)
    wt<-sum(refr.cat[[count.variable.r]][refr.index],na.rm=T)/sum(train.cat[[count.variable.t]][train.index])
    return=cbind(i,wt)
  }
} 
#Define the SOM weights {{{
train.cat[["SOMweight"]]<-wtind[train.cat$GroupFactor,2]
refr.cat[["SOMweight"]]<-1/wtind[refr.cat$GroupFactor,2]
if (train.flag) { 
  train.cat[["SOMweight"]][which(!is.finite(train.cat[["SOMweight"]]))]<-0
  train.cat[["SOMweight"]]<-ifelse(train.cat[["SOMweight"]]>0,1,0)
}
if (refr.flag) { 
  refr.cat[["SOMweight"]][which(!is.finite(refr.cat[["SOMweight"]]))]<-0
  refr.cat[["SOMweight"]]<-ifelse(refr.cat[["SOMweight"]]>0,1,0)
}
#}}}
#Plot the SOM weights /*fold*/ {{{
if (plot>0) {
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_train_SOMweights.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=5*res,res=res)
  dens<-density(log10(train.cat$SOMweight),kernel='rect',from=-3,to=3,bw=0.1/sqrt(12),na.rm=TRUE)
  magplot(dens,xlab='SOMweight',ylab='PDF',unlog='x')
  dev.off()
}
#/*fend*/}}}
#Plot the SOM weights /*fold*/ {{{
if (plot>0) {
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_refr_SOMweights.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=5*res,res=res)
  dens<-density(log10(refr.cat$SOMweight),from=-5,to=1,bw=0.1/sqrt(12),kernel='rect',na.rm=TRUE)
  magplot(dens,xlab='SOMweight',ylab='PDF',unlog='x')
  dev.off()
}
#Paint the SOM weights /*fold*/ {{{
if (plot>0) {
  wrefr.cell<-wtrain.cell<-rep(NA,som.dim[1]*som.dim[2])
  for (i in 1:factor.nbins) { 
    wtrain.cell[which(train.som$cell.clust==i)]<-train.cat$SOMweight[which(train.cat$GroupFactor==i)[1]]
    wrefr.cell[which(train.som$cell.clust==i)]<-refr.cat$SOMweight[which(refr.cat$GroupFactor==i)[1]]
  }
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_SOMweights_paint.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=7*res,res=res)
  layout(cbind(1,2))
  plot(refr.som, type = "property", property = wrefr.cell,#zlim=c(0,1.5),
       main="SOMweight_refr", palette.name=scale.palette(wrefr.cell,palette=BlRd),shape='straight',border=NA,heatkeywidth=som.dim[1]/20,ncol=1e3,heatkeyborder=NA)
  plot(train.som, type = "property", property = wtrain.cell,#zlim=c(0,1.5),
       main="SOMweight_train", palette.name=scale.palette(wtrain.cell,palette=BlRd),shape='straight',border=NA,heatkeywidth=som.dim[1]/20,ncol=1e3,heatkeyborder=NA)
  dev.off()
}
#/*fend*/}}}
#/*fend*/}}}
#/*fend*/}}} 
#Perform SOM z calibration {{{
if (do.zcalib) { 
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Performing SOM z calibration ")
  }
  #/*fend*/}}}
  #Check the zcalib expression for errors {{{ 
  if (grepl("count.variable.t",zcalib.expr)&count.variable.t=='') { 
    stop("count.variable.t is not defined but is used in the cluster z calibration")
  }
  if (grepl("count.variable.r",zcalib.expr)&count.variable.r=='') { 
    stop("count.variable.r is not defined but is used in the cluster z calibration")
  }
  if (count.variable.t!=''&!count.variable.t%in%colnames(train.cat)) { 
    stop("count.variable.t is defined but is not present in the training catalogue!")
  }
  if (count.variable.r!=''&!count.variable.r%in%colnames(refr.cat)) { 
    stop("count.variable.r is defined but is not present in the reference catalogue!")
  }
  if (zt.label!=''&!zt.label%in%colnames(train.cat)) { 
    stop("zt.label is defined but is not present in the training catalogue!")
  }
  if (zr.label!=''&!zr.label%in%colnames(refr.cat)) { 
    stop("zr.label is defined but is not present in the reference catalogue!")
  }
  #}}}
  #Run the zcalib expression per-group {{{
  #Seperate the zcalib expression into individual terms {{{
  split.expr<-split.expr(zcalib.expr,ignore=c('abs','na.rm','cbind','rbind','matrixStats::rowMaxs','colMaxs'))
  split.names<-names(split.expr$components)
  keep<-rep(TRUE,length(split.names))
  for (ind in 1:length(split.names)) { 
    keep[ind]<-grepl(split.names[ind],split.expr$replace.expr)
  }
  split.expr$components<-split.expr$components[keep]
  #Check for expression errors {{{
  if (any(grepl("train.cat",split.expr$components)&grepl("refr.cat",split.expr$components))) { 
    stop("zcalib expression attempts to combine train.cat and refr.cat catalogues (of different lengths!)")
  }
  #}}}
  train.expression<-split.expr$components[which(grepl("train.cat",split.expr$components)&
                                               !grepl("full.train",split.expr$components))]
  refr.expression<-split.expr$components[which(grepl("refr.cat",split.expr$components)&
                                              !grepl("full.refr",split.expr$components))]
  full.train.expression<-split.expr$components[which(grepl("full.train",split.expr$components))]
  full.refr.expression<-split.expr$components[which(grepl("full.refr",split.expr$components))]
  #}}}
  #Convert the train & refr cat names to 'data' {{{
  train.expression<-gsub('train.cat','data',train.expression)
  refr.expression<-gsub('refr.cat','data',refr.expression)
  full.train.expression<-gsub('full.train.cat','train.cat',full.train.expression)
  full.refr.expression<-gsub('full.refr.cat','refr.cat',full.refr.expression)
  #}}}
  #Compute the training cat and reference cat z calibration {{{
  zcalib.frame<-NULL
  if (any(grepl("train.cat",split.expr$components))) {
    #Run the training cat z calibration components {{{
    train.zcalib.vals<-generate.kohgroup.property(som=train.som,data=train.cat,
                                              expression=train.expression,
                                              expr.label=names(train.expression),
                                              n.cores=som.cores,n.cluster.bins=factor.nbins,quiet=TRUE)
    #}}}
    #Add to the z calibration frame {{{
    zcalib.frame<-as.data.table(train.zcalib.vals$property)
    #}}}
  }
  if (any(grepl("refr.cat",split.expr$components))) {
    #Run the reference cat z calibration components {{{
    refr.zcalib.vals<-generate.kohgroup.property(som=refr.som,data=refr.cat,
                                             expression=refr.expression,
                                             expr.label=names(refr.expression),
                                             n.cores=som.cores,n.cluster.bins=factor.nbins,quiet=TRUE)
    #}}}
    #Add to the z calibration frame {{{
    if (is.null(zcalib.frame)) { 
      zcalib.frame<-as.data.table(refr.zcalib.vals$property)
    } else { 
      zcalib.frame<-as.data.table(cbind(zcalib.frame,refr.zcalib.vals$property[,-which(colnames(refr.zcalib.vals$property)=='group.id'),drop=F]))
    }
    #}}}
  }
  if (any(grepl("full.train",split.expr$components))) {
    #Run the training cat z calibration components {{{
    full.train.zcalib.vals<-NULL
    for (expr in full.train.expression) {
      full.train.zcalib.vals<-cbind(full.train.zcalib.vals,eval(parse(text=expr)))
    }
    colnames(full.train.zcalib.vals)<-names(split.expr$components)[which(grepl("full.train",split.expr$components))]
    #}}}
    #Add to the z calibration frame {{{
    zcalib.frame<-as.data.table(cbind(zcalib.frame,full.train.zcalib.vals))
    #}}}
  }
  if (any(grepl("full.refr",split.expr$components))) {
    #Run the reference cat z calibration components {{{
    full.refr.zcalib.vals<-NULL
    for (expr in full.refr.expression) {
      full.refr.zcalib.vals<-cbind(full.refr.zcalib.vals,eval(parse(text=expr)))
    }
    colnames(full.refr.zcalib.vals)<-names(split.expr$components)[which(grepl("full.refr",split.expr$components))]
    #}}}
    #Add to the z calibration frame {{{
    zcalib.frame<-as.data.table(cbind(zcalib.frame,full.refr.zcalib.vals))
    #}}}
  }
  #}}}
  #Run the zcalib expression {{{
  zcalib.res.expr<-paste0("data.table(group.id=group.id,zcalib=",split.expr$replace.expr,")")
  zcalib.result<-zcalib.frame[,eval(parse(text=zcalib.res.expr))]
  #}}}
  #}}}
  #Check that the zcalib value is valid {{{
  if (all(is.na(zcalib.result$zcalib))) { 
    cat("\n")
    print(summary(zcalib.frame))
    print(zcalib.res.expr)
    print(str(zcalib.result))
    stop("All zcalib results are NA! There are no good associations remaining!")
  }
  #}}}
  #Check that the zcalib is a logical {{{ 
  if (any(!is.numeric(zcalib.result$zcalib))) { 
    #The zcalib expression does not return a numeric value! {{{
    stop("Result of zcalib expression was not interpretable as numeric!")
  } else { 
    train.cat$zt.calib.factor<-NA
    for (i in 1:factor.nbins) { 
      train.cat$zt.calib.factor[which(train.cat$GroupFactor==i)]<-zcalib.result$zcalib[which(zcalib.result$group.id==i)]
    }
    #}}}
    #Apply the calibrations {{{
    train.cat[[paste0(zt.label,"_calib")]]<-train.cat[[zt.label]]+train.cat$zt.calib.factor
    #}}}
  }
  #}}}
}
#}}}
#Perform SOM QC {{{
if (do.QC) { 
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Performing SOM Cluster QC ")
  }
  #/*fend*/}}}
  #Check the qc expression for errors {{{ 
  if (grepl("count.variable.t",qc.expr)&count.variable.t=='') { 
    stop("count.variable.t is not defined but is used in the cluster QC")
  }
  if (grepl("count.variable.r",qc.expr)&count.variable.r=='') { 
    stop("count.variable.r is not defined but is used in the cluster QC")
  }
  if (count.variable.t!=''&!count.variable.t%in%colnames(train.cat)) { 
    stop("count.variable.t is defined but is not present in the training catalogue!")
  }
  if (count.variable.r!=''&!count.variable.r%in%colnames(refr.cat)) { 
    stop("count.variable.r is defined but is not present in the reference catalogue!")
  }
  if (zt.label!=''&!zt.label%in%colnames(train.cat)) { 
    stop("zt.label is defined but is not present in the training catalogue!")
  }
  if (zr.label!=''&!zr.label%in%colnames(refr.cat)) { 
    stop("zr.label is defined but is not present in the reference catalogue!")
  }
  #}}}
  #Run the QC expression per-group {{{
  #Seperate the QC expression into individual terms {{{
  split.expr<-split.expr(qc.expr,ignore=c('abs','na.rm','cbind','rbind','matrixStats::rowMaxs','colMaxs'))
  split.names<-names(split.expr$components)
  keep<-rep(TRUE,length(split.names))
  for (ind in 1:length(split.names)) { 
    keep[ind]<-grepl(split.names[ind],split.expr$replace.expr)
  }
  split.expr$components<-split.expr$components[keep]
  #Check for expression errors {{{
  if (any(grepl("train.cat",split.expr$components)&grepl("refr.cat",split.expr$components))) { 
    stop("QC expression attempts to combine train.cat and refr.cat catalogues (of different lengths!)")
  }
  #}}}
  train.expression<-split.expr$components[which(grepl("train.cat",split.expr$components)&
                                               !grepl("full.train",split.expr$components))]
  refr.expression<-split.expr$components[which(grepl("refr.cat",split.expr$components)&
                                              !grepl("full.refr",split.expr$components))]
  full.train.expression<-split.expr$components[which(grepl("full.train",split.expr$components))]
  full.refr.expression<-split.expr$components[which(grepl("full.refr",split.expr$components))]
  #}}}
  #Convert the train & refr cat names to 'data' {{{
  train.expression<-gsub('train.cat','data',train.expression)
  refr.expression<-gsub('refr.cat','data',refr.expression)
  full.train.expression<-gsub('full.train.cat','train.cat',full.train.expression)
  full.refr.expression<-gsub('full.refr.cat','refr.cat',full.refr.expression)
  #}}}
  #Compute the training cat and reference cat QC {{{
  qc.frame<-NULL
  if (any(grepl("train.cat",split.expr$components))) {
    #Run the training cat QC components {{{
    train.qc.vals<-generate.kohgroup.property(som=train.som,data=train.cat,
                                              expression=train.expression,
                                              expr.label=names(train.expression),
                                              n.cores=som.cores,n.cluster.bins=factor.nbins,quiet=TRUE)
    #}}}
    #Add to the QC frame {{{
    qc.frame<-as.data.table(train.qc.vals$property)
    #}}}
  }
  if (any(grepl("refr.cat",split.expr$components))) {
    #Run the reference cat QC components {{{
    refr.qc.vals<-generate.kohgroup.property(som=refr.som,data=refr.cat,
                                             expression=refr.expression,
                                             expr.label=names(refr.expression),
                                             n.cores=som.cores,n.cluster.bins=factor.nbins,quiet=TRUE)
    #}}}
    #Add to the QC frame {{{
    if (is.null(qc.frame)) { 
      qc.frame<-as.data.table(refr.qc.vals$property)
    } else { 
      qc.frame<-as.data.table(cbind(qc.frame,refr.qc.vals$property[,-which(colnames(refr.qc.vals$property)=='group.id'),drop=F]))
    }
    #}}}
  }
  if (any(grepl("full.train",split.expr$components))) {
    #Run the training cat QC components {{{
    full.train.qc.vals<-NULL
    for (expr in full.train.expression) {
      full.train.qc.vals<-cbind(full.train.qc.vals,eval(parse(text=expr)))
    }
    colnames(full.train.qc.vals)<-names(split.expr$components)[which(grepl("full.train",split.expr$components))]
    #}}}
    #Add to the QC frame {{{
    qc.frame<-as.data.table(cbind(qc.frame,full.train.qc.vals))
    #}}}
  }
  if (any(grepl("full.refr",split.expr$components))) {
    #Run the reference cat QC components {{{
    full.refr.qc.vals<-NULL
    for (expr in full.refr.expression) {
      full.refr.qc.vals<-cbind(full.refr.qc.vals,eval(parse(text=expr)))
    }
    colnames(full.refr.qc.vals)<-names(split.expr$components)[which(grepl("full.refr",split.expr$components))]
    #}}}
    #Add to the QC frame {{{
    qc.frame<-as.data.table(cbind(qc.frame,full.refr.qc.vals))
    #}}}
  }
  #}}}
  #Run the QC expression {{{
  qc.res.expr<-paste0("data.table(group.id=group.id,QCeval=",split.expr$replace.expr,")")
  qc.result<-qc.frame[,eval(parse(text=qc.res.expr))]
  #}}}
  #}}}
  #Check that the QC returned something {{{
  if (length(qc.result$QCeval)==0) { 
    stop("QC result is length 0! There are no good associations remaining!")
  }
  #}}}
  #Check that the QC value is valid {{{
  if (all(is.na(qc.result$QCeval))) { 
    print(str(split.expr))
    print(summary(qc.frame))
    stop("All QC results are NA! There are no good associations remaining!")
  }
  #}}}
  #Check that the QC is a logical {{{ 
  if (any(!is.logical(qc.result$QCeval))) { 
    #Write QC values to catalogue {{{
    warning("Result of QC was not interpretable as logical!\nQC values will be output but not applied!!")
    if (!quiet) { 
      cat(" - WARNING: QC output is not logical! QC will not be applied!!")
      cat("\n    -> Writing SOM QC results to output catalogues ")
    }
    refr.cat$GroupQC<-train.cat$GroupQC<-NA
    for (i in 1:factor.nbins) { 
      train.cat$GroupQC[which(train.cat$GroupFactor==i)]<-qc.result$QCeval[which(qc.result$group.id==i)]
      refr.cat$GroupQC[which(refr.cat$GroupFactor==i)]<-qc.result$QCeval[which(qc.result$group.id==i)]
    }
    #}}}
  } else { 
    #Apply the QC to catalogue {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> Assigning QC results to sources ")
    }
    goodgroups<-qc.result$group.id[which(qc.result$QCeval)]
    meanz.orig<-weighted.mean(train.cat[[zt.label]],train.cat$SOMweight)
    if (count.variable.r!=""){ 
      count<-refr.cat[[count.variable.r]]
      mask<-ifelse(refr.cat$SOMweight>0,1,0)
      neff.orig<-(sum(count*mask)^2/sum((count*mask)^2))/(sum(count)^2/sum(count^2))
    } else { 
      neff.orig<-length(which(refr.cat$SOMweight>0))/nrow(refr.cat)
    } 
    train.cat$SOMweight[which(!train.cat$GroupFactor%in%goodgroups)]<-0
    train.cat$QCFlag<-ifelse(train.cat$GroupFactor%in%goodgroups,0,1)
    refr.cat$SOMweight[which(!refr.cat$GroupFactor%in%goodgroups)]<-0
    refr.cat$QCFlag<-ifelse(refr.cat$GroupFactor%in%goodgroups,0,1)
    meanz.new<-weighted.mean(train.cat[[zt.label]],train.cat$SOMweight)
    if (count.variable.r!=""){ 
      count<-refr.cat[[count.variable.r]]
      mask<-ifelse(refr.cat$SOMweight>0,1,0)
      neff.new<-(sum(count*mask)^2/sum((count*mask)^2))/(sum(count)^2/sum(count^2))
    } else { 
      neff.new<-length(which(refr.cat$SOMweight>0))/nrow(refr.cat)
    } 
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat(paste0("\n       -Non-QC'd mu_z: <z> = ",round(digits=3,meanz.orig)," @ phot repr = ",round(digits=1,neff.orig*100),"%"))
      cat(paste0("\n       -Post-QC mu_z:  <z> = ",round(digits=3,meanz.new)," @ phot repr = ",round(digits=1,neff.new*100),"%"))
      cat("\n    -> Finishing QC ")
    }
    #}}}
  }
  #}}}
}
#}}}
#Generate the reference catalogue plots  /*fold*/ {{{
if (plot>0) { 
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Plotting the Reference data ")
  }
  #/*fend*/}}}
  #Plot the counts zoomed in /*fold*/ {{{
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_counts.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=10*res,res=res)
  layout(cbind(1,2))
  train.count.tab<-table(train.som$unit.classif)
  refr.count.tab<-table(refr.som$unit.classif)
  plot(train.som, type="counts",shape='straight',border=NA,heatkeywidth=som.dim[1]/20,main='Training Sample',palette=scale.palette(log10(train.count.tab)),ncol=1e3,heatkeyborder=NA,zlog=TRUE)
  add.cluster.boundaries(train.som,train.som$cell.clust,lwd=1,col=hsv(1,0,0,0.3))
  plot(refr.som, type="counts",shape='straight',border=NA,heatkeywidth=som.dim[1]/20,main='Reference Sample',palette=scale.palette(log10(refr.count.tab)),ncol=1e3,heatkeyborder=NA,zlog=TRUE)
  add.cluster.boundaries(refr.som,train.som$cell.clust,lwd=1,col=hsv(1,0,0,0.3))
  dev.off()
  #/*fend*/}}}
  #Measure the per-cell zrefr and z_B /*fold*/ {{{
  zrefr.cell<-ztrain.cell<-zrefr.sd.cell<-ztrain.sd.cell<-rep(NA,som.dim[1]*som.dim[2])
  for (i in 1:factor.nbins) { 
    #Assign the data to the cluster if it's SOMcell belongs to the cluster/*fold*/ {{{
    if (length(which(refr.cat$GroupFactor==i))!=0) {
      zrefr.cell[which(train.som$cell.clust==i)]<-median(refr.cat[[zr.label]][which(refr.cat$GroupFactor==i)],na.rm=TRUE) 
      zrefr.sd.cell[which(train.som$cell.clust==i)]<-mad(refr.cat[[zr.label]][which(refr.cat$GroupFactor==i)],na.rm=TRUE) 
    }
    if (length(which(train.cat$GroupFactor==i))!=0) {
      ztrain.cell[which(train.som$cell.clust==i)]<-median(train.cat[[zt.label]][which(train.cat$GroupFactor==i)],na.rm=TRUE) 
      ztrain.sd.cell[which(train.som$cell.clust==i)]<-mad(train.cat[[zt.label]][which(train.cat$GroupFactor==i)],na.rm=TRUE) 
    }
    #/*fend*/}}}
  }
  #/*fend*/}}}
  if (plot>1) { 
    #Plot the inferred z distributions /*fold*/ {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_zz.png',output.file[loop.num],fixed=TRUE)),height=5*res,width=10*res,res=res)
    layout(cbind(1,2))
    par(mar=c(2,2,0,0),oma=c(2,2,1,1))
    #Generate the implied z distributions /*fold*/ {{{
    syn.refr.z<-rep(NA,refr.cat.len)
    syn.train.z<-rep(NA,train.cat.len)
    for (i in 1:factor.nbins) { 
      refr.ind<-which(refr.cat$GroupFactor==i)
      train.ind<-which(train.cat$GroupFactor==i)
      syn.refr.z[refr.ind]<-sample(refr.cat[[zr.label]][refr.ind],size=length(refr.ind),replace=T)
      syn.train.z[train.ind]<-sample(train.cat[[zt.label]][train.ind],size=length(train.ind),replace=T)
    }
    #/*fend*/}}}
    hist2D(train.cat[[zt.label]],syn.train.z,nbins=101,zlog=T,xlim=c(0,2),ylim=c(0,2),xlab='z input',ylab='z sampled',main='Training Sample',asp=0,flip=T)
    hist2D(refr.cat[[zr.label]],syn.refr.z,nbins=101,zlog=T,xlim=c(0,2),ylim=c(0,2),xlab='z input',ylab='z sampled',main='Reference Sample',asp=0,flip=T)
    dev.off()
    #/*fend*/}}}
    #Plot the z distributions zoomed in /*fold*/ {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_zpaint.png',output.file[loop.num],fixed=TRUE)),height=10*res,width=10*res,res=res)
    layout(matrix(1:4,2,2,byrow=T))
    plot(refr.som, type = "property", property = zrefr.cell,#zlim=c(0,1.5),
         main="Z_refr", palette.name=scale.palette(zrefr.cell,palette=BlRd),shape='straight',border=NA,heatkeywidth=som.dim[1]/20,ncol=1e3,heatkeyborder=NA)
    plot(train.som, type = "property", property = ztrain.cell,#zlim=c(0,1.5),
         main="Z_B", palette.name=scale.palette(ztrain.cell,palette=BlRd),shape='straight',border=NA,heatkeywidth=som.dim[1]/20,ncol=1e3,heatkeyborder=NA)
    plot(refr.som, type = "property", property = zrefr.sd.cell,#zlim=c(0,0.5),
         main=paste("mad","Z_refr"), palette.name=scale.palette(zrefr.sd.cell,palette=BlRd),shape='straight',border=NA,heatkeywidth=som.dim[1]/20,ncol=1e3,heatkeyborder=NA)
    plot(train.som, type = "property", property = ztrain.sd.cell,#zlim=c(0,0.5),
         main=paste("mad","Z_B"), palette.name=scale.palette(ztrain.sd.cell,palette=BlRd),shape='straight',border=NA,heatkeywidth=som.dim[1]/20,ncol=1e3,heatkeyborder=NA)
    dev.off()
    #/*fend*/}}}
  }
  #Plot the effect of the SOM weights /*fold*/ {{{
  pdf(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),'_wgt_colours.pdf',output.file[loop.num],fixed=TRUE)),height=5,width=10)
  par(mar=c(2,2,2,2),oma=c(3,3,3,1))
  suppressWarnings(lay.mat<-matrix(c(1:(length(factor.label)+1),rep(0,100)),
                            ncol=4,nrow=min(c(floor(length(factor.label)/4)+2,3)),byrow=T))
  #Remove any useless rows /*fold*/ {{{
  if (any(rowSums(lay.mat)==0)) { 
    lay.mat<-lay.mat[-which(rowSums(lay.mat)==0),]
  }#/*fend*/}}}

  layout(lay.mat)
  for (factor.expr in factor.label) { 
    #Calculate the expression result
    train.value<-train.cat[,eval(parse(text=factor.expr))]
    refr.value<-refr.cat[,eval(parse(text=factor.expr))]
    #seperate out the components
    seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",factor.expr),' ')))
    refr.good<-rep(TRUE,length=refr.cat.len)
    train.good<-rep(TRUE,length=train.cat.len)
    for (label in seperated.labels) { 
      refr.good<- refr.good & (abs(refr.cat[[label]])<90)
      train.good<- train.good & (abs(train.cat[[label]])<90)
    }
    lims<-quantile(refr.value[refr.good],probs=c(0.05,0.95))
    train.dens<-density(train.value[train.good],bw=0.20/sqrt(12),kern='rect',
                  from=floor(min(lims)),to=ceiling(max(lims)))
    refr.dens<-density(refr.value[refr.good],bw=0.20/sqrt(12),
                  kern='rect',from=floor(min(lims)),to=ceiling(max(lims)))
    magplot(refr.dens,col='darkred',lwd=2,lty=1,ylim=c(0,(max(c(train.dens$y,refr.dens$y)))),
            side=1:4,label=c(T,T,F,F),xlab=gsub("MAG_GAAP_","",factor.expr),ylab='PDF')
    lines(train.dens,col='darkblue',lwd=2,lty=1)
  }
  plot(1,1,type='n',axes=F,xlab='',ylab='')
  legend('left',legend=c("Raw Reference","Raw Training"),lwd=2,lty=c(1,1),col=c('darkred','darkblue'),bty='n',cex=1.2)

  layout(lay.mat)
  for (factor.expr in factor.label) { 
    #Calculate the expression result
    train.value<-train.cat[,eval(parse(text=factor.expr))]
    refr.value<-refr.cat[,eval(parse(text=factor.expr))]
    #seperate out the components
    seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",factor.expr),' ')))
    refr.good<-rep(TRUE,length=refr.cat.len)
    train.good<-rep(TRUE,length=train.cat.len)
    for (label in seperated.labels) { 
      refr.good<- refr.good & (abs(refr.cat[[label]])<90)
      train.good<- train.good & (abs(train.cat[[label]])<90)
    }
    lims<-quantile(refr.value[refr.good],probs=c(0.05,0.95))
    suppressWarnings(
    train.wgt.dens<-density(train.value[train.good],bw=0.20/sqrt(12),
            weight=train.cat$SOMweight[train.good]/sum(train.cat$SOMweight,na.rm=T),
            kern='rect',from=floor(min(lims)),to=ceiling(max(lims)))
    )
    refr.dens<-density(refr.value[refr.good],bw=0.20/sqrt(12),
                  kern='rect',from=floor(min(lims)),to=ceiling(max(lims)))
    magplot(refr.dens,col='darkred',lwd=2,lty=1,ylim=c(0,(max(c(train.wgt.dens$y,refr.dens$y)))),
            side=1:4,label=c(T,T,F,F),xlab=gsub("MAG_GAAP_","",factor.expr),ylab='PDF')
    lines(train.wgt.dens,col='darkblue',lwd=2,lty=1)
  }
  plot(1,1,type='n',axes=F,xlab='',ylab='')
  legend('left',legend=c("Raw Reference","Weighted Training"),lwd=2,lty=c(1,1),col=c('darkred','darkblue'),bty='n',cex=1.2)

  if (plot>1) { 
    layout(lay.mat)
    for (factor.expr in factor.label) { 
      #Calculate the expression result
      train.value<-train.cat[,eval(parse(text=factor.expr))]
      refr.value<-refr.cat[,eval(parse(text=factor.expr))]
      #seperate out the components
      seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",factor.expr),' ')))
      refr.good<-rep(TRUE,length=refr.cat.len)
      train.good<-rep(TRUE,length=train.cat.len)
      for (label in seperated.labels) { 
        refr.good<- refr.good & (abs(refr.cat[[label]])<90)
        train.good<- train.good & (abs(train.cat[[label]])<90)
      }
      lims<-quantile(refr.value[refr.good],probs=c(0.05,0.95))
      suppressWarnings(
      refr.wgt.dens<-density(refr.value[refr.good],bw=0.20/sqrt(12),
              weight=refr.cat$SOMweight[refr.good]/sum(refr.cat$SOMweight,na.rm=T),
              kern='rect',from=floor(min(lims)),to=ceiling(max(lims)))
      )
      train.dens<-density(train.value[train.good],bw=0.20/sqrt(12),kern='rect',
                    from=floor(min(lims)),to=ceiling(max(lims)))
      magplot(refr.wgt.dens,col='darkred',lwd=2,lty=1,ylim=c(0,(max(c(train.dens$y,refr.wgt.dens$y)))),
              side=1:4,label=c(T,T,F,F),xlab=gsub("MAG_GAAP_","",factor.expr),ylab='PDF')
      lines(train.dens,col='darkblue',lwd=2,lty=1)
    }
    plot(1,1,type='n',axes=F,xlab='',ylab='')
    legend('left',legend=c("Weighted Reference","Raw Training"),lwd=2,lty=c(1,1),col=c('darkred','darkblue'),bty='n',cex=1.2)
    dev.off()
    #/*fend*/}}}
  }
}
#/*fend*/}}}
#/*fend*/}}}

#Do we want to write reduced size catalogues? /*fold*/ {{{
if (short.write) { 
  seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",factor.label),' ')))
  if (any(sapply(train.cat[,seperated.labels,with=F],class)=='function')) { 
    seperated.labels[-which(sapply(train.cat[1,seperated.labels,with=F],class)=='')]<-0
  }
  shortcol<-c(seperated.labels,count.variable.t,count.variable.r,zt.label,
              zr.label,"GroupFactor","SOMweight")
  if (do.zcalib) { 
    shortcol<-c(shortcol,"zt.calib.factor",paste0(zt.label,"_calib"))
  } 
  if (do.QC) { 
    shortcol<-c(shortcol,"GroupQC","QCFlag")
  } 
  train.cat<-train.cat[,which(colnames(train.cat)%in%shortcol),with=F]
  refr.cat<-refr.cat[,which(colnames(refr.cat)%in%shortcol),with=F]
}
#/*fend*/}}}
#Output the DIR weighted catalogues /*fold*/ {{{
if(!quiet) { 
  cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
  cat(paste("\n  # DIR_weights Training catalogue name:",sub(paste0('.',output.ending[loop.num]),paste0('_DIRsom',addstr[loop.num],'.',output.ending[loop.num]),output.file[loop.num],fixed=TRUE),'\n')) 
  cat(paste("  > Outputing DIR weights catalogue")) 
  timer<-proc.time()[3]
}
#Output the training catalogue /*fold*/ {{{
print(paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_DIRsom',addstr[loop.num],'.',output.ending[loop.num]),output.file[loop.num],fixed=TRUE)))
write.file(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_DIRsom',addstr[loop.num],'.',output.ending[loop.num]),output.file[loop.num],fixed=TRUE)),train.cat,quote=F,row.names=F,extname='OBJECTS')
#/*fend*/}}}
#Notify /*fold*/ {{{
if(!quiet) { 
  io.clock<-io.clock+proc.time()[3]-timer
  cat(paste(" - Done",as.time(proc.time()[3]-timer,digits=0),"\n")) 
}
#/*fend*/}}}
if(!quiet) { 
  cat(paste("\n  # DIR_weights Reference catalogue name:",sub(paste0('.',output.ending[loop.num]),paste0('_refr_DIRsom',addstr[loop.num],'.',output.ending[loop.num]),output.file[loop.num],fixed=TRUE),'\n')) 
  cat(paste("  > Outputing DIR weights catalogue")) 
  timer<-proc.time()[3]
}
#Output the reference catalogue /*fold*/ {{{
write.file(file=paste0(output.path,'/',sub(paste0('.',output.ending[loop.num]),paste0('_refr_DIRsom',addstr[loop.num],'.',output.ending[loop.num]),output.file[loop.num],fixed=TRUE)),refr.cat,quote=F,row.names=F,extname='OBJECTS')
#/*fend*/}}}
#Notify /*fold*/ {{{
if(!quiet) { 
  io.clock<-io.clock+proc.time()[3]-timer
  cat(paste(" - Done",as.time(proc.time()[3]-timer,digits=0),"\n")) 
}
#/*fend*/}}}
#/*fend*/}}}

#Notify & Loop /*fold*/ {{{
if(!quiet & loop.num!=loop.length) { 
  cat(paste("Loop completed. Elapsed time:",as.time(proc.time()[3]-start.timer,digits=0),
            "; Approximate time remaining:",
            as.time((proc.time()[3]-start.timer)/(loop.num-loop.start+1)*(loop.length-loop.num),digits=0),"\n")) 
}
#/*fend*/}}}
#/*fend*/}}}
}

#Notify and End /*fold*/ {{{
if(!quiet) { cat(paste("DIR weight Generation completed in",as.time(proc.time()[3]-start.timer,digits=0),"(",as.time(io.clock,digits=0),"catalogue IO)\n\n")) }
#/*fend*/}}}

