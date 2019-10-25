#!/usr/bin/Rscript
#
# Script to generate DIR weights from an input refr-z catalogue
# Created by: A.H. Wright (16-01-2019)
#

#Function prints the help document {{{
.help.print<-function() { 
  cat(paste("\nDIR_som.R: script for creating DIR weights for an arbitrary survey\n\n",
            "Script calling Syntax:\n",
            "  Rscript DIR_som.R [options] -i <input catalogues> -k <catalogue keywords>\n",
            "For full example running syntax, run with the '-hd' switch.\n",
            "Required command-line parameters:\n",
            "                -i : (list)   input catalogues. Should be first the reference catalogue, then training catalogue(s).\n",
            "                -k : (list)   catalogue keywords. Should be first the redshift variable names for the reference & training catalogues, and then the catalogue keywords used in SOM training.\n",
            "                              Importantly, the training keywords can be expressions (i.e. combinations of magnitudes to create colours). \n",
            "Example: \n",
            "  Rscript DIR_som.R -i MyPhotometricGalaxyCat.csv MySpectroscopicTrainingCat.csv -k z_phot z_spec u_mag-g_mag g_mag-i_mag r_mag Z_mag\n",
            "Available options:\n",
            "\nTraining data manipulation:\n",
            "                -r : [switch] rescale (i.e. whiten) the input data prior to training the SOM\n",
            "               -nr : [switch] do not rescale (i.e. whiten) the input data prior to training the SOM\n",
            "\nSOM construction and training:\n",
            "                -t : [switch] use a toroidal (i.e. non-flat) SOM topology\n",
            "               -nt : [switch] use a flat (i.e. non-toroidal) SOM topology\n",
            "            --topo : (string) pixel topology for the SOM; must be either 'rectangular' or 'hexagonal'\n",
            "      --sparse.som : (float)  fraction of the input training catalogue to use (sparse sampled) when training the SOM\n",
            "      --sparse.var : (string) variable to sample over when constructing the `sparse' training sample for the SOM \n",
            "     --som.dim -sd : (2x int) the dimension of the SOM \n",
            "    --som.rate -sr : (2x float) the learning rate of the SOM (linear decrease from n to m)\n",
            "  --som.method -sm : (string) the method for SOM iteration \n",
            "   --som.cores -sc : (int)    the number of cores used during SOM computation \n",
            "    --som.iter -si : (int)    the number of SOM iterations \n",
            "         --old.som : (path)   path to a pre-trained SOM for the pipeline to use \n",
            "        --only.som : [switch] only output the SOM, then quit. \n",
            "\nCalibration Weight Generation:\n",
            "--factor.nbins -fn : (int)    the number of hierarchical clusters generated from the FACTORs \n",
            "               -ct : (string) the counting/weighting variable for the training catalogue(s) \n",
            "               -cr : (string) the counting/weighting variable for the reference catalogue\n",
            "\nFile handling and Scripting:\n",
            "                -o : (path)   the path to the output folder for the script \n",
            "               -of : (string) the base output filename for the script\n",
            "               -as : (string) string to append the the output filenames \n",
            "               -ls : (int) the starting loop number to count from when running multiple catalogues (useful for scripting) \n",
            "           --force : [switch] force the construction of output catalogues, even if they already exist \n",
            "\nLDAC tools (if using LDAC filetypes):\n",
            "                -l : (path)   the path to the LDAC binary file ldactoasc. Default: `which ldactoasc`\n",
            "              -lor : (list)   the list of LDAC options to use when reading the reference catalogue\n",
            "              -lot : (list)   the list of LDAC options to use when reading the training catalogue(s)\n",
            "\nPlotting & Runtime:\n",
            "                -p : [switch] create plots of the data and SOMs\n",
            "               -pp : [switch] create even more plots of the data and SOMs\n",
            "               -np : [switch] do not create plots of the data and SOMs\n",
            "            --seed : (int)    seed number to use for Randoms generation. Default: 666\n",
            "            --test : [switch] run the pipeline in testing mode; uses a low-res defaults and a thinned input catalogue\n",
            "                -q : [switch] execute quietly (no prompts)\n",
            "                -h : [switch] print this help\n",
            "               -hd : [switch] print the default parameter values\n\n\n"))
  quit()
}
.default.print<-function() { 
  cat(paste("\ngenerate_randoms_2D.R: script for creating 2D randoms for an arbitrary survey geometry\n",
            "Script calling Syntax:\n",
            "  Rscript DIR_som.R [options] -i <InputReferenceCat> <InputTrainingCat1> <InputTrainingCat2>...\n",
            "Running without specifying any additional options (see the help by running with '-h') causes\n",
            "internal default parameters to be used for a number of settings. This means that running:\n",
            "  Rscript DIR_som.R -i ... -k ... \n",
            "is equivalent to running: \n",
            "  Rscript DIR_som.R -i ... -k ... --seed 666 -r -t --topo hexagonal -sd 55 55 -sr 0.05 0.01 -sm pbatch -sc -1 \\\n",
            "                    -si 100 -fn 100 \n",
            "These internal defaults are:\n",
            "            --seed : (666)       seed number to use for Randoms generation. Default: 666\n",
            "                -r : [switch]    rescale (i.e. whiten) the input data prior to training the SOM\n",
            "                -t : [switch]    use a toroidal (i.e. non-flat) SOM topology\n",
            "            --topo : (hexagonal) pixel topology for the SOM; must be either 'rectangular' or 'hexagonal'\n",
            "     --som.dim -sd : (55 55)     the dimension of the SOM \n",
            "    --som.rate -sr : (0.05 0.01) the learning rate of the SOM (linear decrease from n to m)\n",
            "  --som.method -sm : (pbatch)    the method for SOM iteration \n",
            "   --som.cores -sc : (-1)        the number of cores used during SOM computation (-1 == All available)\n",
            "    --som.iter -si : (100)       the number of SOM iterations \n",
            "--factor.nbins -fn : (100)       the number of hierarchical clusters generated from the FACTORs\n",
            "Running with the --test switch causes these settings to change:\n",
            "                       --> -sd : 12 12 \n",
            "                       --> -si : 20 \n",
            "                       --> reduces input training and reference catalogues to 20% original size\n",
            "\n\nNB: all parameters are defined in the order of specification. I.e. running with:\n",
            "Rscript DIR_som.R -p -pp -np -i ... -k ...\n",
            "will cause no plots to be output (because we set: some plots (-p), then all plots (-pp), and then no plots (-np)).\n",
            "The _exception_ to this rule is '--test': the test parameters are set to their lower values BEFORE reading all other parameters.\n",
            "This is so one can always override the lower 'test' settings, if desired. \n\n"))
#}}}
  quit()
}
#}}}

#Read the command line options {{{
inputs<-commandArgs(TRUE)
if (length(inputs)==0) { .help.print() } 
#}}}

#Start the timer and load the required Packages {{{
start.timer<-proc.time()[3]
suppressWarnings(suppressPackageStartupMessages(require(kohonen)))
suppressWarnings(suppressPackageStartupMessages(require(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(require(helpRfuncs)))
suppressWarnings(suppressPackageStartupMessages(require(LAMBDAR)))
suppressWarnings(suppressPackageStartupMessages(require(KernSmooth)))
suppressWarnings(suppressPackageStartupMessages(require(itertools)))
#}}}

#Functions {{{
#Function to scale a colour palette to desired quantiles {{{
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
#}}}

#Create the function to plot the Density Bar {{{
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
#}}}

#Function sets a new subplot within the current plot area {{{
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
#}}}

#Function sets the mainplot panel for a new figure assuming a {{{
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
#}}}

#Function plots the segments for each hierarchical cluster {{{
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
#}}}

#Print Rounding functions {{{
fsignif<-function(x,digits) gsub('-','\U2212',sapply(signif(x,digits), sprintf, fmt=paste0("%#.",digits,"g")))
fround<-function(x,digits) gsub('-','\U2212',gsub('-0.00','0.00',(sapply(round(x,digits), sprintf, fmt=paste0("%#.",digits,"f")))))
#}}}
#Show bandwidth function {{{
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
#}}}

#The colour palette {{{
BlRd<-function(n,alpha=1){rainbow(n,end=2/3,alpha=alpha)[n:1]}
#}}}
#}}}

#Define the colour palette and figure placement{{{
frac.fig<-c(0.1,0.9,0.1,0.2)
BlBuRd<-colorRampPalette(c('black',rev(brewer.pal(10,"RdBu"))))
WtBuRd<-colorRampPalette(c('white',rev(brewer.pal(10,"RdBu")[-(5)])))
#}}}

#Read the options {{{
#Default parameter values {{{
only.som<-force<-sparse.som<-reuse<-useMult<-quiet<-FALSE
input.cat.opt<-input.key.opt<-FALSE
loop.start<-plot<-0
sparse.min.num<-1
seed<-666
res<-200
maxNAfrac=1
detect.limit<-80
missing.val<--99
count.variable.1<-count.variable.2<-''
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
catalogues<-NULL
rescale<-TRUE
testing<-FALSE
#}}}
#Loop through the command arguments {{{
while (length(inputs)!=0) {
  #Check the options syntax {{{
  while (length(inputs)!=0 && inputs[1]=='') { inputs<-inputs[-1] }  
  if (!grepl('^-',inputs[1])) {
    print(inputs)
    stop(paste("Incorrect options provided!",
               "Check the lengths for each option!\n",
               "Only -i and -k parameters can have more than 1 item"))
  }
  #Check for test variable {{{
  if (any(inputs=='--test')) { 
    inputs<-c(inputs[which(inputs=='--test')],inputs[which(inputs!='--test')])
  }
  #}}}
  #}}}
  if (inputs[1]=='-p') { 
    #Create the plots {{{
    plot<-1
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-pp') {
    #Do not create the plots {{{
    plot<-2
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-np') {
    #Do not create the plots {{{
    plot<-0
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-r') {
    #rescale input parameters {{{
    rescale<-TRUE
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-nr') {
    #do not rescale input parameters {{{
    rescale<-FALSE
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-as') {
    #read the addstring  {{{
    inputs<-inputs[-1]
    addstr<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-i') {
    #Read the input catalogues {{{
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) { 
      refr.catpath<-inputs[1:(which(grepl('^-',inputs))[1]-1)][1]
      catalogues<-inputs[1:(which(grepl('^-',inputs))[1]-1)][-1]
      inputs<-inputs[-(1:(which(grepl('^-',inputs))[1]-1))]
    } else { 
      refr.catpath<-inputs[1]
      catalogues<-inputs[-1]
      inputs<-NULL
    } 
    input.cat.opt<-TRUE
    #}}}
  } else if (inputs[1]=='-cr') {
    #Define the count variable for the reference cat {{{
    inputs<-inputs[-1]
    count.variable.1<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-ct') {
    #Define the count variable for the training cat(s) {{{
    inputs<-inputs[-1]
    count.variable.2<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-k') {
    #Read the catalogue keywords {{{
    if (any(grepl('^-',inputs[-1]))) {
      key.ids<-1:(which(grepl('^-',inputs[-1]))[1])
    } else { 
      key.ids<-1:length(inputs)
    }
    zr.label<-inputs[key.ids[2]]
    zt.label<-inputs[key.ids[3]]
    keys<-inputs[key.ids[c(-1,-2,-3)]]
    inputs<-inputs[-(key.ids)]
    if (length(key.ids)<4) { 
      stop("There are no columns to group the data on!")
    } else if (length(key.ids)>=4) { 
      suppressWarnings(suppressPackageStartupMessages(require(kohonen)))
      #factor.label<-keys[-(1:3)]
      factor.label<-keys
    } else { 
      stop(paste0("Catalogue Keywords must have length 3 or greater. It is ",length(key.ids),":-> ",keys," <-\n",
                  "i.e.:\n",
                  "Rscript DIR_som.R -k z_ref Z_B MAG_GAAP_r -i <InputReferenceCat> <InputTrainingCat1> ...\n",
                  "Rscript DIR_som.R -k z_ref Z_B MAG_GAAP_u MAG_GAAP_g MAG_GAAP_r MAG_GAAP_i -i <InputReferenceCat> <InputTrainingCat1> ...\n",
                  "Rscript DIR_som.R -k z_ref Z_B MAG_GAAP_u MAG_GAAP_g MAG_GAAP_r MAG_GAAP_i MAG_GAAP_Z MAG_GAAP_Y MAG_GAAP_J MAG_GAAP_H MAG_GAAP_Ks -i <InputReferenceCat> <InputTrainingCat1> ...\n"))
    }
    #Create the full keys string
    #keys<-paste(keys,collapse=' ')
    keys<-paste('-k',paste(unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",keys),' '))),collapse=' '))
    #}}}
    input.key.opt<-TRUE
  } else if (inputs[1]=='-of') {
    #Define the output directory  {{{
    inputs<-inputs[-1]
    output.file<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-o') {
    #Define the output directory  {{{
    inputs<-inputs[-1]
    output.path<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-ls') {
    #Define the loop.start counter {{{
    inputs<-inputs[-1]
    loop.start<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-l') {
    #Define the ldactoasc binary  {{{
    inputs<-inputs[-1]
    ldactoasc<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-lor') {
    #Define the ldacoptions for reference catalogue {{{
    inputs<-inputs[-1]
    ldac.options.1<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-lot') {
    #Define the ldacoptions for training catalogue(s)  {{{
    inputs<-inputs[-1]
    ldac.options.2<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-nt') {
    #Define whether the SOM is toroidal {{{
    inputs<-inputs[-1]
    som.toroidal<-FALSE
    #}}}
  } else if (inputs[1]=='--topo') {
    #Define whether the SOM is toroidal {{{
    inputs<-inputs[-1]
    som.topo<-inputs[1]
    inputs<-inputs[-1]
    if (!som.topo%in%c("rectangular","hexagonal")) { 
      stop(paste0('SOM Topology must be "rectangular" or "hexagonal", not: ',som.topo))
    }
    #}}}
  } else if (inputs[1]=='-t') {
    #Define whether the SOM is toroidal {{{
    inputs<-inputs[-1]
    som.toroidal<-TRUE
    #}}}
  } else if (inputs[1]=='--factor.nbins'|inputs[1]=='-fn') {
    #Define the number of factor bins {{{
    inputs<-inputs[-1]
    factor.nbins<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='--sparse.var') {
    inputs<-inputs[-1]
    sparse.var<-inputs[1]
    inputs<-inputs[-1]
  } else if (inputs[1]=='--sparse.som') {
    #Load an already calculated som from file {{{
    inputs<-inputs[-1]
    sparse.frac<-as.numeric(inputs[1])
    if (sparse.frac>=1) { 
      stop("Sparse sampling fraction must be less than 1!") 
    } else if (sparse.frac <=0) { 
      stop("Sparse sampling fraction must be greater than 0!") 
    }
    inputs<-inputs[-1]
    sparse.som<-TRUE
    #}}}
  } else if (inputs[1]=='--old.som') {
    #Load an already calculated som from file {{{
    inputs<-inputs[-1]
    som.data.file<-inputs[1]
    inputs<-inputs[-1]
    reuse<-TRUE
    #}}}
  } else if (inputs[1]=='--som.dim'|inputs[1]=='-sd') {
    #Define the SOM dimension {{{
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) {
      dim.ids<-1:(which(grepl('^-',inputs))[1]-1)
    } else { 
      dim.ids<-1:length(inputs)
    }
    if (length(dim.ids)!=2) { stop("som.dim must be of length 2; i.e. --som.dim 55 55") } 
    som.dim<-as.numeric(inputs[dim.ids])
    inputs<-inputs[-dim.ids]
    #}}}
  } else if (inputs[1]=='--som.rate'|inputs[1]=='-sr') {
    #Define the SOM rate {{{
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) {
      rate.ids<-1:(which(grepl('^-',inputs))[1]-1)
    } else { 
      rate.ids<-1:length(inputs)
    }
    if (length(rate.ids)!=2) { stop("som.rate must be of length 2; i.e. --som.rate 0.05 0.01") } 
    som.rate<-as.numeric(inputs[rate.ids])
    inputs<-inputs[-rate.ids]
    #}}}
  } else if (inputs[1]=='--som.method'|inputs[1]=='-sm') {
    #Define the SOM method {{{
    inputs<-inputs[-1]
    som.method<-inputs[1]
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='--som.cores'|inputs[1]=='-sc') {
    #Define the number of SOM cores {{{
    inputs<-inputs[-1]
    som.cores<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='--som.iter'|inputs[1]=='-si') {
    #Define the number of SOM iterations {{{
    inputs<-inputs[-1]
    som.iter<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='--seed') {
    #Define the seed for SOM generation {{{
    inputs<-inputs[-1]
    seed<-as.numeric(inputs[1])
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='--force') { 
    #Force catalogue creation {{{
    force<-TRUE
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='--only.som') { 
    #Only output the SOM {{{
    only.som<-TRUE
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='--test') { 
    #Run in testing mode {{{
    testing<-TRUE
    inputs<-inputs[-1]
    #Set the testing defaults {{{
    som.iter<-20
    som.dim<-c(12,12)
    #}}}
    #}}}
  } else if (inputs[1]=='-q') {
    #Run Quietly {{{
    quiet<-TRUE
    inputs<-inputs[-1]
    #}}}
  } else if (inputs[1]=='-hd') {
    #Print the default values help function{{{
    .default.print()
    return(NULL)
    #}}}
  } else if (inputs[1]=='-h') {
    #Print the help function{{{
    .help.print()
    return(NULL)
    #}}}
  } else {
    stop(paste("Unknown option",inputs[1]))
  }
}
#}}}
#}}}

#Check for the required options {{{
if (!(input.cat.opt & input.key.opt)) { 
  stop("You must specify the input catalogues and the input keywords. Run with -h and -hd for help.")
}
#}}}

#Check for the output path {{{
if (!dir.exists(output.path)) { 
  dir.create(output.path,recursive=TRUE)
}
#}}}

#Check for the input catalogues {{{
if (length(catalogues)==0) { 
  stop("No catalogues provided! Call Syntax:\nRscript DIR_som.R [options] -i <InputReferenceCat> <InputTrainingCat1> ...") 
}
#}}}

#Prompt and start {{{
if (!quiet) { 
  cat("Creating DIR weights with input training catalogue(s):\n ") 
  if (length(catalogues) > 5) { 
    cat(paste0('  --> ',catalogues[1:3],'\n'))
    cat(paste0('  ....\n'))
    cat(paste0('  --> ',rev(catalogues)[3:1],'\n'))
  } else {
    cat(paste0('  --> ',catalogues,'\n'))
  }
}
#}}}

#Set the seed for randomisation {{{
set.seed(seed)
#}}}

#Initialise loop counters {{{
io.clock<-0
loop.num<-loop.start-1 #this is just here for name convenience when the pipe errors
loop.length<-length(catalogues)+loop.start-1
#}}}

#Read in the input catalogue {{{
if (!quiet) { 
  cat(paste("  > Reading Reference Catalogue")) 
  timer<-proc.time()[3]
}
if (grepl('.cat',refr.catpath,fixed=T)) { 
  if (!quiet) { cat(" (with ldactoasc):\n\n") }
  refr.cat<-try(read.ldac(refr.catpath,ldactoasc=ldactoasc,data.table=TRUE,
                     diagnostic=TRUE ,options=paste('-s',keys,zr.label,count.variable.1,ldac.options.1),force=TRUE,clean=TRUE,showProgress=FALSE),silent=TRUE)
  if (suppressWarnings(class(refr.cat)=='try-error')) { 
    cat(" - ")
    stop(paste0("Failure to read input catalogue. Check your ldactoasc binary?\n",
                "  NB: We assume that input \'.cat\' files are in LDAC format. Table details can be specified at input.\n",
                "  For other formats we use: .fits (FITS standard); .csv/.asc (CSV,TSV,ASCII); .Rdata (R data file)\n"))
  }
  if (!quiet) { cat("\n  -----> Completed!\n") }
} else if (grepl('.Rdata',refr.catpath,fixed=T)) {
  nam<-try(load(refr.catpath),silent=TRUE)
  if (suppressWarnings(class(nam)=='try-error')) { 
    cat(" - ")
    stop(paste0("Failure to load input Rdata catalogue.\n"))
  }
  if (nam=="nam") {
    #Whoops, we overwrote the catalogue while checking its name! {{{
    load(refr.catpath)
    refr.cat<-nam
    rm('nam')
    #}}}
  } else if (nam!="refr.cat") { 
    refr.cat<-get(nam)
    rm(nam)
  } 
} else if (grepl('.fits',refr.catpath,fixed=T)) {
  refr.cat<-try(read.fits.cat(refr.catpath,data.table=TRUE),silent=TRUE)
  if (suppressWarnings(class(refr.cat)=='try-error')) { 
    cat(" - ")
    stop(paste0("Failure to read input fits catalogue.\n"))
  }
} else if (grepl('.csv',refr.catpath,fixed=T) | grepl('.asc',refr.catpath,fixed=T)) {
  refr.cat<-try(fread(refr.catpath,data.table=TRUE,showProgress=FALSE),silent=TRUE)
  if (suppressWarnings(class(refr.cat)=='try-error')) { 
    cat(" - ")
    stop(paste0("Failure to read input CSV catalogue.\n"))
  }
} else { 
  cat(" - ")
  stop(paste0("Failure to read input catalogue. Unknown file type. Expects .cat for ldac, .fits for binary, or .csv\n"))
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

#}}}

for (train.catpath in catalogues) { 
#Loop through catalogues {{{
loop.num<-loop.num+1
#Get the catalogue name and ending from the catalogue path given #{{{
if (!quiet) { cat(paste("Working on Catalogue:",train.catpath,"\n")) }
train.catnam<-rev(strsplit(train.catpath,'/')[[1]])[1]
train.ending<-rev(strsplit(train.catnam,'.',fixed=TRUE)[[1]])[1]
#Check for the output catalogues {{{
if (!exists('output.file')) { 
  output.file<-train.catnam
  #Define the filename addition {{{
  addstr<-paste0('_SOM_',paste(factor.label,collapse='_'))
  
  if (nchar(addstr>25)) { 
    print.length<-length(which(cumsum(nchar(gsub("[aeiou_]","",factor.label,perl=T,ignore.case=T)))<25))
    addstr<-paste0('_SOM_',paste(gsub("[aeiou_]","",factor.label[1:print.length],
                                      perl=T,ignore.case=T),collapse='_'),
                   '_pl_',length(factor.label)-print.length)
  }
  #}}}
} else if (length(output.file)!=length(catalogues)) { 
  addstr<-paste0("_",loop.num)
}
output.ending<-rev(strsplit(output.file,'.',fixed=TRUE)[[1]])[1]
#}}}
#}}}

#Check if the output catalogue already exists {{{
if (grepl('.cat',output.file,fixed=T)) {
  #Output an LDAC catalogue {{{
  file =paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_DIRsom',addstr,'.cat'),output.file,fixed=TRUE))
  file2=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_refr_DIRsom',addstr,'.cat'),output.file,fixed=TRUE))
  #}}}
} else if (grepl('.Rdata',output.file,fixed=T)) {
  #Output an Rdata catalogue {{{
  file =paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_DIRsom',addstr,'.Rdata'),output.file,fixed=TRUE))
  file2=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_refr_DIRsom',addstr,'.Rdata'),output.file,fixed=TRUE))
  #}}}
} else { 
  #Output a CSV catalogue {{{
  file =paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_DIRsom',addstr,'.csv'),output.file,fixed=TRUE))
  file2=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_refr_DIRsom',addstr,'.csv'),output.file,fixed=TRUE))
  #}}}
}
if (file.exists(file) & file.exists(file2)) { 
  if (!force) { 
    cat("  ### Output Catalogues already exist! Not using [--force], so skipping! ###\n")
    next
  } else { 
    cat("  ### Output Catalogues already exist, but continuing because of --force ###\n")
  } 
}
#}}}

#Read in the input catalogue {{{
if (!quiet) { 
  cat(paste("  > Reading Training Catalogue")) 
  timer<-proc.time()[3]
}
if (grepl('.cat',train.catpath,fixed=T)) { 
  if (!quiet) { cat(" (with ldactoasc):\n\n") }
  train.cat<-try(read.ldac(train.catpath,ldactoasc=ldactoasc,data.table=FALSE,
                     diagnostic=TRUE ,options=paste('-s',keys,zt.label,count.variable.2,ldac.options.2),force=TRUE,clean=TRUE,showProgress=FALSE),silent=TRUE)
  if (suppressWarnings(class(train.cat)=='try-error')) { 
    cat(" - ")
    stop(paste0("Failure to read input catalogue. Check your ldactoasc binary?\n",
                "  NB: We assume that input \'.cat\' files are in LDAC format. Table details can be specified at input.\n",
                "  For other formats we use: .fits (FITS standard); .csv (CSV,TSV,ASCII); .Rdata (R data file)\n"))
  }
  if (!quiet) { cat("\n  -----> Completed!\n") }
} else if (grepl('.Rdata',train.catpath,fixed=T)) {
  nam<-try(load(train.catpath),silent=TRUE)
  if (suppressWarnings(class(nam)=='try-error')) { 
    cat(" - ")
    stop(paste0("Failure to load input Rdata catalogue.\n"))
  }
  if (nam=="nam") {
    #Whoops, we overwrote the catalogue while checking its name! {{{
    load(train.catpath,data.table=FALSE)
    train.cat<-nam
    rm('nam')
    #}}}
  } else if (nam!="train.cat") { 
    train.cat<-get(nam)
    rm(nam)
  } 
} else if (grepl('.fits',train.catpath,fixed=T)) {
  train.cat<-try(read.fits.cat(train.catpath,data.table=FALSE),silent=TRUE)
  if (suppressWarnings(class(train.cat)=='try-error')) { 
    cat(" - ")
    stop(paste0("Failure to read input fits catalogue.\n"))
  }
} else if (grepl('.csv',train.catpath,fixed=T) | grepl('.asc',refr.catpath,fixed=T)) {
  train.cat<-try(fread(train.catpath,data.table=FALSE,showProgress=FALSE),silent=TRUE)
  if (suppressWarnings(class(train.cat)=='try-error')) { 
    cat(" - ")
    stop(paste0("Failure to read input CSV catalogue.\n"))
  }
} else { 
  cat(" - ")
  stop(paste0("Failure to read input catalogue. Unknown file type. Expects .cat for ldac, .fits for binary, or .csv\n"))
}
#Make sure that the table is a data.table
if (!is.data.table(train.cat)){ 
  train.cat<-as.data.table(train.cat)
}
train.catnam<-rev(strsplit(train.catpath,'/')[[1]])[1]
train.ending<-rev(strsplit(train.catnam,'.',fixed=TRUE)[[1]])[1]
if (testing) {
  train.cat<-train.cat[runif(nrow(train.cat))<0.2,]
}
#Catalogue length

train.cat.len<-nrow(train.cat)
if (!quiet) { 
  io.clock<-io.clock+proc.time()[3]-timer
  cat(paste(" - Done",as.time(proc.time()[3]-timer,digits=0),"\n"))
  if (testing) { 
    cat("  # Running in Testing mode! Many parameters degraded to improve speed (--som.dim) #\n")
  }
}

#}}}

#Generate the DIR wieights {{{
#Check the number of HCs {{{
if (factor.nbins>=som.dim[1]*som.dim[2]) { 
  if (!quiet) { 
    cat("  # WARNING: Number of hierarchical clusters is >= number of SOM cells!           #\n")
    cat("  #          We will not use the hierarchical clustering, just split by SOM cell! #\n")
  }
  factor.nbins<-som.dim[1]*som.dim[2]
}
#}}}
#Notify {{{
if (!quiet) { 
  cat(paste("  > Generating the SOM from Training catalogue "))
  short.timer<-timer<-proc.time()[3]
}
#}}}
#Measure the SOM {{{
if (length(factor.label)==1) { 
  #There is a single parameter factor {{{
  #Notify {{{
  if (!quiet) { 
    cat(" (single numeric factor) {\n")
    cat(paste("    -> binning the factor in",factor.nbins,"chunks"))
  }#}}}
  #Bin the factor {{{
  cat.tab<-cut(train.cat[,eval(parse(text=factor.label))],breaks=seq(min(train.cat[,eval(parse(text=factor.label))]),max(train.cat[,eval(parse(text=factor.label))]),length=factor.nbins+1),include.lowest=TRUE)
  train.cat$GroupFactor<-as.numeric(cat.tab)
  cat.tab<-table(as.numeric(cat.tab))
  cat.tab.names<-names(cat.tab)
  if (any(grepl(" ",cat.tab.names))) {
    warning("Stripping whitespace from Factor names") 
    cat.tab.names<-gsub(" ","",cat.tab.names)
  }
  #}}}
  #}}}
} else { 
  #There are multiple factors: Cluster the data with self organising map {{{
  #Notify {{{
  if (!quiet) { 
    cat(" (multiple numeric factors) {\n")
    cat("    -> constructing scaled data vector")
  }#}}}
  #Check for character columns {{{
  seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",factor.label),' ')))
  if (any(sapply(train.cat[,seperated.labels,with=F],class)=='character')) { 
    if (!quiet) { cat(" (converting character cols to numeric!)") }
    for (i in sperated.labels[which(sapply(train.cat[1,seperated.labels,with=F],class)=='character')]) { 
      train.cat[,i,with=F]<-as.numeric(train.cat[,i,with=F])
    }
  }
  #}}}
  #Scale the data for use in the SOM {{{
  train.sc<-matrix(NA,nrow=train.cat.len,ncol=length(factor.label))
  refr.sc<-matrix(NA,nrow=refr.cat.len,ncol=length(factor.label))
  if (reuse && !file.exists(som.data.file)) {
    reuse<-FALSE
    warning("SOM data file does not exist! Cannot use previous SOM. Constructing a new one!")
  }
  if (rescale & !reuse) { 
    rescale.param<-matrix(NA,nrow=2,ncol=length(factor.label))
    colnames(rescale.param)<-factor.label
  } else if (rescale) { 
    load(file=som.data.file)
    rescale.param<-train.som$rescale.param
  }
  #Assign names 
  colnames(train.sc)<-colnames(refr.sc)<-factor.label
  #Loop through the factor expressions (they could be columns or expressions!)
  for (factor.expr in factor.label) { 
    #Calculate the expression result
    train.value<-train.cat[,eval(parse(text=factor.expr))]
    refr.value<-refr.cat[,eval(parse(text=factor.expr))]
    #seperate out the components
    seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",factor.expr),' ')))
    #Start with everything being OK, and then...
    train.not.detected<-train.missing<-rep(FALSE,train.cat.len)
    refr.not.detected<-refr.missing<-rep(FALSE,refr.cat.len)
    #...Loop through the components 
    for (label in seperated.labels) { 
      #Get the detected and non detected sources
      train.missing<-train.missing | (train.cat[[label]]==missing.val)
      refr.missing<-refr.missing | (refr.cat[[label]]==missing.val)
      train.not.detected<-train.not.detected | (train.cat[[label]]>=detect.limit & !train.missing)
      refr.not.detected<-refr.not.detected | (refr.cat[[label]]>=detect.limit & !refr.missing)
    }
    #Set the values with missing parts to dummy values 
    train.value[train.missing]<-NA
    refr.value[refr.missing]<-NA
    train.value[train.not.detected]<-NA
    refr.value[refr.not.detected]<-NA
    if (rescale & !reuse) { 
      #Calculate the detected median and mad
      med.tmp<-median(train.value,na.rm=T)
      mad.tmp<-mad(train.value,na.rm=T)
      #Scale the refr z and train z data consistently
      train.value<-(train.value-med.tmp)/mad.tmp
      refr.value<-(refr.value-med.tmp)/mad.tmp
      rescale.param[,factor.expr]<-c(med.tmp,mad.tmp)
    } else if (rescale) { 
      #Scale the refr z and train z data consistently
      train.value<-(train.value-rescale.param[1,factor.expr])/rescale.param[2,factor.expr]
      refr.value<-(refr.value-rescale.param[1,factor.expr])/rescale.param[2,factor.expr]
    } 
    train.sc[,factor.expr]<-train.value
    refr.sc[,factor.expr]<-refr.value
  }

  #Plot the resulting parameter distributions {{{ 
  if (plot>1) { 
    #Distributions before scaling {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_SOM_inputs.png',output.file,fixed=TRUE)),height=10*res,width=10*res,res=res)
    suppressWarnings(lay.mat<-matrix(c(1:length(factor.label),rep(0,100)),
                              ncol=4,nrow=floor(length(factor.label))/4+1,byrow=T))
    layout(lay.mat)
    par(mar=c(2,2,1,1))
    for (factor in factor.label) { 
      dens<-density(train.cat[,eval(parse(text=factor))],kernel='rect',na.rm=TRUE)
      magplot(dens,xlab=factor,ylab='')
      showbw(dens,kernel='rect',loc='topright',scale=0.2,col='blue')
      label('topleft',paste0('frac missing: ',round(length(which(!is.finite(train.cat[,eval(parse(text=factor))])))/length(train.cat[,eval(parse(text=factor))])*100,digits=2),"%"))
    }
    dev.off()
    #}}}
    #Distributions after scaling {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_SOM_scaled.png',output.file,fixed=TRUE)),height=10*res,width=10*res,res=res)
    layout(lay.mat)
    par(mar=c(2,2,1,1))
    for (factor in factor.label) { 
      dens<-density(train.sc[,factor],kernel='rect',na.rm=TRUE)
      magplot(dens,xlab=factor,ylab='')
      showbw(dens,kernel='rect',loc='topright',scale=0.2,col='blue')
      label('topleft',paste0('frac missing: ',round(length(which(!is.finite(train.cat[,eval(parse(text=factor))])))/length(train.cat[,eval(parse(text=factor))])*100,digits=2),"%"))
    }
    dev.off()
    #}}}
  }
  #}}}
  #Notify {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat(paste0("\n    -> constructing som grid (",som.dim[1],'x',som.dim[2],')'))
  }#}}}
  #}}}
  #Create the SOM grid {{{
  data.grid<-somgrid(xdim = som.dim[1], ydim=som.dim[2], topo=som.topo, toroidal=som.toroidal)
  #}}}
  #Generate the SOM {{{
  #if (reuse) { 
  if (reuse && file.exists(som.data.file)) {
    #Notify {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> loading previously calculated SOM (!!)")
    }#}}}
    #Load the SOM {{{
    name<-load(som.data.file)
    if (name!="train.som") { 
      train.som<-get(name)
    }
    #}}}
    #Notify {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> passing training set into previous SOM (!!)")
    }#}}}
    #Get the data positions from the trained SOM {{{
    if (som.cores<1) {
      num_splits<-try(as.numeric(system("grep -c ^processor /proc/cpuinfo",intern=T)))
    } else { 
      num_splits<-som.cores
    }
    if (class(num_splits)!='try-error') { 
      require(foreach)
      require(doParallel)
      registerDoParallel(cores=num_splits)
      train.som$training.classif<-train.som$unit.classif
      if (!quiet) { 
        cat(paste(' [parallel:',round(nrow(train.sc)/num_splits),'gal/core]'))
      }
      train.som$unit.classif <-
       foreach(d=isplitRows(train.sc, chunks=num_splits),
       .combine=c, .packages=c("stats")) %dopar% {
          return=predict(train.som, newdata=d)$unit.classif
      }
    } else { 
      if (!quiet) { 
        cat(paste(' [serial: 1 core only]'))
      }
      pred.som<-predict(train.som,newdata=train.sc)
      train.som$training.classif<-train.som$unit.classif
      train.som$unit.classif<-pred.som$unit.classif
    } 
    #}}}
  } else { 
    if (reuse) { 
      warning("SOM data file does not exist! Cannot use previous SOM. Constructing a new one!")
      reuse<-FALSE
    }
    if (sparse.som) { 
      #Notify {{{
      if (!quiet) { 
        cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
        short.timer<-proc.time()[3]
      }#}}}
      #Calculate the SOM using a sparse sampling of the data {{{
      if (sparse.frac*nrow(train.sc)<som.dim[1]*som.dim[2]*sparse.min.num) { 
        if (!quiet) { 
          cat(paste0("- WARNING: Sparse sampling creates fewer than ",sparse.min.num," galaxies per SOM cell!\n"))
          cat(paste0("       |-> forcing at least ",sparse.min.num," galaxies per cell.\n"))
        }
        sparse.frac<-som.dim[1]*som.dim[2]*sparse.min.num/nrow(train.sc)
        if (!quiet) { 
          if (sparse.frac > 1) { 
            cat("           |-> ERROR: the SOM grid is too fine to sparse sample. Using full data vector\n")
            sparse.frac<-1
          } else { 
            cat("    -> continuing constructing SOM from sparse sampling of data vector")
          }
        }
      }
      if (exists("sparse.var") && length(train.cat[[sparse.var]])!=0) { 
        cat(paste("\n    -> constructing SOM from sampling of",sparse.var,"in data vector"))
        sparse.index<-sample(nrow(train.sc),size=ceiling(sparse.frac*nrow(train.sc)),prob=train.cat[[sparse.var]])
      } else {
        cat("\n    -> constructing SOM from sparse sampling of data vector")
        sparse.index<-sample(nrow(train.sc),size=ceiling(sparse.frac*nrow(train.sc)))
      }
      #Distributions after scaling {{{
      if (plot>1) { 
        png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_SOM_scaled_sparse.png',output.file,fixed=TRUE)),height=10*res,width=10*res,res=res)
        layout(lay.mat)
        par(mar=c(2,2,1,1))
        for (factor in factor.label) { 
          dens<-density(train.sc[,factor],kernel='rect',na.rm=TRUE)
          magplot(dens,xlab=factor,ylab='')
          lines(density(train.sc[sparse.index,factor],kernel='rect',bw=dens$bw,na.rm=TRUE),col='red')
          showbw(dens,kernel='rect',loc='topright',scale=0.2,col='blue')
      label('topleft',paste0('frac missing: ',round(length(which(!is.finite(train.cat[,eval(parse(text=factor))])))/length(train.cat[,eval(parse(text=factor))])*100,digits=2),"%"))
        }
        dev.off()
        #}}}
      }
      train.som<-(som(train.sc[sparse.index,], grid=data.grid, rlen=som.iter, alpha=som.rate, cores=som.cores, mode=som.method,maxNA=maxNAfrac))
      if (class(train.som)=='try-error') { 
        #Notify & Loop {{{
        cat("Error: SOM Training Failed! Skipping\n")
        if(!quiet & loop.num!=loop.length) { 
          cat(paste("Loop ERROR. Elapsed time:",as.time(proc.time()[3]-start.timer,digits=0),
                    "; Approximate time remaining:",
            as.time((proc.time()[3]-start.timer)/(loop.num-loop.start+1)*(loop.length-loop.num),digits=0),"\n")) 
          next
        }
        #}}}
      }
      if (rescale) { 
        train.som$rescale.param<-rescale.param
      }
      #Save the SOM to file {{{
      save(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0(addstr,'_SOMdata.Rdata'),output.file,fixed=TRUE)),train.som)
      #}}}
      #Skip the rest of the measurements?! {{{
      if (only.som) { 
        #Notify {{{
        if (!quiet) { 
          cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
          short.timer<-proc.time()[3]
          cat("\n    -> Skipping remaining loop because '--only.som' flag was set!!")
          next
        }#}}}
      }
      #}}}
      #Notify {{{
      if (!quiet) { 
        cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
        short.timer<-proc.time()[3]
        cat("\n    -> dispersing data into trained SOM")
      }#}}}
      #Get the data positions from the trained SOM {{{
      if (som.cores<1) {
        num_splits<-try(as.numeric(system("grep -c ^processor /proc/cpuinfo",intern=T)))
      } else { 
        num_splits<-som.cores
      }
      if (class(num_splits)!='try-error') { 
        require(foreach)
        require(doParallel)
        registerDoParallel(cores=num_splits)
        if (!quiet) { 
          cat(paste(' [parallel:',round(nrow(train.sc)/num_splits),'gal/core]'))
        }
        train.som$training.classif<-train.som$unit.classif
        train.som$unit.classif <-
         foreach(d=isplitRows(train.sc, chunks=num_splits),
         .combine=c, .packages=c("stats")) %dopar% {
            return=predict(train.som, newdata=d)$unit.classif
        }
      } else { 
        if (!quiet) { 
          cat(paste(' [serial: 1 core only]'))
        }
        pred.som<-predict(train.som,newdata=train.sc)
        train.som$training.classif<-train.som$unit.classif
        train.som$unit.classif<-pred.som$unit.classif
      } 
      #}}}
      #}}}
    } else { 
      #Notify {{{
      if (!quiet) { 
        cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
        short.timer<-proc.time()[3]
        cat("\n    -> constructing SOM from full data vector")
      }#}}}
      #Calculate the SOM using the full data vector {{{
      train.som<-(som(train.sc, grid=data.grid, rlen=som.iter, alpha=som.rate, cores=som.cores, mode=som.method,maxNA=maxNAfrac))
      if (class(train.som)=='try-error') { 
        #Notify & Loop {{{
        cat("Error: SOM Training Failed! Skipping\n")
        if(!quiet & loop.num!=loop.length) { 
          cat(paste("Loop ERROR. Elapsed time:",as.time(proc.time()[3]-start.timer,digits=0),
                    "; Approximate time remaining:",
            as.time((proc.time()[3]-start.timer)/(loop.num-loop.start+1)*(loop.length-loop.num),digits=0),"\n")) 
          next
        }
        #}}}
      }
      if (rescale) { 
        train.som$rescale.param<-rescale.param
      }
      #Save the SOM to file {{{
      save(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0(addstr,'_SOMdata.Rdata'),output.file,fixed=TRUE)),train.som)
      #}}}
      #Skip the rest of the measurements?! {{{
      if (only.som) { 
        #Notify {{{
        if (!quiet) { 
          cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
          short.timer<-proc.time()[3]
          cat("\n    -> Skipping remaining loop because '--only.som' flag was set!!")
          next
        }#}}}
      }
      #}}}
      #}}}
    } 
  } 
  #Notify {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat(paste("\n    -> constructing",factor.nbins,"hierarchical clusters"))
  }#}}}
  #}}}
  #Split into factor.nbins hierarchical clusters {{{
  if (factor.nbins!=som.dim[1]*som.dim[2]) { 
    train.hc = cutree(hclust(dist(train.som$codes[[1]])), factor.nbins)
  } else { 
    train.hc = 1:factor.nbins
  } 
  #}}}
  #Get the makeup of each of the HCs, in terms of components {{{
  cluster.spread<-cluster.median<-matrix(NA,nrow=factor.nbins,ncol=ncol(train.som$codes[[1]]))
  for (i in 1:factor.nbins) { 
    if (length(which(train.hc==i))>=1) { 
      cluster.median[i,]<-colMedians(train.som$codes[[1]][which(train.hc==i),,drop=FALSE])
      cluster.spread[i,]<-colMads(train.som$codes[[1]][which(train.hc==i),,drop=FALSE])
    }
  }
  colnames(cluster.median)<-colnames(cluster.spread)<-colnames(train.som$codes[[1]])
  #}}}
  #Plot the SOM and additional information {{{
  if (plot>0) { 
    #Notify {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> Plotting SOM figures")
    }#}}}
    #Plot the som codes {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_SOM_codes.png',output.file,fixed=TRUE)),height=10*res,width=10*res,res=res)
    plot(train.som, type="codes", bgcol=rainbow(factor.nbins)[train.hc],shape='straight',codeRendering='lines',border=NA)
    add.cluster.boundaries(train.som,train.hc,lwd=3,col=hsv(1,0,0,0.8))
    dev.off()
    #}}}
    #Plot the cluster makeup {{{
    #Plot the changes at each iteration {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_SOM_changes.png',output.file,fixed=TRUE)),height=5*res,width=5*res,res=res)
    plot(train.som, type="changes")
    dev.off()
    #}}}
    if (plot>1) {
      #Plot the various SOM components {{{
      png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_SOM_props%02d.png',output.file,fixed=TRUE)),height=5*res,width=5*res,res=res)
      #Set the margins {{{
      par(mar=c(0,0,0,0))
      #}}}
      ##Define the plot layout {{{
      #suppressWarnings(lay.mat<-matrix(c(1,2,3,4,c(4+1:length(train.som$codes[[1]][1,]),rep(0,100))),
      #                          ncol=4,nrow=floor(length(train.som$codes[[1]][1,])/4)+2,byrow=T))
      ##Remove any useless rows {{{
      #if (any(rowSums(lay.mat)==0)) { 
      #  lay.mat<-lay.mat[-which(rowSums(lay.mat)==0),]
      #}#}}}
      #layout(lay.mat)
      ##}}}
      #Get the colour palette {{{
      BlRd<-function(n,alpha=1){rainbow(n,end=2/3,alpha=alpha)[n:1]}
      #}}}
      #Plot the counts per SOM cell {{{
      cell.counts<-plot(train.som, type="count",shape='straight',border=NA,heatkeywidth=1,zlog=TRUE)
      add.cluster.boundaries(train.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
      #}}}
      #Plot the neighbour distances per SOM cell {{{
      plot(train.som, type="dist.neighbours",shape='straight',border=NA,heatkeywidth=1)
      add.cluster.boundaries(train.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
      #}}}
      #Plot the quality per SOM cell {{{
      plot(train.som, type="quality",shape='straight',border=NA,heatkeywidth=1)
      add.cluster.boundaries(train.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
      #}}}
      #Plot the number density per SOM cell {{{
      plot(train.som, type="count",shape='straight',border=NA,heatkeywidth=1)
      add.cluster.boundaries(train.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
      #}}}
      for (i in 1:length(train.som$codes[[1]][1,])) { 
        #Plot the i^th data property per SOM cell {{{
        plot(train.som, type = "property", property = train.som$codes[[1]][,i],
             main=paste0(colnames(train.som$data[[1]])[i],'(scaled)'), palette.name=BlRd,shape='straight',border=NA,heatkeywidth=1)
        add.cluster.boundaries(train.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
        #}}}
      }
      #Close the plot device {{{
      dev.off()
      #}}}
      #}}}
    }
    #}}}
  }
  #}}}
  #Notify {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Dispersing data into cluster bins")
  }#}}}
  #Sort the data into the heirarchical clusters {{{
  #Get the individual data-to-SOMcell IDs {{{
  somclust<-train.som$unit.classif
  #}}}
  #Convert the IDs from data-to-SOMcell into data-to-cluster {{{
  if (factor.nbins!=som.dim[1]*som.dim[2]) { 
    for (i in 1:factor.nbins) { 
      #Assign the data to the cluster if it's SOMcell belongs to the cluster
      somclust[train.som$unit.classif%in%which(train.hc==i)]<-i 
    } 
  }
  #}}}
  cat.tab<-table(somclust)
  train.cat$GroupFactor<-somclust
  #}}}
  #Get the data-cluster names {{{
  cat.tab.names<-names(cat.tab)
  if (any(grepl(" ",cat.tab.names))) {
    warning("Stripping whitespace from Factor names") 
    cat.tab.names<-gsub(" ","",cat.tab.names)
  }
  #}}}
  #}}}
} 
#}}}
#Notify {{{
if (!quiet) { 
  cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
  short.timer<-proc.time()[3]
  cat("\n    -> Passing reference data vector into SOM")
}#}}}
#Get the data positions from the trained SOM {{{
refr.som<-train.som
if (som.cores<1) {
  num_splits<-try(as.numeric(system("grep -c ^processor /proc/cpuinfo",intern=T)))
} else { 
  num_splits<-som.cores
}
if (class(num_splits)!='try-error') { 
  require(foreach)
  require(doParallel)
  registerDoParallel(cores=num_splits)
  if (!quiet) { 
    cat(paste(' [parallel:',round(nrow(refr.sc)/num_splits),'gal/core]'))
  }
  refr.som$unit.classif <-
   foreach(d=isplitRows(refr.sc, chunks=num_splits),
   .combine=c, .packages=c("stats")) %dopar% {
      return=predict(train.som, newdata=d)$unit.classif
  }
} else { 
  pred.som<-predict(train.som,newdata=refr.sc)
  refr.som$training.classif<-refr.som$unit.classif
  refr.som$unit.classif<-pred.som$unit.classif
} 
save(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0(addstr,'_refr_SOMdata.Rdata'),output.file,fixed=TRUE)),refr.som)
#}}}
#Plot the SOM and additional information {{{
if (plot>0) { 
  #Notify {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Plotting SOM figures")
  }#}}}
  #Plot the som codes {{{
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_refr_SOM_codes.png',output.file,fixed=TRUE)),height=10*res,width=10*res,res=res)
  plot(refr.som, type="codes", bgcol=rainbow(factor.nbins)[train.hc],shape='straight',codeRendering='lines',border=NA)
  add.cluster.boundaries(refr.som,train.hc,lwd=3,col=hsv(1,0,0,0.8))
  dev.off()
  #}}}
  #Plot the changes at each iteration {{{
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_refr_SOM_changes.png',output.file,fixed=TRUE)),height=5*res,width=5*res,res=res)
  plot(refr.som, type="changes")
  dev.off()
  #}}}
  if (plot>1) { 
    #Plot the various SOM components {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_refr_SOM_props%02d.png',output.file,fixed=TRUE)),height=5*res,width=12*res,res=res)
    #Set the margins {{{
    par(mar=c(0,0,0,0))
    #}}}
    ##Define the plot layout {{{
    #suppressWarnings(lay.mat<-matrix(c(1,2,3,4,c(4+1:length(refr.som$codes[[1]][1,]),rep(0,100))),
    #                          ncol=4,nrow=floor(length(refr.som$codes[[1]][1,])/4)+2,byrow=T))
    ##Remove any useless rows {{{
    #if (any(rowSums(lay.mat)==0)) { 
    #  lay.mat<-lay.mat[-which(rowSums(lay.mat)==0),]
    #}#}}}
    #layout(lay.mat)
    #}}}
    #Plot the counts per SOM cell {{{
    cell.counts<-plot(refr.som, type="count",shape='straight',border=NA,heatkeywidth=1,zlog=TRUE)
    add.cluster.boundaries(refr.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
    #}}}
    #Plot the neighbour distances per SOM cell {{{
    plot(refr.som, type="dist.neighbours",shape='straight',border=NA,heatkeywidth=1)
    add.cluster.boundaries(refr.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
    #}}}
    #Plot the quality per SOM cell {{{
    plot(refr.som, type="quality",shape='straight',border=NA,heatkeywidth=1)
    add.cluster.boundaries(refr.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
    #}}}
    #Plot the number density per SOM cell {{{
    plot(refr.som, type="count",shape='straight',border=NA,heatkeywidth=1)
    add.cluster.boundaries(refr.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
    #}}}
    for (i in 1:length(refr.som$codes[[1]][1,])) { 
      #Plot the i^th data property per SOM cell {{{
      plot(refr.som, type = "property", property = refr.som$codes[[1]][,i],
           main=paste0(colnames(refr.som$data[[1]])[i],'(scaled)'), palette.name=BlRd,shape='straight',border=NA,heatkeywidth=1)
      add.cluster.boundaries(refr.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
      #}}}
    }
    #Close the plot device {{{
    dev.off()
    #}}}
    #}}}
  }
}
#}}}
#Notify {{{
if (!quiet) { 
  cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
  short.timer<-proc.time()[3]
  cat("\n    -> Sorting Reference Data into HCs ")
}
#}}}
#Sort the refr data into the heirarchical clusters {{{
#Get the individual data-to-SOMcell IDs {{{
somclust<-refr.som$unit.classif
#}}}
#Convert the IDs from data-to-SOMcell into data-to-cluster {{{
if (factor.nbins!=som.dim[1]*som.dim[2]) { 
  require(foreach)
  require(doParallel)
  registerDoParallel(cores=num_splits)
  somind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=FALSE)%dopar%{ 
    #Assign the data to the cluster if it's SOMcell belongs to the cluster
    ind<-which(refr.som$unit.classif%in%which(train.hc==i))
    return=cbind(ind,rep(i,length(ind)))
  } 
  somclust[somind[,1]]<-somind[,2]
}
#}}}
cat.tab<-table(somclust)
refr.cat$GroupFactor<-somclust
#}}}
#Notify {{{
if (!quiet) { 
  cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
  short.timer<-proc.time()[3]
  cat("\n    -> Generating the SOM DIR weights ")
}
#}}}
#Generate the weights {{{
#train.cat$SOMweight<-NA
#refr.cat$SOMweight<-NA
#if (count.variable.1=='' & count.variable.2=='') { 
#  for (i in 1:factor.nbins) { 
#    train.index<-which(train.cat$GroupFactor==i)
#    refr.index<-which(refr.cat$GroupFactor==i)
#    train.cat[["SOMweight"]][train.index]<-length(refr.index)/length(train.index)
#    refr.cat[["SOMweight"]][refr.index]<-length(train.index)/length(refr.index)
#  }
#} else if (count.variable.1=='') { 
#  for (i in 1:factor.nbins) { 
#    train.index<-which(train.cat$GroupFactor==i)
#    refr.index<-which(refr.cat$GroupFactor==i)
#    train.cat[["SOMweight"]][train.index]<-length(refr.index)/sum(train.cat[[count.variable.2]][train.index])
#    refr.cat[["SOMweight"]][refr.index]<-sum(train.cat[[count.variable.2]][train.index],na.rm=T)/length(refr.index)
#  }
#} else if (count.variable.2=='') { 
#  for (i in 1:factor.nbins) { 
#    train.index<-which(train.cat$GroupFactor==i)
#    refr.index<-which(refr.cat$GroupFactor==i)
#    train.cat[["SOMweight"]][train.index]<-sum(refr.cat[[count.variable.1]][refr.index],na.rm=T)/length(train.index)
#    refr.cat[["SOMweight"]][refr.index]<-length(train.index)/sum(refr.cat[[count.variable.1]][refr.index])
#  }
#} else { 
#  for (i in 1:factor.nbins) { 
#    train.index<-which(train.cat$GroupFactor==i)
#    refr.index<-which(refr.cat$GroupFactor==i)
#    train.cat[["SOMweight"]][train.index]<-sum(refr.cat[[count.variable.1]][refr.index],na.rm=T)/sum(train.cat[[count.variable.2]][train.index])
#    refr.cat[["SOMweight"]][refr.index]<-sum(train.cat[[count.variable.2]][train.index],na.rm=T)/sum(refr.cat[[count.variable.1]][refr.index])
#  }
#} 
if (count.variable.1=='' & count.variable.2=='') { 
  wtind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=TRUE)%dopar%{
    train.index<-which(train.cat$GroupFactor==i)
    refr.index<-which(refr.cat$GroupFactor==i)
    wt<-length(refr.index)/length(train.index)
    return=cbind(i,wt)
  }
} else if (count.variable.1=='') { 
  wtind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=TRUE)%dopar%{
    train.index<-which(train.cat$GroupFactor==i)
    refr.index<-which(refr.cat$GroupFactor==i)
    wt<-length(refr.index)/sum(train.cat[[count.variable.2]][train.index])
    return=cbind(i,wt)
  }
} else if (count.variable.2=='') { 
  wtind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=TRUE)%dopar%{
    train.index<-which(train.cat$GroupFactor==i)
    refr.index<-which(refr.cat$GroupFactor==i)
    wt<-sum(refr.cat[[count.variable.1]][refr.index],na.rm=T)/length(train.index)
    return=cbind(i,wt)
  }
} else { 
  wtind<-foreach(i=1:factor.nbins,.combine=rbind,.inorder=TRUE)%dopar%{
    train.index<-which(train.cat$GroupFactor==i)
    refr.index<-which(refr.cat$GroupFactor==i)
    wt<-sum(refr.cat[[count.variable.1]][refr.index],na.rm=T)/sum(train.cat[[count.variable.2]][train.index])
    return=cbind(i,wt)
  }
} 
train.cat[["SOMweight"]]<-wtind[train.cat$GroupFactor,2]
refr.cat[["SOMweight"]]<-1/wtind[refr.cat$GroupFactor,2]
#train.cat$SOMweight<-train.cat$SOMweight/sum(train.cat$SOMweight,na.rm=TRUE)
#Plot the SOM weights {{{
if (plot>0) {
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_train_SOMweights.png',output.file,fixed=TRUE)),height=5*res,width=5*res,res=res)
  dens<-density(log10(train.cat$SOMweight),kernel='rect',from=-3,to=3,bw=0.1/sqrt(12),na.rm=TRUE)
  magplot(dens,xlab='SOMweight',ylab='PDF',unlog='x')
  dev.off()
}
#}}}
#Plot the SOM weights {{{
if (plot>0) {
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_refr_SOMweights.png',output.file,fixed=TRUE)),height=5*res,width=5*res,res=res)
  dens<-density(log10(refr.cat$SOMweight),from=-5,to=1,bw=0.1/sqrt(12),kernel='rect',na.rm=TRUE)
  magplot(dens,xlab='SOMweight',ylab='PDF',unlog='x')
  dev.off()
}
#Paint the SOM weights {{{
if (plot>0) {
  wrefr.cell<-wtrain.cell<-rep(NA,som.dim[1]*som.dim[2])
  for (i in 1:factor.nbins) { 
    wtrain.cell[which(train.hc==i)]<-train.cat$SOMweight[which(train.cat$GroupFactor==i)[1]]
    wrefr.cell[which(train.hc==i)]<-refr.cat$SOMweight[which(refr.cat$GroupFactor==i)[1]]
  }
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_SOMweights_paint.png',output.file,fixed=TRUE)),height=5*res,width=7*res,res=res)
  layout(cbind(1,2))
  plot(refr.som, type = "property", property = wrefr.cell,#zlim=c(0,1.5),
       main="SOMweight_refr", palette.name=scale.palette(wrefr.cell,palette=BlRd),shape='straight',border=NA,heatkeywidth=som.dim[1]/20,ncol=1e3,heatkeyborder=NA)
  plot(train.som, type = "property", property = wtrain.cell,#zlim=c(0,1.5),
       main="SOMweight_train", palette.name=scale.palette(wtrain.cell,palette=BlRd),shape='straight',border=NA,heatkeywidth=som.dim[1]/20,ncol=1e3,heatkeyborder=NA)
  dev.off()
}
#}}}
#}}}
#}}}
#Generate the reference catalogue plots {{{ 
if (plot>0) { 
  #Notify {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    short.timer<-proc.time()[3]
    cat("\n    -> Plotting the Reference data ")
  }
  #}}}
  #Plot the counts zoomed in {{{
  png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_counts.png',output.file,fixed=TRUE)),height=5*res,width=10*res,res=res)
  layout(cbind(1,2))
  train.count.tab<-table(train.som$unit.classif)
  refr.count.tab<-table(refr.som$unit.classif)
  plot(train.som, type="counts",shape='straight',border=NA,heatkeywidth=som.dim[1]/20,main='Training Sample',palette=scale.palette(log10(train.count.tab)),ncol=1e3,heatkeyborder=NA,zlog=TRUE)
  #plot(train.som, type="counts",shape='straight',border=NA,heatkeywidth=som.dim[1]/20,main='Training Sample',palette=scale.palette(train.count.tab),ncol=1e3,heatkeyborder=NA)
  #     #zlim=quantile(train.count.tab,probs=c(0.1,0.9),na.rm=TRUE))
  add.cluster.boundaries(train.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
  plot(refr.som, type="counts",shape='straight',border=NA,heatkeywidth=som.dim[1]/20,main='Reference Sample',palette=scale.palette(log10(refr.count.tab)),ncol=1e3,heatkeyborder=NA,zlog=TRUE)
  #plot(refr.som, type="counts",shape='straight',border=NA,heatkeywidth=som.dim[1]/20,main='Reference Sample',palette=scale.palette(refr.count.tab),ncol=1e3,heatkeyborder=NA)
  #     #zlim=quantile(refr.count.tab,probs=c(0.1,0.9),na.rm=TRUE))
  add.cluster.boundaries(refr.som,train.hc,lwd=1,col=hsv(1,0,0,0.3))
  dev.off()
  #}}}
  #Measure the per-cell zrefr and z_B {{{
  zrefr.cell<-ztrain.cell<-zrefr.sd.cell<-ztrain.sd.cell<-rep(NA,som.dim[1]*som.dim[2])
  for (i in 1:factor.nbins) { 
    #Assign the data to the cluster if it's SOMcell belongs to the cluster
    if (length(which(refr.cat$GroupFactor==i))!=0) {
      zrefr.cell[which(train.hc==i)]<-median(refr.cat[[zr.label]][which(refr.cat$GroupFactor==i)],na.rm=TRUE) 
      zrefr.sd.cell[which(train.hc==i)]<-mad(refr.cat[[zr.label]][which(refr.cat$GroupFactor==i)],na.rm=TRUE) 
    }
    if (length(which(train.cat$GroupFactor==i))!=0) {
      ztrain.cell[which(train.hc==i)]<-median(train.cat[[zt.label]][which(train.cat$GroupFactor==i)],na.rm=TRUE) 
      ztrain.sd.cell[which(train.hc==i)]<-mad(train.cat[[zt.label]][which(train.cat$GroupFactor==i)],na.rm=TRUE) 
    }
  }
  #}}}
  if (plot>1) { 
    #Plot the SEDs of each cluster {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_clustSEDs.png',output.file,fixed=TRUE)),height=10*res,width=20*res,res=res)
    magplot(NA,xlim=c(0,length(factor.label)),ylim=range(cluster.median-colMedians(cluster.median)),xlab='Factor Index',ylab='Factor value + Cluster Index')
    count<-1
    for (i in order(cluster.median[,1])) { 
      lines(1:length(factor.label),cluster.median[i,]-colMedians(cluster.median),col=hsv(2/3*count/factor.nbins,alpha=0.3),lwd=2)
      count<-count+1
    }
    dev.off()
    #}}}
    #Plot the inferred z distributions {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_zz.png',output.file,fixed=TRUE)),height=5*res,width=10*res,res=res)
    layout(cbind(1,2))
    par(mar=c(2,2,0,0),oma=c(2,2,1,1))
    #Generate the implied z distributions 
    syn.refr.z<-rep(NA,refr.cat.len)
    syn.train.z<-rep(NA,train.cat.len)
    for (i in 1:factor.nbins) { 
      refr.ind<-which(refr.cat$GroupFactor==i)
      train.ind<-which(train.cat$GroupFactor==i)
      syn.refr.z[refr.ind]<-sample(refr.cat[[zr.label]][refr.ind],size=length(refr.ind),replace=T)
      syn.train.z[train.ind]<-sample(train.cat[[zt.label]][train.ind],size=length(train.ind),replace=T)
    }
    hist2D(train.cat[[zt.label]],syn.train.z,nbins=101,zlog=T,xlim=c(0,2),ylim=c(0,2),xlab='z input',ylab='z sampled',main='Training Sample',asp=0,flip=T)
    hist2D(refr.cat[[zr.label]],syn.refr.z,nbins=101,zlog=T,xlim=c(0,2),ylim=c(0,2),xlab='z input',ylab='z sampled',main='Reference Sample',asp=0,flip=T)
    dev.off()
    #}}}
    #Plot the z distributions zoomed in {{{
    png(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_zpaint.png',output.file,fixed=TRUE)),height=10*res,width=10*res,res=res)
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
    #}}}
  }
  #Plot the effect of the SOM weights {{{
  pdf(file=paste0(output.path,'/',sub(paste0('.',output.ending),'_wgt_colours.pdf',output.file,fixed=TRUE)),height=5,width=10)
  par(mar=c(2,2,2,2),oma=c(3,3,3,1))
  suppressWarnings(lay.mat<-matrix(c(1:(length(factor.label)+1),rep(0,100)),
                            ncol=4,nrow=min(c(floor(length(factor.label)/4)+2,3)),byrow=T))
  #Remove any useless rows {{{
  if (any(rowSums(lay.mat)==0)) { 
    lay.mat<-lay.mat[-which(rowSums(lay.mat)==0),]
  }#}}}

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
    #}}}
  }
}
#}}}
#}}}

#Output the DIR weighted catalogues {{{
if(!quiet) { 
  cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
  cat(paste("\n  # DIR_weights Training catalogue name:",sub(paste0('.',output.ending),paste0('_DIRsom',addstr,'.',output.ending),output.file,fixed=TRUE),'\n')) 
  cat(paste("  > Outputing DIR weights catalogue")) 
  timer<-proc.time()[3]
}
if (grepl('.cat',output.file,fixed=T)) {
  #Output an LDAC catalogue {{{
  if (is.data.table(train.cat)) { 
    train.cat<-as.data.frame(train.cat)
  }
  write.ldac(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_DIRsom',addstr,'.cat'),output.file,fixed=TRUE)),train.cat,force=TRUE,clean=TRUE,
  asctoldac=sub('ldactoasc','asctoldac',ldactoasc))
  #}}}
} else if (grepl('.Rdata',output.file,fixed=T)) {
  #Output an Rdata catalogue {{{
  save(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_DIRsom',addstr,'.Rdata'),output.file,fixed=TRUE)),train.cat)
  #}}}
} else { 
  #Output a CSV catalogue {{{
  write.csv(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_DIRsom',addstr,'.csv'),output.file,fixed=TRUE)),train.cat,quote=F,row.names=F)
  #}}}
}
if(!quiet) { 
  cat(paste("\n  # DIR_weights Reference catalogue name:",sub(paste0('.',output.ending),paste0('_refr_DIRsom',addstr,'.',output.ending),output.file,fixed=TRUE),'\n')) 
  cat(paste("  > Outputing DIR weights catalogue")) 
  timer<-proc.time()[3]
}
if (grepl('.cat',output.file,fixed=T)) {
  #Output an LDAC catalogue {{{
  if (is.data.table(refr.cat)) { 
    refr.cat<-as.data.frame(refr.cat)
  }
  write.ldac(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_refr_DIRsom',addstr,'.cat'),output.file,fixed=TRUE)),refr.cat,force=TRUE,clean=TRUE,
  asctoldac=sub('ldactoasc','asctoldac',ldactoasc))
  #}}}
} else if (grepl('.Rdata',output.file,fixed=T)) {
  #Output an Rdata catalogue {{{
  save(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_refr_DIRsom',addstr,'.Rdata'),output.file,fixed=TRUE)),refr.cat)
  #}}}
} else { 
  #Output a CSV catalogue {{{
  write.csv(file=paste0(output.path,'/',sub(paste0('.',output.ending),paste0('_refr_DIRsom',addstr,'.csv'),output.file,fixed=TRUE)),refr.cat,quote=F,row.names=F)
  #}}}
}
#Notify {{{
if(!quiet) { 
  io.clock<-io.clock+proc.time()[3]-timer
  cat(paste(" - Done",as.time(proc.time()[3]-timer,digits=0),"\n")) 
}
#}}}
#}}}

#Notify & Loop {{{
if(!quiet & loop.num!=loop.length) { 
  cat(paste("Loop completed. Elapsed time:",as.time(proc.time()[3]-start.timer,digits=0),
            "; Approximate time remaining:",
            as.time((proc.time()[3]-start.timer)/(loop.num-loop.start+1)*(loop.length-loop.num),digits=0),"\n")) 
}
#}}}
#}}}
}

#Notify and End {{{
if(!quiet) { cat(paste("DIR weight Generation completed in",as.time(proc.time()[3]-start.timer,digits=0),"(",as.time(io.clock,digits=0),"catalogue IO)\n\n")) }
#}}}

