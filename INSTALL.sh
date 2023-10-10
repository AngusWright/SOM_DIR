#
# Install script for the SOM_DIR repository 
# Written by A.H. Wright (2018-10-22)
#

#Attempt to get the paths of R & Rscript from $PATH
R=`which R`
Rscript=`which Rscript`

#If R doesn't exist, it should be installed 
if [ "$R" == "" ]; 
then
  echo "R must be manually installed on your machine & placed in your \$PATH"
  exit 1
fi 

#If Rscript doesn't exist, it should be installed 
if [ "$Rscript" == "" ]; 
then
  echo "Rscript must be manually installed on your machine & placed in your \$PATH"
  echo "NB: You appear to have R, so Rscript is probably just missing from your \$PATH"
  exit 1
fi 

#Check that the relevant R libraries are installed 
$R --no-save --slave  <<EOF 

#Function to quickly check for and install libraries
require.and.load<-function(name,githubrep) { 
  if (!missing(githubrep)) { 
    remotes::install_github(paste(githubrep,name,sep='/'),upgrade='always')
    if (grepl('/',name)) {
      name<-rev(strsplit(name,'/')[[1]])[1]
    }
  } else { 
    install.packages(name,repos='https://cloud.r-project.org/')
  }
  if (!require(name,character.only=TRUE)) { 
    stop(paste("Failed to install package",name))
  }
}

#Look for and install the important packages for the script
require.and.load('remotes') 
require.and.load('RColorBrewer') 
require.and.load('kohonen') 
require.and.load('KernSmooth') 
require.and.load('itertools') 
require.and.load('matrixStats')
require.and.load('data.table') 
require.and.load('helpRfuncs','AngusWright') 
require.and.load('FITSio') 
#require.and.load('LAMBDAR','AngusWright') 
require.and.load('kohonen/kohonen','AngusWright') 
EOF


