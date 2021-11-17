

#========================================================================================#
#  ======== CONNECT TO CLUSTER  =========================================================#
#========================================================================================#

ssh -o ConnectTimeout=10 ucjucho@myriad.rc.ucl.ac.uk

# =======Define Home directory ==============
HOME="/lustre/scratch/scratch/ucjucho/TS2021_CommonLiabAddiction"
# ======= Define local directory ==============
LOCAL="githubProjects/TabeaSchoeler/TS2021_CommonLiabAddiction"
# ======= Pathname where all gwas sumstats are stored  ==============
gwasDir="/lustre/scratch/scratch/ucjucho/GWAS_SumStats"

# =========== LOAD REQUIRED SOFTWARES ================
# Load R module
module unload -f compilers mpi gcc-libs
module load beta-modules
module load r/r-4.0.2_bc-3.11
# Load plink and python
module load plink/1.90b3.40
module load java/1.8.0_92


# Create script containing required libraries and directories
cat > $HOME/analysis/input.R <<EOF
gwasDir='$gwasDir'
HOME='$HOME'
LOCAL='$LOCAL'
R_libPaths='$HOME/programs/R'
.libPaths(R_libPaths)
load.lib=c('devtools', 'dplyr','GenomicSEM', 'rdrop2', 'data.table','R.utils', 'parallel', 'phenoscanner',  'biomaRt', 'gprofiler2', 'VariantAnnotation', 'plyr', 'CMplot', 'openxlsx')
sapply(load.lib,require,character=TRUE)
EOF




#========================================================================================#
#  ======== PREPARE R LIBRARIES =========================================================#
#========================================================================================#
R
source("input.R")
.libPaths(R_libPaths)
update.packages(ask=FALSE, checkBuilt=TRUE) # Update packages after installing new version of R

# Check if all required packages are installed
install.lib <-load.lib[!load.lib %in% rownames(installed.packages())]

# Install missing packages
for(lib in install.lib) {
    install.packages(lib,dependencies=TRUE, lib=R_libPaths)
  }

# Load all packages
sapply(load.lib,require,character=TRUE)

# install older lavaan version
sessionInfo()
#require(devtools)
packageVersion("lavaan") # new version 0.6.8
install_version("lavaan", version = "0.6.7", repos = "http://cran.us.r-project.org")
q(save="no")



#========================================================================================#
#  ======== PREPARE DATA FOR GWAS =======================================================#
#========================================================================================#

# Download LD scores and reference panel data
cd $HOME/data/processed
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -xvjf eur_w_ld_chr.tar.bz2
rm eur_w_ld_chr.tar.bz2
# go to https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v/folder/48313001768
# Download the following files and upload them into $HOME/data/processed
# Unizip folder
gunzip reference.1000G.maf.0.005.txt.gz


# Run Script to munge files
cat > $HOME/analysis/mungeSumstats.sh <<EOF
#!/bin/bash -l

#$ -l h_rt=8:30:0
#$ -l mem=1G
#$ -N mungeSumstats.sh
#$ -wd $HOME/output/munge
#$ -pe smp 12
#$ -l tmpfs=15G

# Load packages
module unload -f compilers mpi gcc-libs
module load beta-modules
module load r/r-4.0.2_bc-3.11

R --no-save < $HOME/analysis/mungeSumstats.R --args $HOME > $HOME/output/munge/mungeSumstats.log
EOF
qsub $HOME/analysis/mungeSumstats.sh
tail $HOME/output/mungeSumstats.log


#========================================================================================#
#  ======== RUN GWAS ====================================================================#
#========================================================================================#


GenSEMrun () {

echo "use $nCores cores"
echo "create script for GenSEM when testing model: $model"

cd $HOME/output/batchlogs
for input in $(sort -u $HOME/output/gwa_input/gwa_input_${model})
do
echo Start $input
echo "save output file as $output"


cat > $HOME/output/batchlogs/${input}_${output}.R <<EOF
#!/usr/bin/env Rscript

HOME="${HOME}"
# Load libraries
source(paste0(HOME, "/analysis/input.R"))
libfolder=paste0(HOME,"/programs/R")
.libPaths(libfolder)
sapply(load.lib,require,character=TRUE)
library(GenomicSEM)


# Function to write log file

writeOut=function(line){
    write(line,
    file=paste0(HOME, "/output/batchlogs/${input}_${output}", "_GenSEMlog.txt"),
    append=TRUE)
    }

writeOut("load ldsc regression output")
print(paste0("Filename for rds file: ",HOME,"/output/rds/LDSCoutput_${model}_SUD.rds"))
writeOut("Check if correct file is loaded")
LDSCoutput_SUD=readRDS(paste0(HOME,"/output/rds/LDSCoutput_${model}_SUD.rds"))


dir(paste0(HOME, "/output/gwa_input/"))
writeOut("Select subset of SNPs")
inputPath=paste0(HOME, "/output/gwa_input/$input")
Model_D_SS_trunc=read.table(inputPath, header = TRUE)

writeOut(paste0("Number of SNPs included for this batch job: ", length(Model_D_SS_trunc[["SNP"]])))


GenSEMpar=function(SNPdata){
.libPaths(libfolder)

writeOut(paste0("Number of SNPs included: ", length(SNPdata[["SNP"]])))
writeOut("Define SEM model")

source(paste0(HOME, "/output/batchlogs/model_${model}_${output}.R"))
print(model)

start_time = Sys.time()
writeOut(paste0("Start time user GWAS: ",  start_time))


con <- file("$input.${output}.console.userGWA.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

# Run the model
results = userGWAS(covstruc=LDSCoutput_SUD,
                   SNPs=SNPdata,
                   estimation = ${estimation}, #"DWLS"
                   model = model,
                   printwarn = TRUE,
                   toler = 1e-200,
                   parallel = FALSE,
                   SNPSE = 0.0005,
                   sub=${extract},
                   GC="none")
print(str(results))


# Restore output to console
sink()
sink(type="message")

results=as.data.frame(results)

writeOut("I'm done with the commonFactor!")
writeOut(paste0("End time user commonFactor: ",  Sys.time()))

writeOut("Finished analysis for subset SNPs - return output")

return(results)

}

library(parallel)

writeOut("Setup parallel processing")
int <- $nCores
cl <- makeCluster(int)

print("Display info about each process in the cluster")
print(clusterCall(cl, function() Sys.info()))

writeOut(paste0("Create ", int, " subsets"))
Model_D_SS_trunc_list=split(Model_D_SS_trunc, sample(1:int, nrow(Model_D_SS_trunc), replace=T))
writeOut(paste0("Length of list: ", length(Model_D_SS_trunc_list), " subsets"))

nSNPsIteration <- do.call("rbind", lapply(Model_D_SS_trunc_list , function(x) length(x[["SNP"]])))
writeOut("Number of SNPs per subset included in clusterApply")
writeOut( nSNPsIteration)

writeOut("Export functions")
writeOut("Export packages")
clusterExport(cl, varlist=c("GenSEMpar", "libfolder", "HOME", "writeOut", "LDSCoutput_SUD"))
clusterEvalQ(cl, libfolder )
clusterEvalQ(cl, .libPaths(libfolder) )
clusterEvalQ(cl, library(GenomicSEM) )
clusterEvalQ(cl, c("GenSEMpar", "libfolder", "HOME", "writeOut", "LDSCoutput_SUD", library(GenomicSEM) ))

writeOut("Run genomicSEM for $input")
resultList=parLapply(cl, Model_D_SS_trunc_list,function(x) GenSEMpar(SNPdata=x))

writeOut("Finished genomicSEM for $input")

writeOut("bind output together")
results <- do.call("rbind", resultList)

writeOut("Save output on the cluster")
write.csv(results, file= paste0(HOME, "/output/csv/${input}_${output}", ".csv"))
writeOut("All output is saved on the cluster")

EOF


cat > $HOME/output/batchlogs/${input}_${output}.sh <<EOF
#!/bin/bash -l

#$ -l h_rt=$time:00:0
#$ -l mem=0.5G
#$ -N ${input}_${output}.sh
#$ -wd $HOME/output/batchlogs
#$ -pe smp $nCores
#$ -l tmpfs=15G

# Load packages
module unload -f compilers mpi gcc-libs
module load beta-modules
module load r/r-4.0.2_bc-3.11

export R_LIBS=$HOME/programs/R:$R_LIBS

R --no-save < $HOME/output/batchlogs/${input}_${output}.R > $HOME/output/batchlogs/${input}_${output}.batch.log

EOF

echo "Submit ${input}_${output} as batch job"

qsub $HOME/output/batchlogs/${input}_${output}.sh

done

}

# ======= Run function =======
# Common factor model
model="cancigalc"
output="constr_corRes"
estimation='"DWLS"'
extract='c("F1~SNP")'
cat > $HOME/output/batchlogs/model_${model}_${output}.R <<EOF
model <- 'F1 =~ 1*CigaretteDependency + x*CigaretteDependency + y*DrinksPerWeek + y*CigarettesPerDay + x*AlcoholDependency + y*CannabisUseFrequency + x*CannabisUseDisorder
      DrinksPerWeek ~~ a*DrinksPerWeek
      CigarettesPerDay ~~ c*CigarettesPerDay
      CigaretteDependency ~~ d*CigaretteDependency
      AlcoholDependency ~~ e*AlcoholDependency
      CannabisUseFrequency ~~ f*CannabisUseFrequency
      CannabisUseDisorder ~~ g*CannabisUseDisorder
      CannabisUseFrequency ~~ CannabisUseDisorder
      CigarettesPerDay ~~ CigaretteDependency
      AlcoholDependency ~~ DrinksPerWeek
      a > .0001
      c > .0001
      d > .0001
      e > .0001
      f > .0001
      g > .0001
      F1 ~ SNP
      DrinksPerWeek ~ 0*SNP
      CigarettesPerDay ~ 0*SNP
      CigaretteDependency ~ 0*SNP
      AlcoholDependency ~ 0*SNP
      CannabisUseFrequency ~ 0*SNP
      CannabisUseDisorder ~ 0*SNP'
EOF
nCores=16 # define number of cors
time=47 # running time in hours
GenSEMrun "nCores" "model" "output" "estimation" "extract" "time"


# Common factor model (heterogeneity test)
model="cancigalc"
output="constr_corRes_HT"
estimation='"DWLS"'
extract='c("F1~~F1")'
cat > $HOME/output/batchlogs/model_${model}_${output}.R <<EOF
model <- 'F1 =~ 1*CigaretteDependency + x*CigaretteDependency + y*DrinksPerWeek + y*CigarettesPerDay + x*AlcoholDependency + y*CannabisUseFrequency + x*CannabisUseDisorder
      DrinksPerWeek ~~ a*DrinksPerWeek
      CigarettesPerDay ~~ c*CigarettesPerDay
      CigaretteDependency ~~ d*CigaretteDependency
      AlcoholDependency ~~ e*AlcoholDependency
      CannabisUseFrequency ~~ f*CannabisUseFrequency
      CannabisUseDisorder ~~ g*CannabisUseDisorder
      CannabisUseFrequency ~~ CannabisUseDisorder
      CigarettesPerDay ~~ CigaretteDependency
      AlcoholDependency ~~ DrinksPerWeek
      a > .0001
      c > .0001
      d > .0001
      e > .0001
      f > .0001
      g > .0001
      CigaretteDependency + DrinksPerWeek + CigarettesPerDay + AlcoholDependency + CannabisUseFrequency + CannabisUseDisorder ~ SNP'
EOF
nCores=16 # define number of cors
time=47 # running time in hours
GenSEMrun "nCores" "model" "output" "estimation" "extract" "time"


#========================================================================================#
#  ======== PROCESS GWAS OUTPUT =========================================================#
#========================================================================================#

modelName="cancigalc_constr_corRes"
cat > $HOME/analysis/processingMultiGWA.sh <<EOF
#!/bin/bash -l

#$ -l h_rt=3:30:0
#$ -l mem=1G
#$ -N processingMultiGWA.sh
#$ -wd $HOME/output/processingMultiGWA
#$ -pe smp 10
#$ -l tmpfs=15G

# Load packages
module unload -f compilers mpi gcc-libs
module load beta-modules
module load r/r-4.0.2_bc-3.11

R --no-save < $HOME/analysis/processingMultiGWA.R --args $HOME $modelName > $HOME/output/processingMultiGWA/processingMultiGWA.log
EOF
qsub $HOME/analysis/processingMultiGWA.sh
tail $HOME/output/processingMultiGWA/processingMultiGWA.log



#========================================================================================#
#  ======== LD SCORE REGRESSION WITH OTHER TRAITS =======================================#
#========================================================================================#

modelName="cancigalc_constr_corRes"
cat > $HOME/analysis/ldscoreregressionAnalysis.sh <<EOF
#!/bin/bash -l

#$ -l h_rt=8:30:0
#$ -l mem=1G
#$ -N ldscoreregressionAnalysis.sh
#$ -wd $HOME/output/ldsc
#$ -pe smp 14
#$ -l tmpfs=15G

# Load packages
module unload -f compilers mpi gcc-libs
module load beta-modules
module load r/r-4.0.2_bc-3.11

R --no-save < $HOME/analysis/ldscoreregressionAnalysis.R --args $HOME $modelName > $HOME/output/ldsc/ldscoreregressionAnalysis.log
EOF
qsub $HOME/analysis/ldscoreregressionAnalysis.sh
tail $HOME/output/ldsc/ldscoreregressionAnalysis.log


#========================================================================================#
#  ======== PATHWAY ANALYSIS: PASCAL ====================================================#
#========================================================================================#
# Download PASCAL
# wget http://www2.unil.ch/cbg/images/3/3d/PASCAL.zip -P $HOME/programs
# unzip $HOME/programs/PASCAL.zip

# Open and change the following
# cd xianyi-OpenBLAS-48f06dd/
# make
# make install PREFIX="../lib/openBLASlib"

# To
# make PREFIX="../lib/openBLASlib" TARGET=HASWELL
# make install PREFIX="../lib/openBLASlib" TARGET=HASWELL

# make file executable
# chmod +x $HOME/programs/PASCAL/installScript.sh
# cd $HOME/programs/PASCAL/
# ./installScript.sh

# download depict dataset
# https://drive.google.com/file/d/1IiIRrmHUHoQB7Ve9P1N660AIhbzTh9kw/view


# Prepare summary statistics files
cd $HOME/analysis
R
source("input.R")
libfolder=paste0(HOME,"/programs/R")
.libPaths(libfolder)

modelName="cancigalc_constr_corRes"
GWASsumStats=readRDS(paste0(HOME,"/output/rds/" ,modelName, "_Sumstats.rds"))

# Remove heterogenous SNPs from common liability GWAS
sumStatsPascal=GWASsumStats
sumStatsPascal_clean=subset(GWASsumStats[[1]], Q_chisq_pval >= 5e-08)

# Add filtered GWAS to sumstats list
nlist=length(sumStatsPascal)+1
sumStatsPascal[[nlist]]=sumStatsPascal_clean
names(sumStatsPascal)[nlist]=paste0(names(GWASsumStats)[1], "_filtered")

# Formatting of summary statistic files
funcSelectCol=function(df){
  df_sub=subset(df, select=c(SNP, P))
  return(df_sub)
}
sumStatsPascalSub=lapply(sumStatsPascal, function(x) funcSelectCol(x))
names(sumStatsPascalSub)=names(sumStatsPascal)

# Save output
setwd(paste0(HOME, "/data/pascal"))
sapply(names(sumStatsPascalSub), function(x)
  write.table(sumStatsPascalSub[[x]], file=paste(x, "pascal", sep="."),
              col.names=F, row.names=F, quote=F, sep="\t") )

paste0(names(sumStatsPascalSub), ".pascal")

# Create list for loop
write.table(paste0(names(sumStatsPascalSub), ".pascal"),
            file=paste0(HOME, "/data/pascal/pascalInput"),col.names=F, row.names=F, quote=F, sep=" " )

q(save="no")

# ===== NOTES
# Two options are available:
#     1) the max of chi-squared statistics => based on the most significant SNP
#     2) sum of chi-squared statistics => based on the average association signal across the region
# --genescoring=sum // # --genescoring=max

# Submit in loop
for input in $(cat $HOME/data/pascal/pascalInput)
do
echo Start $input
#input="commonLiability_filtered.pascal"
cat > $HOME/programs/PASCAL/output/$input.sh <<EOF
#!/bin/bash -l

#$ -l h_rt=30:00:0
#$ -l mem=5G
#$ -pe smp 13
#$ -N $input.sh
#$ -wd $HOME/programs/PASCAL

# Load packages
module unload compilers mpi
module load gcc-libs/4.9.2
module load compilers/gnu/4.9.2
module load java/1.8.0_92

cd $HOME/programs/PASCAL # important: working directory has to be the pascal folder (otherwise it won't run)

# REACTOME
$HOME/programs/PASCAL/Pascal \
--pval=$HOME/data/pascal/$input \
--genesetfile=$HOME/programs/PASCAL/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt \
--runpathway=on \
--genescoring=sum \
--outsuffix=.msigdb.$input

# DEPICT
$HOME/programs/PASCAL/Pascal \
--pval=$HOME/data/pascal/$input \
--genesetfile=$HOME/programs/PASCAL/resources/genesets/depict/depict_discretized_cutoff3.2.gmt \
--runpathway=on \
--genescoring=sum \
--outsuffix=.depict.$input

EOF

echo Submit btach job for $input
qsub $HOME/programs/PASCAL/output/$input.sh

done


# =========== Process PASCAL results =================
cd $HOME/analysis
R
# Load libraries
source("input.R")

# Go to PASCAL directory
setwd(paste0(HOME, "/programs/PASCAL/output"))

# Get file names (REACTOME)
filesReac = dir(recursive=F, full.names=F, pattern="*msigBIOCARTA_KEGG_REACTOME*")
filesReacName = unlist(lapply(filesReac, function(x) strsplit(x, ".PathwaySet--msigBIOCARTA_KEGG_REACTOME--")[[1]][1]))
# Read in REACTOME
listPascalReact= lapply(filesReac, function(x) fread(paste0(HOME, "/programs/PASCAL/output/", x),header=T,data.table=T))
names(listPascalReact)=filesReacName

# Get file names (DEPICT)
filesDepict = dir(recursive=F, full.names=F, pattern="*depict_discretized_cutoff*")
filesDepictName = unlist(lapply(filesDepict, function(x) strsplit(x, ".PathwaySet--depict_discretized_cutoff3.2")[[1]][1]))
# Read in DEPICT
listPascalDepict= lapply(filesDepict, function(x) fread(paste0(HOME, "/programs/PASCAL/output/", x),header=T,data.table=T))
names(listPascalDepict)=filesDepictName

# Export output to cluster
saveRDS(listPascalReact, paste0(HOME, "/output/rds/PascalReact.rds"))
saveRDS(listPascalDepict, paste0(HOME, "/output/rds/PascalDepict.rds"))

# Export to dropbox
drop_auth(rdstoken = paste0(HOME, "/token.rds"))  # authentication for dropbox
drop_upload(paste0(HOME,"/output/rds/PascalReact.rds"), path = paste0(LOCAL, "/output/rds"), mode = "overwrite")
drop_upload(paste0(HOME,"/output/rds/PascalDepict.rds"), path = paste0(LOCAL, "/output/rds"), mode = "overwrite")

q(save="no")




cd $HOME/analysis
R
# Load libraries
source("input.R")
# Upload scripts to dropbos
drop_auth(rdstoken = paste0(HOME, "/token.rds"))  # authentication for dropbox
drop_upload(paste0(HOME,"/analysis/mungeSumstats.R"), path = paste0(LOCAL, "/analysis"), mode = "overwrite")
drop_upload(paste0(HOME,"/analysis/processingMultiGWA.R"), path = paste0(LOCAL, "/analysis"), mode = "overwrite")
drop_upload(paste0(HOME,"/analysis/ldscoreregressionAnalysis.R"), path = paste0(LOCAL, "/analysis"), mode = "overwrite")


