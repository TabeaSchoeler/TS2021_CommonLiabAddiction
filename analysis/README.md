Analysis
================

</br></br></br></br>

## Table of Contents

1.  [Install R packages](#packages)
2.  [Download summary statistic files](#model)
3.  [Pre-processing of summary statistic files](#pre)
4.  [Munging of statistic files](#munge)
5.  [Specification of the structural model](#strucModel)
6.  [Run the common liability GWA analysis](#gwaRun)
7.  [Processing the common liability GWA output](#gwaOutput)

</br></br></br></br>

# [Install R packages](#packages)

``` r
# Load and install libraries
.libPaths( c( .libPaths(), "~/Dropbox/progs/R/library") )
load.lib=c('pheatmap', 'ggthemes', 'ggpubr', 'Qtlizer', 'devtools', 'GenomicSEM', 'data.table', 'Matrix',
           'sem', 'Matrix', 'stats', 'semTools', 'ggcorrplot', 'grex', 'openxlsx', 'rJava', 'qqman', 'gmodels', 'ggplot2',
           'ggcorrplot', 'tidyverse', 'reshape2', 'pdftools', 'plyr', 'phenoscanner','ggpubr',
           'biomaRt', 'viridis', 'rdrop2', 'lavaan', 'gprofiler2', 'knitr',
           'gridExtra', 'CMplot', 'configr', 'purrr', 'TwoSampleMR', 'GOexpress', 'ieugwasr', 'xlsx')


install.lib<-load.lib[!load.lib %in% installed.packages()]

# Install missing packages
for(lib in install.lib) install.packages(lib, dependencies=TRUE)
# Load all packages
sapply(load.lib,require,character=TRUE)
```

</br></br>

# [Download summary statistic files](#model)

| Link to publication                                                                             | Link to summary statistic file                                                                                                                                                       |
|:------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [Cannabis use (frequency) (Publication)](http://www.nealelab.is/uk-biobank)                     | [Cannabis use (frequency) (Sumstats)](http://www.nealelab.is/uk-biobank)                                                                                                             |
| [Cannabis use (dependence) (Publication)](https://doi.org/10.1016/S2215-0366(20)30339-4)        | [Cannabis use (dependence) (Sumstats)](https://www.med.unc.edu/pgc/download-results/sud)                                                                                             |
| [Cocaine use (dependence) (Publication)](https://doi.org/10.1016/j.pnpbp.2019.109667)           | [Cocaine use (dependence) (Sumstats)](personal%20correspondence)                                                                                                                     |
| [Risk tolerance (Publication)](https://doi.org/10.1038/s41588-018-0309-3)                       | [Risk tolerance (Sumstats)](https://www.thessgac.org/data)                                                                                                                           |
| [Car speeding propensity (Publication)](https://doi.org/10.1038/s41588-018-0309-3)              | [Car speeding propensity (Sumstats)](https://www.thessgac.org/data)                                                                                                                  |
| [Irritability (Publication)](http://www.nealelab.is/uk-biobank)                                 | [Irritability (Sumstats)](http://www.nealelab.is/uk-biobank)                                                                                                                         |
| [Miserableness (Publication)](http://www.nealelab.is/uk-biobank)                                | [Miserableness (Sumstats)](http://www.nealelab.is/uk-biobank)                                                                                                                        |
| [Mood swings (Publication)](http://www.nealelab.is/uk-biobank)                                  | [Mood swings (Sumstats)](http://www.nealelab.is/uk-biobank)                                                                                                                          |
| [Nervous feelings (Publication)](http://www.nealelab.is/uk-biobank)                             | [Nervous feelings (Sumstats)](http://www.nealelab.is/uk-biobank)                                                                                                                     |
| [Coffee intake (Publication)](http://www.nealelab.is/uk-biobank)                                | [Coffee intake (Sumstats)](http://www.nealelab.is/uk-biobank)                                                                                                                        |
| [ADHD (Publication)](https://doi.org/10.1038/s41588-018-0269-7)                                 | [ADHD (Sumstats)](https://www.med.unc.edu/pgc/download-results/)                                                                                                                     |
| [Anxiety (Publication)](https://doi.org/10.1038/mp.2015.197)                                    | [Anxiety (Sumstats)](https://www.med.unc.edu/pgc/download-results/)                                                                                                                  |
| [Neuroticism (Publication)](https://doi.org/10.1038/s41588-018-0151-7)                          | [Neuroticism (Sumstats)](https://ctg.cncr.nl/software/summary_statistics)                                                                                                            |
| [Depression (diagnosis) (Publication)](https://doi.org/10.1038/s41588-018-0090-3)               | [Depression (diagnosis) (Sumstats)](https://www.med.unc.edu/pgc/download-results/)                                                                                                   |
| [Intelligence (Publication)](https://doi.org/10.1038/s41588-018-0152-6)                         | [Intelligence (Sumstats)](https://ctg.cncr.nl/software/summary_statistics)                                                                                                           |
| [Insomnia (Publication)](https://doi.org/10.1038/ng.3888)                                       | [Insomnia (Sumstats)](http://ctg.cncr.nl/software/summary_statistics)                                                                                                                |
| [Sleep duration (Publication)](https://doi.org/10.1038/ng.3888)                                 | [Sleep duration (Sumstats)](http://ctg.cncr.nl/software/summary_statistics)                                                                                                          |
| [Schizophrenia (Publication)](https://doi.org/10.1101/2020.09.12.20192922)                      | [Schizophrenia (Sumstats)](https://www.med.unc.edu/pgc/download-results/)                                                                                                            |
| [BMI (Publication)](https://doi.org/10.1093/hmg/ddy271)                                         | [BMI (Sumstats)](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GIANT_Consortium_2016_Exome_Array_Data_is_Available_Here_for_Download) |
| [Educational attainment (Publication)](https://doi.org/10.1038/s41588-018-0147-3)               | [Educational attainment (Sumstats)](https://www.thessgac.org/data)                                                                                                                   |
| [Autism (diagnosis) (Publication)](https://doi.org/10.1038/s41588-019-0344-8)                   | [Autism (diagnosis) (Sumstats)](https://ipsych.dk/en/research/downloads/)                                                                                                            |
| [Depressive symptoms (Publication)](https://doi.org/10.1038/ng.3552)                            | [Depressive symptoms (Sumstats)](https://www.thessgac.org/data)                                                                                                                      |
| [Extraversion (IRT) (Publication)](https://doi.org/10.1007/s10519-014-9654-x)                   | [Extraversion (IRT) (Sumstats)](https://tweelingenregister.vu.nl/gpc)                                                                                                                |
| [Openness (Publication)](https://doi.org/10.1038/mp.2010.128)                                   | [Openness (Sumstats)](https://tweelingenregister.vu.nl/gpc)                                                                                                                          |
| [Sexual partners (number) (Publication)](https://doi.org/10.1038/s41588-018-0309-3)             | [Sexual partners (number) (Sumstats)](https://www.thessgac.org/data)                                                                                                                 |
| [Anorexia (diagnosis) (Publication)](https://doi.org/10.1038/s41588-019-0439-2)                 | [Anorexia (diagnosis) (Sumstats)](https://www.med.unc.edu/pgc/download-results/)                                                                                                     |
| [Height (Publication)](https://doi.org/10.1093/hmg/ddy271)                                      | [Height (Sumstats)](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files)                                                                    |
| [Birth weight (Publication)](http://www.nealelab.is/uk-biobank)                                 | [Birth weight (Sumstats)](http://www.nealelab.is/uk-biobank)                                                                                                                         |
| [OCD (diagnosis) (Publication)](https://doi.org/10.1038/mp.2017.154)                            | [OCD (diagnosis) (Sumstats)](https://www.med.unc.edu/pgc/download-results/)                                                                                                          |
| [Loneliness (Publication)](doi.org/10.17863/CAM.23511)                                          | [Loneliness (Sumstats)](https://www.repository.cam.ac.uk/handle/1810/277812)                                                                                                         |
| [Social isolation (Publication)](doi.org/10.17863/CAM.23511)                                    | [Social isolation (Sumstats)](https://www.repository.cam.ac.uk/handle/1810/277812)                                                                                                   |
| [Income (Publication)](doi.org/10.1038/s41467-019-13585-5)                                      | [Income (Sumstats)](ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/HillWD_31844048_GCST009524)                                                                            |
| [Cigarette use (age of onset) (Publication)](https://www.nature.com/articles/s41588-018-0307-5) | [Cigarette use (age of onset) (Sumstats)](https://conservancy.umn.edu/handle/11299/201564)                                                                                           |
| [Opioid use (dependence) (Publication)](NA)                                                     | [Opioid use (dependence) (Sumstats)](https://www.med.unc.edu/pgc/download-results/sud)                                                                                               |
| [Cigarette use (dependence) (Publication)](https://doi.org/10.1038/s41588-018-0248-z)           | [Cigarette use (dependence) (Sumstats)](https://atlas.ctglab.nl/traitDB/3689)                                                                                                        |
| [Alcohol use (dependence) (Publication)](https://doi.org/10.1038/s41593-018-0275-1)             | [Alcohol use (dependence) (Sumstats)](https://www.med.unc.edu/pgc/download-results/sud)                                                                                              |
| [Cigarette use (frequency) (Publication)](https://www.nature.com/articles/s41588-018-0307-5)    | [Cigarette use (frequency) (Sumstats)](https://conservancy.umn.edu/handle/1129+O38+N38+M3+C38)                                                                                       |
| [Alcohol use (frequency) (Publication)](https://www.nature.com/articles/s41588-018-0307-5)      | [Alcohol use (frequency) (Sumstats)](https://conservancy.umn.edu/handle/11299/201564)                                                                                                |
| [Bipolar (diagnosis) (Publication)](https://doi.org/10.1101/2020.09.17.20187054)                | [Bipolar (diagnosis) (Sumstats)](https://figshare.com/s/53af8c444bb949be273c)                                                                                                        |
| [Cortical thickness (Publication)](https://doi.org/10.1126/science.aay6690)                     | [Cortical thickness (Sumstats)](http://enigma.ini.usc.edu/protocols/genetics-protocols/)                                                                                             |
| [Cortical surface area (Publication)](https://doi.org/10.1126/science.aay6690)                  | [Cortical surface area (Sumstats)](http://enigma.ini.usc.edu/protocols/genetics-protocols/)                                                                                          |

</br></br>

# [Pre-processing of summary statistic files](#pre)

-   Below is the pipeline that runs all the pre-processing steps of the
    summary statistic files and then performs the LD score regression
    analysis
-   The code below is copied from the shell script
    [‘CommonLiabAddiction.sh’](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/CommonLiabAddiction.sh),
    which relies on the R script
    [‘CommonLiabAddiction.R’](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/CommonLiabAddiction.R)

</br>

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

</br>

-   the script
    [‘CommonLiabAddiction.R’](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/CommonLiabAddiction.R)
    will read in all the summary statistic files from
    [GenSEM\_GWAS.xlsx](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/)
    that are selected to be included in the multivariate genome-wide
    association analysis. Selected are those that those with
    `directory == "gensem"`

-   all summary statistic files have to to be stored in a parent
    directory (e.g. `gwasDir="path/to/sumstats"`), with all
    sub-directories being specified in the column “fileName” in the
    excel sheet
    [GenSEM\_GWAS.xlsx](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/)

-   for each summary statistic file, a couple of details have to be
    included in the excel sheet, including the column names (e.g. label
    used for chromosome, effect allele, effect estimate etc.), the
    outcome format (binary vs. continuous) and the sample size

-   for binary summary statistic files, two additional estimates have to
    be included, namely the sample prevalence (number of cases / all
    participants) and the population prevalence (i.e. the assumed
    prevalence of affected individuals in the general population).
    Importantly, the columns indexing the effect estimates have to
    reflect log odds, and the column including the standard errors have
    to reflect the standard errors of the log odds. If studies report
    Odds ratios, these have to be converted to log odds before running
    the analyses

</br>

``` r
# Create path name
gwasSumStastAll$DirPath=paste0(gwasDir,gwasSumStastAll$fileName)

# Import data from dropbox 
drop_auth(rdstoken = paste0(HOME, "/token.rds"))  # authentication for dropbox
drop_acc()
drop_download(
  local_path = paste0(HOME,"/analysis/GenSEM_GWAS.xlsx"),
  path = paste0(LOCAL,"/analysis/GenSEM_GWAS.xlsx"), overwrite=TRUE)
gwasSumStastAll=read.xlsx(paste0(HOME,"/analysis/GenSEM_GWAS.xlsx"), na.strings="NA")

gwasSumStastAll$DirPathgz=paste0(HOME,"/data/processed/", gwasSumStastAll$label, ".sumstats.gz") # Include original GWAS as well as munged 

# GWAS files
# Select only substance use traits
gwasSumStast = subset(gwasSumStastAll, directory == "gensem" )
```

</br>

-   next, the selected summary statistic files are read into R and
    prepared for munging
-   from all summary statistic files, the columns indexing the SNP, A1
    (affect allele), A2 (alternative allele), BETA (the effect
    estimate), SE (the standard error of the effect estimate), P (the
    corresponding p-value) and N (the sample size) are selected and
    relabeled
-   after removing duplicate SNPs, the summary statistic files are
    stored

``` r
# Prepare clean datasets for munging
for ( i in 1:length(gwasSumStast$DirPath) ) { 
  print(gwasSumStast$DirPath[i])
  gwa_in=fread(gwasSumStast$DirPath[i])
  
  gwa_select=subset(gwa_in , select=c(gwasSumStast$SNP[i], 
                                      gwasSumStast$A1[i],
                                      gwasSumStast$A2[i],
                                      gwasSumStast$effect[i],
                                      gwasSumStast$se[i],
                                      gwasSumStast$P[i],
                                      gwasSumStast$N_snp[i]))
  
  colnames(gwa_select)=c("SNP", "A1", "A2", "BETA", "SE", "P", "N")
  print("Make numeric")
  gwa_select$BETA=as.numeric(gwa_select$BETA)  
  gwa_select$SE=as.numeric(gwa_select$SE)  
  gwa_select$P=as.numeric(gwa_select$P)  
  
  gwa_unique=gwa_select[!duplicated(gwa_select$SNP)]
  print(paste0("Remove ", NROW(gwa_select)-NROW(gwa_unique), " duplicate SNPs"))
  # Save output
  print("Save output")
  write.table(gwa_unique, 
              file=paste0(HOME, "/data/processed/", gwasSumStast$label[i], "_clean"), 
              sep="\t", 
              row.names = FALSE, 
              col.names = TRUE, 
              quote=F) 
  head(gwa_unique)
  gwa_unique=NULL
}
```

</br></br>

# [Munging of statistic files](#munge)

-   the files are now munged using the `munge()` function implemented in
    `GenomicSEM`
-   the documentation of this function can be found
    [here](https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects)

``` r
for ( i in 1:length(gwasSumStast$DirPath) ) { 
  setwd(paste0(HOME,"/data/processed"))
  print("Remove munged file if already in folder")
  system(paste0("rm ", paste0(HOME, "/data/processed/", gwasSumStast$label[i], ".sumstats.gz") ))
  
  munge(files = paste0(HOME, "/data/processed/", gwasSumStast$label[i], "_clean"), 
        hm3 = paste0(HOME,"/data/processed/w_hm3.snplist") ,
        trait.names=gwasSumStast$label[i],
        info.filter = 0.9, 
        maf.filter = 0.01) 
}
```

</br>

-   in the next steps, multivariate LD score regression analysis is
    performed
-   Genomic control was applied to all summary statistics showing
    evidence of uncontrolled confounding (LD score intercept &gt; 1), by
    multiplying standard errors by the LD score intercept

``` r
# ======================================= LD SCORE REGRESSION
LDSCoutput_SUD <- ldsc(gwasSumStast$DirPathgz, 
                       samplePrevalence, 
                       gwasSumStast$population.prev, 
                       ld=paste0(HOME,"/data/processed/eur_w_ld_chr/"), ## set the directory where the LD scores used the used in the analysis are located
                       wld=paste0(HOME,"/data/processed/eur_w_ld_chr/"), 
                       gwasSumStast$label, 
                       stand= TRUE)
# Format decimals
correlationSUD=round(LDSCoutput_SUD$S_Stand, 3)

# Replace cor=1 with heritability estimates
h2out=data.frame(trait=colnames(correlationSUD),
                 est= rep(NA, length(colnames(correlationSUD))))

for (i in 1:length(colnames(correlationSUD))) {
  h2out$est[i]=round(LDSCoutput_SUD$S[i,i],3)
}

h2out
rownames(correlationSUD)=colnames(correlationSUD)
correlationSUD
#colnames(LDSCoutput_SUD$S_Stand)=recodeName(colnames(LDSCoutput_SUD$S_Stand))
#rownames(LDSCoutput_SUD$S_Stand)=colnames(LDSCoutput_SUD$S_Stand)
LDSCoutput_SUD

# ==== Check and replace GWAs where genomic control has to be applied 
str(LDSCoutput_SUD)
intercept=LDSCoutput_SUD$I
colnames(intercept)=colnames(LDSCoutput_SUD$S)
rownames(intercept)=colnames(LDSCoutput_SUD$S)
checkGC=data.frame(pheno = colnames(intercept),
                   intercept = diag(intercept)) # get the LDSC intercept per sumstats

checkGC=subset(checkGC, intercept >=1 ) # select problematic phenotypes


for ( i in 1:length(checkGC$pheno) ) {
  print(paste0("Read in data for ", checkGC$pheno[i]))
  phenoDirmunge=paste0(HOME,"/data/processed/", checkGC$pheno[i], ".sumstats.gz") 
  phenoDir=paste0(HOME,"/data/processed/", checkGC$pheno[i], "_clean") 
  data <- fread(phenoDir,header=T,data.table=F)
  
  N <- fread(phenoDirmunge,header=T,data.table=F)$N[1]
  print(paste0("Multiply SE by LDSC intercept of ", checkGC$intercept[i]))
  data$SE = data$SE * sqrt(checkGC$intercept[i])
  data$Zscore=as.numeric( data$BETA)/as.numeric( data$SE) # estimate adjusted p-values
  data$P=2*pnorm(-abs(data$Zscore))
  data$Zscore=NULL
  head(data)
  
  print("Save file on cluster")
  write.table(data, 
              file=paste0(HOME, "/data/processed/", checkGC$pheno[i], "_clean"), 
              sep="\t", 
              row.names = FALSE, 
              col.names = TRUE, 
              quote=F) 
  print("Remove non-GC controlled munged file and munge again")
  setwd(paste0(HOME,"/data/processed"))
  system(paste0("rm ", phenoDirmunge))
  
  munge(files = paste0(HOME, "/data/processed/", checkGC$pheno[i], "_clean"), 
        hm3 = paste0(HOME,"/data/processed//w_hm3.snplist") ,
        trait.names=checkGC$pheno[i],
        N = N,
        info.filter = 0.9, 
        maf.filter = 0.01) 
}

# ===== Re-run LDSC regression including GCed sumstats

LDSCoutput_SUD <- ldsc(gwasSumStast$DirPathgz, 
                       samplePrevalence, 
                       gwasSumStast$population.prev, 
                       ld=paste0(HOME,"/data/processed/eur_w_ld_chr/"), ## set the directory where the LD scores used the used in the analysis are located
                       wld=paste0(HOME,"/data/processed/eur_w_ld_chr/"), 
                       gwasSumStast$label, 
                       stand= TRUE)

data.frame(pheno = colnames(LDSCoutput_SUD$S_Stand),
           intercept = diag(LDSCoutput_SUD$I)) 

# Export and upload onto Dropbox
modelName="cancigalc"
saveRDS(LDSCoutput_SUD, paste0(HOME,"/output/rds/LDSCoutput_", modelName, "_SUD.rds"))
# Export to dropbox
drop_auth(rdstoken = paste0(HOME, "/token.rds"))  # authentication for dropbox
drop_acc()
drop_upload(paste0(HOME,"/output/rds/LDSCoutput_", modelName, "_SUD.rds"), path = paste0(LOCAL, "/output/rds")) 
```

</br></br>

# [Specification of the structural model](#strucModel)

-   the common liability model is specified using lavaan syntax
-   Equality constrains were imposed on paths belonging to the same
    pattern of substance use, i.e., equal weights across measures of
    dependence, and equal weights across measures of frequency of use
-   Correlated residuals were included to allow for within-substance
    class associations

``` r
# Constrained and correalted residuals
model <- 'F1 =~ 1*CigaretteDependency + x*CigaretteDependency + y*DrinksPerWeek + y*CigarettesPerDay + x*AlcoholDependency + y*CannabisUseFrequency + x*CannabisUseDisorder
      DrinksPerWeek ~~ a*DrinksPerWeek
      CigarettesPerDay ~~ c*CigarettesPerDay
      CigaretteDependency ~~ d*CigaretteDependency
      AlcoholDependency ~~ e*AlcoholDependency
      CannabisUseFrequency ~~ f*CannabisUseFrequency
      CannabisUseDisorder ~~ g*CannabisUseDisorder
      AlcoholDependency ~~ DrinksPerWeek
      CigaretteDependency ~~ CigarettesPerDay
      CannabisUseDisorder ~~ CannabisUseFrequency
      a > .0001
      c > .0001
      d > .0001
      e > .0001
      f > .0001
      g > .0001'

# Build the model
CommonFac_model=usermodel(covstruc=LDSCoutput_SUD, 
                          estimation = "DWLS", 
                          model = model,  
                          CFIcalc = TRUE, 
                          std.lv = FALSE, 
                          imp_cov = FALSE)
#  SRMR: A value less than .08 is generally considered a good fit (Hu & Bentler, 1999). 
# Get model fit indicies
CommonFac_model$modelfit
resultsLDSC=data.frame(pheno=paste0(CommonFac_model$results$lhs, CommonFac_model$results$op, CommonFac_model$results$rhs),
                       loading = round(CommonFac_model$results$STD_All, 2))

# Export and upload onto Dropbox
saveRDS(CommonFac_model, paste0(HOME,"/output/rds/CommonFac_model.rds"))
drop_upload(paste0(HOME,"/output/rds/CommonFac_model.rds"), path = paste0(LOCAL, "/output/rds"))                      
```

-   Prepare the summary statistics for the multivariate GWA analysis

``` r
se.logit=NULL
se.logit=ifelse(gwasSumStast$type=="con", F, NA)
se.logit=ifelse(gwasSumStast$type=="binary", T, se.logit)

# indicate whether trait is continuous or binary
OLS=ifelse(gwasSumStast$type=="con", T, NA) # OLS=T for continuous traits
OLS=ifelse(gwasSumStast$type=="binary", F, OLS) # OLS=F for binary traits

# Hail model for biobank data
lin_probHail=ifelse(gwasSumStast$linprob=="yes", T, F)
propHail=ifelse(gwasSumStast$linprob=="yes", samplePrevalence, NA)

setwd(paste0(HOME,"/output"))

SUD_sumstats <- sumstats(files=paste0(HOME, "/data/processed/", gwasSumStast$label, "_clean"), # should be the same as the name of the files used for the munge function in Step 1 
                         ref=paste0(HOME,"/data/processed/reference.1000G.maf.0.005.txt"), # The reference file used to calculate SNP variance across traits. We use 1000 genomes phase 3 
                         trait.names=gwasSumStast$label,
                         se.logit=se.logit, 
                         info.filter=.6, 
                         maf.filter=0.01,
                         linprob = lin_probHail,
                         prop = propHail,
                         OLS = OLS)
head(SUD_sumstats)

# ===== Save datasets with subsets of SNPs
# Split in subsets of SNPs
baseModel="cancigalc"

n_chunks=500000
n_rows <- nrow(SUD_sumstats)
r  <- rep(1:ceiling(n_rows/n_chunks),each=n_chunks)[1:n_rows]
SUD_sumstats_list <- split(SUD_sumstats,r)
n_splits=seq(1:length(SUD_sumstats_list)) # Number of subsets:75
names(SUD_sumstats_list)=paste0("GenSem_sub_par", n_splits, "_", baseModel)

# Save output on cluster
sapply(names(SUD_sumstats_list), function(x) 
  write.table(SUD_sumstats_list[[x]], file=paste0(HOME, "/output/gwa_input/", x), 
              col.names=T, row.names=F, quote=F, sep="\t") )

# Create list for loop
gwa_input_name=names(SUD_sumstats_list)
write.table(gwa_input_name, 
            file=paste0(HOME, "/output/gwa_input/gwa_input_", baseModel),
            col.names=F, 
            row.names=F, 
            quote=F, 
            sep=" " )
```

</br></br>

# [Run the common liability GWA analysis](#gwaRun)

-   after having executed the script
    `R --no-save < $HOME/analysis/mungeSumstats.R`, the summary
    statistic files and the LD score regression output have been
    prepared for the GWA analysis

-   the analysis is performed using the user specified GWA function
    `userGWAS()` as implemented in GenomicSEM
    (cf. [here](https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects)
    for a detailed description)

-   the script below is used to speed up the analysis, by splitting all
    included SNPs into chunks of 500,000 (cf. script aboved
    `n_chunks=500000`) and running the GWA using parallel processing

<!-- -->

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

-   once the above function is specified, we can use to run the GWA
-   the function can process different structural models and parameters
-   `estimation='"DWLS"'` indexes that diagonally weighted least squares
    is used for estimation (alternatively, maximum likelihood can be
    specified with `estimation='"ML"'`)
-   `extract='c("F1~SNP")'` is used to tell the program that only the
    SNP effects on the common liability should be retuned in the output.

### GWA analysis on the common liability

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

### GWA analysis to obtain the heterogeneity statistics

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

</br></br>

# [Processing the common liability GWA output](#gwaOutput)

-   once we have the GWA results, we can submit the Rscript
    [processingMultiGWA.R](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis)
    to process the results

-   this script will combine all the chunks of SNPs for which the
    effects on the common liability was estimates

-   it also includes the estimates from the heterogeneity test, which
    were obtained through a GWA of a model in which the SNP effects were
    not allowed to operate through the common liability

-   finally, it also include the results from the original summary
    statistic files used to derive the common liability

-   the output is a list of data including the common liability and the
    individual substance use phenotypes

<!-- -->

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

[create an anchor](#anchors-in-markdown)

[](#mendelian-randomization-analysis)

# Mendelian Randomization Analysis
