Analysis
================

</br></br></br></br>

# Install R packages

[](#install-r-packages)

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

# Download summary statistic files

[](#download-summary-statistic-files)

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

Summary statistic files

</br></br>

# Pre-processing of the summary statistic files

[](#pre-processing-of-the-summary-statistic-files)

-   Below is the pipeline that runs all the pre-processing steps of the
    summary statistic files and then performs the LD score regression
    analysis
-   The code below is copied from the shell script
    [CommonLiabAddiction.sh](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/CommonLiabAddiction.sh),
    which relies on the R script
    [CommonLiabAddiction.R](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/CommonLiabAddiction.R)

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
    [CommonLiabAddiction.R](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/CommonLiabAddiction.R)
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

# Munging of statistic files

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

</br></br></br>

# Estimate genetic correlations

[](#estimate-genetic-correlations)

-   in this step, multivariate LD score regression analysis is performed
-   Please note that genomic control was applied to all summary
    statistics showing evidence of uncontrolled confounding (LD score
    intercept &gt; 1), by multiplying standard errors by the LD score
    intercept

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

</br></br></br>

# Specify the structural model

[](#specify-the-structural-model)

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

# Run the multivariate genome-wide association study

[](#run-the-multivariate-genome-wide-association-study)

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

</br></br></br>

# [Processing the common liability GWA output](#ld-score-regression-analysis-including-other-traits)

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

</br></br></br>

# Clumping and gene mapping

[](#clumping-and-gene-mapping)

``` r
##################################################################
# === Read in results  ===========================================
##################################################################


# Function to Merge as rbind
mergeLongFormat=function(list, label, cleanLabel){
  dataframe = mapply(`[<-`, list, 'Label', value = label, SIMPLIFY = FALSE)
  dataframeOut=dplyr::bind_rows(dataframe)
  if (cleanLabel==TRUE) {
    dataframeOut$Label=ifelse(duplicated(dataframeOut$Label)==FALSE, as.character(dataframeOut$Label), " ")
    return(dataframeOut)
  } else {
    return(dataframeOut)
  }
}


# Derive function to relabel
recodeName=function(data){
  data=revalue(data, c("commonLiability"="Common liability",
                       "F1"="Common liability",
                       "commonLiability_filtered"= "Common liability",
                       "commonLiability_cancigalc_corRes" = "Common liability",
                       "commonLiability_filtered_cancigalc_corRes" = "Common liability (Qsnp filtered)",
                       "CigaretteDependency" = "Cigarette (dependence)",
                       "AlcoholDependency" = "Alcohol (dependence)",
                       "DrinksPerWeek" = "Alcohol (frequency)",
                       "CigarettesPerDay" = "Cigarette (frequency)",
                       "CannabisUseFrequency"= "Cannabis (frequency)",
                       "CannabisUseDisorder"= "Cannabis (dependence)"))
  return(data)
}

recodeNameBreak=function(data){
  data=revalue(data, c("commonLiability"="Common \n liability",
                       "F1"="Common \n liability",
                       "commonLiability_filtered"= "Common \n liability",
                       "CigaretteDependency" = "Cigarette \n (dependence)",
                       "AlcoholDependency" = "Alcohol \n (dependence)",
                       "DrinksPerWeek" = "Alcohol \n (frequency)",
                       "CigarettesPerDay" = "Cigarette \n (frequency)",
                       "CannabisUseFrequency"= "Cannabis \n (frequency)",
                       "CannabisUseDisorder"= "Cannabis \n (dependence)"))
  return(data)
}


# Rename label (add line breaks)
swr = function(string, nwrap) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)


# Check of order makes sense
orderLabels=function(variable){
  variableOut <- factor(variable,
                        levels = c("commonLiability",
                                   "commonLiability_filtered",
                                   "CigarettesPerDay",
                                   "CigaretteDependency",
                                   "DrinksPerWeek" ,
                                   "AlcoholDependency",
                                   "CannabisUseFrequency",
                                   "CannabisUseDisorder"))
  return(variableOut)
}


# Filter SNPs based on Q-statistic
filterQtest=function(df, name){
  print(paste0("Iteration: for ",name))
  if (name=="commonLiability") {
    print("Common factor")
    df=subset(df,   Q_chisq_pval >= 5e-8)
  } else {
    print("Other substance use phenotypes")
    df=subset(df,   Q_chisq_pval < 5e-8)
  }
  return(df)
}

# Select top SNPs according to p-value
topSNPs=function(data, NtopSNPs, P=P){
  data=data %>% top_n(NtopSNPs, -(as.numeric(P)))
  return(data)
}



# Function to format p-values and number with many decimals
formatNum=function(vecP=NULL, vecNum=NULL, decimal=3, df){

  if(length(vecP)!=0){

    for (i in 1:length(vecP)) {
      print(paste0("Format ",vecP[i]))
      df[[vecP[i]]]=formatC(df[[vecP[i]]], digits=decimal) # Format p-values
    }
  }
  if(length(vecNum)!=0){
    for (i in 1:length(vecNum)) {
      print(paste0("Round ",vecNum[i]))
      df[[vecNum[i]]]=round( df[[vecNum[i]]], decimal)
    }
  }
  return(df)
}


# Order list according to phenotype list
reorderPheno=function(listIn, listOut=list()){
  for(i in 1:length(phenoNames)){
    print(phenoNames[i])
    listOut[[i]] = listIn[[phenoNames[i]]]
  }
  names(listOut)=phenoNames
  return(listOut)
}


filterGenesEQTL=function(df){
  # Remove duplicates (eQTLs are linked to different tissues, but do not need to include here)
  df_out=subset(df,   duplicated(paste0(df$query_term, "_", df$gene))==FALSE)
  return(df_out)
}



##################################################################
# === Read in results =======================================
##################################################################


# +++++++++ Process output  +++++++++
strucmodel="cancigalc_constr_corRes"
GWASsumStats=readRDS(paste0(HOME,"/output/rds/", strucmodel, "_Sumstats.rds"))
postGWA=list()
clumpedSNPsList=list()
clumpedSNP_eQTLList=list()
phenoNames=names(GWASsumStats)


# Identify phenotype corresponding to the lowest p-value
# Standardize beta
standardBeta=function(df){
  print("Derive standardized beta")
  zscore = df$BETA / df$SE
  N=df$Neff
  df$BETA_STD = zscore / sqrt(N)  # Use effective sample size
  df$BETA_STD_pos=sqrt((df$BETA_STD)^2) # make all beta estimates positive
  print("Get standard error")
  VAR_STD = 1/N
  df$SE_STD = sqrt(VAR_STD)
  return(df)
}
GWASsumStats=lapply(GWASsumStats, function(x) standardBeta(x))
str(GWASsumStats)


# ================ Process GWAS input
for(i in 1:length(phenoNames)){
  print(paste0("***START PROCESSING:  ", phenoNames[i], "****"))
  sumStats=GWASsumStats[[phenoNames[i]]]

  print("Select GWA significant SNPs")
  clump_input=subset(sumStats, P<5e-8)

  if(NROW(clump_input)==0){
    print("No significant SNPs in GWAS")
    clumpedSNPsSelected=clump_input
    SNPs_sig_clumped=NROW(clumpedSNPsSelected)
  } else {
    write.table(clump_input,
                file= paste0(HOME, "/data/clump/",phenoNames[i], ".pvalsfile"),
                sep="\t",
                row.names = FALSE,
                col.names=T,
                quote=F)

    system(paste0(HOME, "/data/clump/plink --bfile ", HOME, "/data/clump/g1000_eur --clump ", HOME, "/data/clump/",  phenoNames[i], ".pvalsfile --clump-snp-field SNP --clump-field P --out ", HOME, "/data/clump/", phenoNames[i], ".pvalsfile --clump-kb 250 --clump-r2 0.1 --clump-p1 1"))
    print("Remove SNP file")
    clumped=NULL
    clumped=read.table(paste0(HOME, "/data/clump/", phenoNames[i], ".pvalsfile.clumped"), header = TRUE)
    system(paste0("rm ",HOME, "/data/clump/",phenoNames[i], ".pvalsfile"))

    clumpedSNPsSelected=subset(sumStats, SNP %in% unique(  clumped$SNP))
    print(paste0("Numbr of SNPs included after clumping ", NROW(clumped) ) )
    SNPs_sig_clumped=NROW(clumped)
  }

  if((NROW(clumpedSNPsSelected)==0)){
    clumpedSNPsSelectedTEMP <- data.frame(matrix(ncol = length(sumStats), nrow = 1))
    colnames(clumpedSNPsSelectedTEMP) <- colnames(sumStats)
    clumpedSNPsSelected=clumpedSNPsSelectedTEMP

  }
  SNPs_tot=NROW(sumStats)
  SNPs_sig_tot=NROW(subset(sumStats, P<5e-8)) # get total number of significant SNPs (before clumping)
  Ntotal=round(mean(sumStats$N, is.na=T), 0)

  # ==== eQTL Mapping
  print("eQTL mapping")
  #) Qtlizer to map genes to eQTL
  clumpedSNP_eQTL_query= get_qtls(clumpedSNPsSelected$SNP, corr = NA, max_terms = 5, ld_method = "r2",
                                  ref_version = "hg19", return_obj = "dataframe")

  if((length(clumpedSNP_eQTL_query$query_type)>0)) {
    if((clumpedSNP_eQTL_query$query_term[1]=='na')) {
      clumpedSNP_eQTL_query=subset(clumpedSNP_eQTL_query, query_term!='na')
    }
  }

  if((length(clumpedSNP_eQTL_query$query_type)==0)  ){
    print("no significant eQTLs identified")
    namesDF=colnames(get_qtls("rs62325470", corr = NA, max_terms = 5, ld_method = "r2",
                              ref_version = "hg19", return_obj = "dataframe"))
    clumpedSNP_eQTL_query <- data.frame(matrix(ncol = length(namesDF), nrow = 1))
    colnames(clumpedSNP_eQTL_query) <- namesDF
  }
  if((length(clumpedSNP_eQTL_query$query_type)>0)){
    print("significant eQTLs identified")
    # Note: Distance = distnace between index/proxy variant and gene in kilobases
    # 1 Mb (megabase) =  threshold for cis-effects
    # 1 MB = 1000 KB
    clumpedSNP_eQTL_query=subset(clumpedSNP_eQTL_query, distance<=1000)
    clumpedSNP_eQTL_query=subset(clumpedSNP_eQTL_query, sign_info=="FDR<5%" | sign_info=="FWER<5%" ) # selecg cis variants and those significant after FDR correction [eQTL mapping will map SNPs to genes which likely affect expression of those genes up to 1 Mb (cis-eQTL).FUMA: only eQTLs with FDR ≤ 0.05 will be used.]
  }

  # Merge clumped with eQTL
  clumpedSNP_eQTL= merge(clumpedSNP_eQTL_query, clumpedSNPsSelected, by.x="query_term",by.y="SNP", all.y=TRUE, suffixes = c("", "y"))

  print("add gene description")
  if( (NROW(na.omit(clumpedSNP_eQTL_query)))==0){ # if no eQTLs, create empty dataframe
    clumpedSNP_eQTL_desc=as.data.frame(matrix(ncol=2, nrow=1))
    colnames(clumpedSNP_eQTL_desc)=c("description", "input")
  } else {
    # Add gene description
    description =gconvert(query = unique(clumpedSNP_eQTL$gene) , organism = "hsapiens",
                          target="ENSG", mthreshold = Inf, filter_na = TRUE)

    clumpedSNP_eQTL_desc=subset(description, select = c(input, description))
  }

  clumpedSNP_eQTL=merge(clumpedSNP_eQTL, clumpedSNP_eQTL_desc, by.x="gene", by.y="input",all.x=TRUE, suffixes = c("", "y"))
  clumpedSNP_eQTL$p_eqtl=clumpedSNP_eQTL$p
  clumpedSNP_eQTL$p=NULL
  # Do not include p-value as this corresponds to expression in specific tissue (results are not included)

  # ===== GProfiler/Phenoscanner: Add description for genes
  print("Gene annotation using gprofiler2") # CHECK TUTORIAL: https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html

  # Phenoscanner
  print("Run phenoscanner")
  phenoQuery=clumpedSNPsSelected$SNP

  if( is.na(phenoQuery)==TRUE){
    print("Skip phenoscanner (no significant hits)")
    phenoQuery="rs62325470"
    pheno <- phenoscanner(snpquery=phenoQuery, build=37, catalogue=c("GWAS"))$snps
    phenoComplete=as.data.frame(matrix(ncol=length(colnames(pheno)), nrow=1))
    colnames(phenoComplete)=colnames(pheno)
    pheno=phenoComplete
    geneEnsembl=data.frame(target=NA, name=NA, description=NA)

  } else {
    pheno <- phenoscanner(snpquery=phenoQuery, build=37, catalogue=c("GWAS"))$snps

    pheno$hgnc=revalue(pheno$hgnc, c("-" = NA) )
    print("Add Ensembl ID and description")
    phenoComplete <- data.frame(hgnc=na.omit(pheno$hgnc)) # remove gene rows with missing values
    # Add description
    geneEnsembl=gconvert(query = phenoComplete$hgnc, organism = "hsapiens",
                         target="ENSG", mthreshold = Inf,
                         filter_na = TRUE)[c("target", "name", "description")]  %>% distinct(name, .keep_all = TRUE)
  }

  # Merge pheno and gprofiler results
  geneEnsemblPheno=merge(pheno, geneEnsembl, by.x = "hgnc", by.y = "name", all.x = T, suffixes = c("", ".y"))

  # Add to first dataframe
  clumpedSNPs=merge(clumpedSNPsSelected, geneEnsemblPheno, by.x = "SNP", by.y = "snp", all.x = TRUE, suffixes = c("", "y"))
  clumpedSNPs$P_logp = -log10(clumpedSNPs$P)
  print(paste0("Number of included SNPs after clumping: ", NROW(na.omit(clumpedSNPs)) ))


  print("Add info for Qsnp statistics")
  clumpedSNPs$chi_dich=ifelse(as.numeric(clumpedSNPs$Q_chisq_pval) < 5e-8, "Q_sig", "Q_ns" )
  clumpedSNPs$SNPs_tot=SNPs_tot
  clumpedSNPs$SNPs_sig_tot=SNPs_sig_tot
  clumpedSNPs$SNPs_sig_clumped=SNPs_sig_clumped
  clumpedSNPs$N=Ntotal

  clumpedSNP_eQTL$chi_dich=ifelse(as.numeric(clumpedSNP_eQTL$Q_chisq_pval) < 5e-8, "Q_sig", "Q_ns" )
  clumpedSNP_eQTL$SNPs_tot=SNPs_tot
  clumpedSNP_eQTL$SNPs_sig_tot=SNPs_sig_tot
  clumpedSNP_eQTL$SNPs_sig_clumped=SNPs_sig_clumped
  clumpedSNP_eQTL$N=Ntotal

  # Combine all datasets (not filtered)
  listOut=list(sumStats, clumpedSNPs, clumpedSNP_eQTL)
  names(listOut)=c("sumStats", "clumpedSNPs", 'clumpedSNP_eQTL')

  # Save loop estimates in list
  postGWA[[i]]=listOut

}
names(postGWA)=phenoNames

# Extract sumstats file
sumStatsList=lapply(postGWA, `[[`, "sumStats")
names(sumStatsList)=phenoNames

##################################################################
# === Quick summary table of GWA results =========================
##################################################################

# Function to select only unique SNPs
processSNPs=function(df){
  # Remove SNPs that exist multiple times (eg when linking to different genes)
  df$ID_GENE=paste0(df$SNP, " (", df$hgnc, ")")
  df_out=subset(df, duplicated(ID_GENE)==FALSE, select = c("SNP","CHR" ,"A1","A2" ,"BETA", "SE", "BETA_STD_pos", "BETA_STD", "SE_STD", "BP","P_logp", "P", "Q_chisq_pval", "chi_dich","hgnc","target" ,"ID_GENE", "consequence", "pos_hg19", "description","SNPs_tot" ,"SNPs_sig_tot", "SNPs_sig_clumped", "N") )
  print(paste0("Number SNPs included following gene annotation: ",   length(df_out$SNP)))
  print(paste0("Number of clumped SNPs: ", df$SNPs_sig_clumped[1]))
  return(df_out)
}

# FUNCTION: Get number of significant hits (shared vs non-shared)
numberShared=function(df){
  df$shared_snps=NROW(subset(df, chi_dich=="Q_ns"))
  df$non_shared_snps=NROW(subset(df, chi_dich=="Q_sig"))
  df$minP=formatC(min(df$P, na.rm = TRUE), digits=2)
  df$GWAS=NA
  df_out= subset( df[1,], select=c(GWAS,
                                   N,
                                   SNPs_tot,
                                   SNPs_sig_clumped,
                                   minP,
                                   shared_snps,
                                   non_shared_snps))
  return(df_out)
}

# Get results for clumped data
postGWA_clump=lapply(postGWA, `[[`, "clumpedSNPs") # Select only results from clumping
postGWA_clump=lapply(postGWA_clump, function(x) processSNPs(x))
GWA_shortSumList=lapply(postGWA_clump, function(x) numberShared(x)) # Check number of shared vs non-shared SNPs
GWA_shortSum=do.call(rbind, reorderPheno(GWA_shortSumList))
GWA_shortSum$GWAS=recodeName(names(postGWA_clump))
GWA_shortSum$minP=ifelse(GWA_shortSum$minP=="Inf", NA, GWA_shortSum$minP)
colnames(GWA_shortSum)=c("GWAS", "N (sample)", "included SNPs", "number of LD-independent genome-wide SNPs", "smallest p-value", "SNPs (shared)", "SNPs (non-shared)")



##################################################################
# === PhenoScanner ===============================================
##################################################################

runPheno=function(data, listDF=list(), nTopSNPsPheno){
  SNPs=levels(as.factor(data$SNP))

  if((length(SNPs)==0)){
    print("No significant hits - skip phenoscanner")
    phenoName=c("SNP", "BETA", "P", "ID_GENE", "trait", "P_pheno")
    pheno_out <- data.frame(matrix(ncol = length(phenoName), nrow = 1))
    colnames(pheno_out) <- phenoName
  }

  if((length(SNPs)>0)){
    print("Run phenoscanner")
    for (i in 1:length(SNPs)) {
      print(paste0("Read in ", SNPs[i]))
      data_sub <-  subset(data, SNP==SNPs[i] )
      pheno=phenoscanner(data_sub$SNP, build=37, catalogue=c("GWAS"))$results
      # If no phenoscanner results
      if((length(pheno)==0)){
        phenoName=c("rsID", "trait", "P_pheno")
        pheno <- data.frame(matrix(ncol = length(phenoName), nrow = 1))
        colnames(pheno) <- phenoName
        pheno$rsID=data_sub$rsID
      } else {
        # Select only trait column
        pheno=subset(pheno, select = c(snp, trait, p))
        # Select according to highest p value estimate
        pheno=pheno %>% top_n(nTopSNPsPheno, -(as.numeric(p)))
        colnames(pheno)=c("rsID", "trait", "P_pheno")
      }
      pheno$snp
      listDF[[i]]=merge(data_sub, pheno, by.x="SNP", by.y="rsID", all.x = TRUE)
    }

    # Combine list
    pheno_out=do.call("rbind", listDF)
    pheno_out=subset(pheno_out, select = c("SNP", "BETA", "P", "ID_GENE", "trait", "P_pheno"))
  }

  return(pheno_out)
}

# Select eqtl data
postGWA_eQTL=lapply(postGWA, `[[`, "clumpedSNP_eQTL") # Select results from clumping
names(postGWA_eQTL)=phenoNames
postGWA_eQTL=lapply(postGWA_eQTL, function(x) filterGenesEQTL(x)) # filter out only genes (ignore tissues)

# Filter SNPs (QNSP)
postGWA_clump_Qfilter=imap(postGWA_clump, function(x, y) filterQtest(x,y))
postGWA_eQTL_Qfilter=imap(postGWA_eQTL, function(x, y) filterQtest(x,y))

# Select only top SNPs
topNGWA=5
postGWA_clump_Qfilter_topSNP=lapply(postGWA_clump_Qfilter, function(x) topSNPs(x, NtopSNPs=topNGWA))
# Run phenoscanner on top SNPs
nTopSNPsPheno=2
postGWA_clumpTopPheno=lapply(postGWA_clump_Qfilter_topSNP, function(x) runPheno(x, nTopSNPsPheno=nTopSNPsPheno))
names(postGWA_clumpTopPheno)=phenoNames
postGWA_clumpPhenoDF=mergeLongFormat(postGWA_clumpTopPheno, label=recodeName(phenoNames), cleanLabel = TRUE)



##################################################################
# === Manhattan plot =============================================
##################################################################

# Function to prepare data for manhattan plot
dataManHplot=function(nameGWA){

  plotDF=data.frame(
    SNP=postGWA[[nameGWA]]$sumStats$SNP,
    Chromosome = postGWA[[nameGWA]]$sumStats$CHR,
    Position = postGWA[[nameGWA]]$sumStats$BP,
    outName = postGWA[[nameGWA]]$sumStats$P)

  colnames(plotDF)[4] = recodeName(nameGWA)

  manHplotInput=postGWA_clump[[nameGWA]]
  # Add row break in resID/gene name variable
  SNPsComFacText=as.character(paste0(manHplotInput$SNP, "\n(" ,manHplotInput$hgnc, ")"))
  # Derive colour scheme for mediated vs non-mediated SNPs
  SNPsComFac=postGWA_clump[[nameGWA]]$SNP
  SNPsComFacCol=SNPsComFac
  SNPsComFacCol[manHplotInput$Q_chisq_pval < 5e-8] = "red4"
  SNPsComFacCol[manHplotInput$Q_chisq_pval >= 5e-8 ] = "navy"
  # Save in list
  listGWA_manhattanOut=list(plotDF=plotDF, SNPs=SNPsComFac, SNPsText=SNPsComFacText, SNPsCol=as.character(SNPsComFacCol))
  return(listGWA_manhattanOut)
}

# Function to plot manhattan plot
ManHplot=function(list, width, height, type, ylimit=NULL){
  plotOut=CMplot(list[["plotDF"]],
                 type="p",
                 plot.type=type,
                 LOG10=TRUE,
                 file="jpg",
                 dpi=300,
                 main="",
                 cex=0.5, # size of the points (all)
                 amplify=TRUE,
                 signal.cex=0.5, # size of the points (significant)
                 threshold=c(5e-8),
                 threshold.col=c("black","grey"),
                 highlight=unlist(list[["SNPs"]]),
                 highlight.text=unlist(list[["SNPsText"]]),
                 highlight.text.cex=0.7,
                 highlight.text.col=unlist(list[["SNPsCol"]]),
                 highlight.col= "black",
                 col=c("grey87","grey60"),
                 ylab=expression(-log[10](italic(p))),
                 highlight.text.font=3,
                 file.output=FALSE,
                 width=width,
                 ylim=ylimit,
                 height=height)
  return(plotOut)
}

height=18
width=37
# Set parameters
ManHplot_width=37
ManHplot_height=18

# Margins
top=2
left=6
right=2
bottom=4

#max(-log10(postGWA[["commonLiability"]]$sumStats$P))

# y lim for main manuscript
ylimit=NULL
i=1
# y lim for supplement
ylimit=c(0,100)

# manhattan plot
for(i in 1:length(phenoNames)){
  # Prepare data
  dataManHplot_com=dataManHplot(phenoNames[i])
  # plot output
  jpeg(file=paste0(HOME,"/results/figures/ManHplot_",phenoNames[i], ".jpeg"),  width=ManHplot_width,height=ManHplot_height, units = "cm", res=1000)
  par(mar = c(bottom, left, top, right))
  dataManHplot_comOut=ManHplot(dataManHplot_com, width=ManHplot_width, height=ManHplot_height, type="m", ylimit=ylimit)
  dev.off()
}

# QQplot
QQplot_width=40
QQplot_height=30
for(i in 1:length(phenoNames)){
  dataQQplot_com=dataManHplot(phenoNames[i])
  jpeg(file=paste0(HOME,"/results/figures/QQplot_",phenoNames[i], ".jpeg"),  width=QQplot_width,height=QQplot_height, units = "cm", res=1000)
  par(mar = c(6, 6, 6, 6))
  QQplot_comOut=ManHplot(dataQQplot_com, width=QQplot_width, height=QQplot_height, type="q", ylimit= c(0,100))
  dev.off()
}

##################################################################
# === Shared and non-shared genetic variants =====================
##################################################################

# Arrange according to standardized beta
bindTogether_clump=lapply(postGWA_clump, function(x) arrange(x, -BETA_STD_pos))
str(bindTogether_clump)

# Bind as dataframe
postGWA_clumpDF=do.call(rbind, bindTogether_clump)

# Match relvant SNPs to estimates from GWAS output
sumStatsRaw=lapply(postGWA, `[[`, "sumStats") # Select raw sumstats files
names(sumStatsRaw)=phenoNames
postGWA_clumpGenes=list()


for(i in 1:length(phenoNames)){
  print(paste0("Start processing of ", phenoNames[[i]]))

  print("Select original files")
  df1=sumStatsRaw[[phenoNames[[i]]]]
  df1$chi_dich[df1$Q_chisq_pval < 5e-8] = "Q_sig"
  df1$chi_dich[df1$Q_chisq_pval >= 5e-8] = "Q_ns"
  df1$P_logp = -log10(df1$P)
  print("Select mapped genes")
  postGWA_clumpDF=subset(postGWA_clumpDF, is.na(SNP)!=T & is.na(hgnc)!=T )
  postGWA_clumpDF$ID_GENE=paste0(postGWA_clumpDF$SNP, " (", postGWA_clumpDF$hgnc, ")")
  postGWA_clumpDF_clean=subset(postGWA_clumpDF, duplicated(ID_GENE)==FALSE)

  print("Merge with original files")
  df1_selected = subset(df1, SNP %in% unique(postGWA_clumpDF$SNP), select=c("SNP","BETA","Q_chisq_pval", "chi_dich", "P_logp", "P", "CHR", "BETA_STD", "BETA_STD_pos", "SE_STD") )
  df1_selected=merge(df1_selected, postGWA_clumpDF_clean, by="SNP", all.y = T, suffixes = c("", ".y"))
  df1_selected=subset(df1_selected, select= c(SNP, BETA, CHR, Q_chisq_pval,chi_dich,  P_logp, P,ID_GENE, hgnc, BETA_STD, BETA_STD_pos, SE_STD))

  # Extract lead SNPs and label in dataframe
  df2_clump=postGWA_clump[[phenoNames[[i]]]]
  df1_selected$leadSNP= ifelse(df1_selected$SNP %in% df2_clump$SNP , "*", "")
  postGWA_clumpGenes[[i]]=df1_selected
}


# Combine in one dataframe
bindTogether = mapply(`[<-`, postGWA_clumpGenes, 'Label', value = names(postGWA_clump), SIMPLIFY = FALSE)
postGWA_clump_genesDF=do.call(rbind, bindTogether)

# Order according to phenotype
postGWA_clump_genesDF=postGWA_clump_genesDF[order(match(postGWA_clump_genesDF$ID_GENE, postGWA_clumpDF$ID_GENE)), ]
postGWA_clump_genesDF$order=seq(1:length(postGWA_clump_genesDF$SNP))

# Add row breaks
postGWA_clump_genesDF$Label_clean=as.factor(recodeName(orderLabels(as.factor(postGWA_clump_genesDF$Label))))

# Get confidence interval
postGWA_clump_genesDF$beta_upper_CI=postGWA_clump_genesDF$BETA_STD_pos + 1.96 * postGWA_clump_genesDF$SE_STD
postGWA_clump_genesDF$beta_lower_CI=postGWA_clump_genesDF$BETA_STD_pos - 1.96 * postGWA_clump_genesDF$SE_STD


# FUNCTION FOR PLOT
comparativePlot=function(df, colGroup){
  plot=ggplot(df, aes(x = reorder(ID_GENE, -order), y = BETA_STD_pos,  group = colGroup, fill =colGroup )) +
    geom_bar(stat = "identity") + coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size=6),
          legend.position="none",
          plot.caption = element_text(hjust = 0))  +
    geom_text(mapping = aes(x=reorder(ID_GENE, -order), y=BETA_STD_pos, label = leadSNP),  size = 3, colour = "black") +
    facet_wrap(~ as.factor(Label_clean),   ncol = length(levels(as.factor(df$Label_clean))), nrow = NULL) +
    geom_linerange(aes(x = ID_GENE, ymin = beta_lower_CI, ymax = beta_upper_CI),
                   lwd = 0.1, position = position_dodge(width = 1/2), colour = "black") +
    scale_y_continuous(limits=c(round(min(postGWA_clump_genesDF$beta_lower_CI),2)-0.01, round(max(postGWA_clump_genesDF$beta_upper_CI),2),
                                breaks=c(0, 0.04, round(max(postGWA_clump_genesDF$beta_upper_CI),2) ), n.breaks = 3 )) +
    labs(y=expression((italic(beta)[std])), x="SNP (annotated gene)")
  return(plot)
}


# Add row breaks
postGWA_clump_genesDF$labelOrdered=orderLabels(postGWA_clump_genesDF$Label)
postGWA_clump_genesDF$Label_clean=as.factor(recodeNameBreak(orderLabels(as.factor(postGWA_clump_genesDF$labelOrdered))))


postGWA_clump_genesDF$chi_dich[postGWA_clump_genesDF$P >= 5e-8 ] = "ns"

#colourScheme=c("aquamarine4","tomato" ,"tomato", "red4","navy", "darkgreen", "grey")
colourScheme=c( "grey","navy", "red4")
clumpedMediatedPlottopPheno = comparativePlot(postGWA_clump_genesDF, colGroup=postGWA_clump_genesDF$chi_dich) +
  scale_fill_manual(name="", values = colourScheme)

print(clumpedMediatedPlottopPheno)


# Export
clumpedMediatedPlot_width=27
clumpedMediatedPlot_height=40
ggsave(paste0(HOME,"/results/figures/clumpedMediatedPlot.pdf"), clumpedMediatedPlottopPheno,  width = clumpedMediatedPlot_width, height = clumpedMediatedPlot_height, units = "cm", limitsize = TRUE)
ggsave(paste0(HOME,"/results/figures/clumpedMediatedPlot.svg"), clumpedMediatedPlottopPheno,  width = clumpedMediatedPlot_width, height = clumpedMediatedPlot_height, units = "cm", limitsize = TRUE, bg="white")





############################################
# ======== eQTL ============================
############################################

# ===== eQTL: SNPs to genes

# combine positional with eQTL
postGWA_pos_eqtl=list()
for(i in 1:length(phenoNames)){
  print(paste0("Merge positional and eQTL mapping for ", phenoNames[i]))
  df_clump=subset(postGWA_clump[[phenoNames[i]]], select=c(SNP, BP, hgnc, target, description, consequence))
  colnames(df_clump)=c("SNP", "BP_pos","gene_pos", "ens_pos", "desc_pos", "consequence_pos")
  df_eqtl=subset(postGWA_eQTL[[phenoNames[i]]], select=c(query_term, gene, distance , ensgid, description, P, Q_chisq_pval, chi_dich))
  colnames(df_eqtl)=c("SNP", "gene_eqtl", "distance","ens_eqtl", "desc_eqtl", "P_GWAS", "Q_chisq_pval", "chi_dich")
  postGWA_pos_eqtl[[i]]=merge(df_clump, df_eqtl, by="SNP", all = T)
}
str(postGWA_pos_eqtl)
names(postGWA_pos_eqtl)=phenoNames
```

# Pasthway analysis in PASCAL and DEPICT

[](#pasthway-analysis-in-PASCAL-and-DEPICT)

### PASCAL

#### Step 1: Format the summary statistic files

``` r
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
```

#### Step 2: Run PASCAL

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

#### Step 3: Process results from PASCAL

-   Script used on the cluster

``` r
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
```

-   Script used on a desctop computer after uploading results stored on
    the cluster:

``` r
################################################
# ======== PASCAL  results =====================
################################################

processPasc=function(df, list=list(), names, source){
  print(names)
  for ( i in 1:length(names) ) {
    dfIn=df[[i]]
    print(paste0("Process ", names[i]))
    dfIn = dfIn %>%
      mutate(p_corr = p.adjust(chi2Pvalue, method = "fdr", n = NROW(dfIn) ),
             p_corr_dich = ifelse(p_corr < 0.05 , "p_cor_sig", "p_cor_ns"),
             p_sig = ifelse(p_corr < 0.05 , "*", ""),
             nTests=NROW(dfIn),
             label = names[i],
             P_logp = -log10(chi2Pvalue),
             P_logp_cor = -log10(p_corr),
             Name = tolower(gsub("_", " ", Name)) ) %>%
      arrange(chi2Pvalue)

    dfIn$nTests=NROW(dfIn)
    dfIn$PathData=source
    dfIn$ID=paste0(dfIn$Name, "_", dfIn$PathData)
    list[[i]]=dfIn
  }
  return(list)
}

# ==== Add GO terms
data("AlvMac_allGO")
goDF=AlvMac_allGO
goDF$namespace_1003=NULL
colnames(goDF)=c("term", "pathAnno")

# ==== Add MP terms
#setwd(paste0(HOME, "/data/"))
#system("wget https://maayanlab.cloud/static/hdfs/harmonizome/data/mgimpo/attribute_list_entries.txt.gz")
mpDF=read.table(file= paste0(HOME, "/data/attribute_list_entries.txt.gz"),
                header = TRUE,
                sep="\t")
mpDF=subset(mpDF, select=c(MPID, Phenotype))
colnames(mpDF)=c("term", "pathAnno")

# ==== Add ENSEMBL terms
ensDF=read.table(file= paste0(HOME, "/programs/depict/data/inweb_mapping.tab"),
                 header = FALSE,
                 sep="\t")
colnames(ensDF)=c("term", "pathAnno")

# Add labels for remaining pathways
otherpfDF <- openxlsx::read.xlsx(
  xlsxFile = paste0(HOME, "/data/pathwayAnno.xlsx"), sheet = 1, skipEmptyRows = TRUE,
  detectDates = TRUE
)
allPath=as.data.frame(rbind(goDF, mpDF, ensDF, otherpfDF))

# Derive function
addGOlabel=function(df){
  df$Name_original=df$Name
  df$Name_rec=toupper(df$Name)
  df=merge(df, allPath, by.x="Name_rec", by.y="term", all.x=TRUE, suffixes=c("",".y") )
  # Get rows for ensemble
  rowNumberENS=grep("ens", df$Name, value=FALSE, perl=TRUE, ignore.case=TRUE)
  df$pathAnno=ifelse(is.na(df$pathAnno)==T,   df$Name_original, df$pathAnno)
  for(i in 1:length(rowNumberENS)){
    ens_term= toupper(df[rowNumberENS[i],]$Name)
    print(    paste0("Check ",i, ": ", ens_term))
    ensDF_select =gconvert(query = ens_term, organism = "hsapiens",
                           target="ENSG", mthreshold = Inf, filter_na = TRUE)

    if(NROW(ensDF_select)==0){
      print("Pathname not available")
      df[rowNumberENS[i],]$pathAnno=  df[rowNumberENS[i],]$Name
    } else {
      if((ensDF_select$target)=="nan" ) {
        print("Pathname not available")
        df[rowNumberENS[i],]$pathAnno=  df[rowNumberENS[i],]$Name
      }
      print("Pathname available")
      df[rowNumberENS[i],]$pathAnno=paste0(ensDF_select$name, " (", ensDF_select$description, ")")
    }


  }

  df$pathAnnoShort = str_trunc(df$pathAnno , 30, "right") # Shorten label
  return(df)
}


# Import rds file (MSigDB)
listPascalReact_raw=readRDS(paste0(HOME, "/output/rds/PascalReact.rds"))
namesPascal=names(listPascalReact_raw)
listPascalReact_all = processPasc(df=listPascalReact_raw, list=list(), names=names(listPascalReact_raw), source="MSigDB")
names(listPascalReact_all)=namesPascal
nSetsReact=listPascalReact_all[[1]]$nTests[1]

# Import rds file (DEPICT)
listPascalDepict_raw=readRDS(paste0(HOME, "/output/rds/PascalDepict.rds"))
listPascalDepict_all = processPasc(df=listPascalDepict_raw, list=list(), names=names(listPascalDepict_raw), source="DEPICT")
names(listPascalDepict_all)=namesPascal
str(listPascalDepict_raw)
nSetsDepict=listPascalDepict_all[[1]]$nTests[1]


# Bind results (all paths)
listPascalAll=list()
for(i in 1:length(namesPascal)){
  print(namesPascal[i])
  df1=as.data.frame(listPascalReact_all[[namesPascal[i]]])
  df2=as.data.frame(listPascalDepict_all[[namesPascal[i]]])
  df=rbind(df1, df2)
  print("Number of significant pathways")
  print( table(df$p_corr_dich))
  listPascalAll[[i]]=df
}
names(listPascalAll)=namesPascal
str(listPascalAll)
listPascalAll[["commonLiability"]]=NULL # remove results for common liability GWA that is not filtered based on Qsnp
str(listPascalAll)

# Select all significant pathways
bindTogetherPASCAL = mapply(`[<-`, listPascalAll, 'Label', value = names(listPascalAll), SIMPLIFY = FALSE)
bindTogetherPASCAL_sig=unique(subset(do.call(rbind, bindTogetherPASCAL), p_corr_dich=="p_cor_sig"))
#bindTogetherPASCAL_sig=unique(subset(do.call(rbind, bindTogetherPASCAL), chi2Pvalue < 0.05))
str(bindTogetherPASCAL_sig)

# Add label
listPascalSelected_anno=addGOlabel(bindTogetherPASCAL_sig)
listPascalSelected_annoSub=subset(listPascalSelected_anno, select=c(ID, pathAnno, pathAnnoShort) )
str(listPascalSelected_annoSub)

# Select the same set of paths from each GWA
listPascalSelected=list()
for (i in 1:length(listPascalAll)) {
  df1=listPascalAll[[i]]
  listPascalSelected[[i]]=merge(listPascalSelected_annoSub, df1, by="ID", all.x=T)
}
names(listPascalSelected)=names(listPascalAll)
str(listPascalSelected)

# Convert to dataframe
library(plyr)
listPascalSelected_DF=mergeLongFormat(listPascalSelected, names(listPascalAll), cleanLabel = FALSE)

# Revalue according to class
listPascalSelected_DF$label=as.factor(orderLabels(as.factor(listPascalSelected_DF$label)))
listPascalSelected_DF$Class=revalue(listPascalSelected_DF$label, c("commonLiability_filtered"="Common \n liability",
                                                                   "CigarettesPerDay" = "Cigarette",
                                                                   "CigaretteDependency" = "Cigarette",
                                                                   "DrinksPerWeek" = "Alcohol",
                                                                   "AlcoholDependency" = "Alcohol",
                                                                   "CannabisUseFrequency" = "Cannabis",
                                                                   "CannabisUseDisorder" = "Cannabis"))
listPascalSelected_DF$Label=recodeNameBreak(listPascalSelected_DF$label)

head(listPascalSelected_DF)
listPascalSelected_DF$Class = factor(listPascalSelected_DF$Class, levels=c('Common \n liability',
                                                                           'Alcohol',
                                                                           'Cigarette',
                                                                           "Cannabis"))
levels(as.factor(listPascalSelected_DF$label))



# Create vector for order (based on common liability p-values)
df_order=listPascalAll[["commonLiability_filtered"]] %>% arrange(-P_logp)
df_order$order=seq(1:length(df_order$Name))
df_order=subset(df_order, select=c(ID, order))
listPascalSelected_DF_ordered=merge(listPascalSelected_DF, df_order, by="ID", all.x = T)

# derive log scale
listPascalSelected_DF_ordered$p_corr_log= -log10(listPascalSelected_DF_ordered$p_corr)

# Remove duplicates
listPascalSelected_DF_ordered=as.data.frame(listPascalSelected_DF_ordered %>%
                                              group_by(label) %>%
                                              distinct(ID, .keep_all= TRUE))

# Select top paths
ntopPascal=15
listPascalSelectedTOP=lapply(listPascalSelected, function(x) top_n(x, ntopPascal,  -(as.numeric(chi2Pvalue))) )
PathsSelection=unique(do.call("rbind",lapply(listPascalSelectedTOP, '[', 'Name'))$Name)
listPascalSelected_DF_topSNP=subset(listPascalSelected_DF_ordered, Name %in% unique(  PathsSelection))
levels(as.factor(listPascalSelected_DF_topSNP$pathAnnoShort))

# Heatmap (MAIN MANUSCRIPT)
PASCAL_heatmapMain=
  listPascalSelected_DF_topSNP %>%
  mutate(name = fct_reorder2(pathAnnoShort, as.factor(label), order)) %>%
  ggplot(
    aes(y=name,
        x= Label,
        fill= P_logp,
        label = p_sig)) +
  geom_tile() +
  facet_grid(cols = vars(Class), scales = "free") +

  scale_fill_gradient(low="white", high="blue",  na.value = "grey50", name=expression(-log[10](italic(p))) ) +
  theme_minimal() + theme(axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y = element_text(size=9)) +
  geom_text(size=5, color="black",  hjust = 0, nudge_x = 0, nudge_y = 0)



# Heatmap (SUPPLEMENT)
PASCAL_heatmapSupp=
  listPascalSelected_DF_ordered %>%
  mutate(name = fct_reorder2(pathAnnoShort, label, order)) %>%
  ggplot(
    aes(y=name,
        x= Label,
        fill= p_corr_log,
        label = p_sig)) +
  geom_tile() +
  facet_grid(cols = vars(Class), scales = "free", space = "free") +

  scale_fill_gradient(low="white", high="blue",  na.value = "grey50", name=expression(-log[10](italic(p))) ) +
  theme_minimal() + theme(axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y = element_text(size=7))  +
  geom_text(size=5, color="black",  hjust = 0, nudge_x = 0, nudge_y = -0.3)
print(PASCAL_heatmapSupp)


# Export plot
pascalPlot_width=40
pascalPlot_height=80
ggsave(paste0(HOME,"/results/figures/PASCAL_suppFig.pdf"), PASCAL_heatmapSupp,  width = pascalPlot_width, height = pascalPlot_height, units = "cm", limitsize = TRUE)


# Combine DEPUCT and Pascal
pascalPlot_comb=ggarrange(DEPICT_heatmap,
                          PASCAL_heatmapMain,
                          common.legend = TRUE,
                          legend = "bottom",
                          ncol=2,
                          labels = c("A", "B"))

# Export plot
pascalPlot_width=44
pascalPlot_height=35
ggsave(paste0(HOME,"/results/figures/pascalPlot_comb.pdf"), pascalPlot_comb,  width = pascalPlot_width, height = pascalPlot_height, units = "cm", limitsize = TRUE)
ggsave(paste0(HOME,"/results/figures/pascalPlot_comb.svg"), pascalPlot_comb,  width = pascalPlot_width, height = pascalPlot_height, units = "cm", limitsize = TRUE)
```

### DEPICT

Note that the analysis done using DEPICT is included in the R Script
[CommonLiabAddiction.R](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/CommonLiabAddiction.R)
and was performed on a desktop computer

    ################################################
    # ======== DEPICT  =============================
    ################################################
    # ======== RUN DEPICT ON SUBSET OF SNPs ===================================

    # Exclude SNPs for common liability based on Qsnp
    prepareDEPdata=function(df, name){
      print(paste0("Iteration: for ",name))
      if (name=="commonLiability") {
        print("Common factor")
        dfout=subset(df,   P<=5e-5 & Q_chisq_pval >= 5e-8, select=c(SNP, P))
      } else {
        print("Other substance use phenotypes")
        dfout=subset(df,   P <= 5e-5, select=c(SNP, P))
      }
      print(paste0("Number of incldued SNPs: ", NROW(dfout)))
      return(dfout)
    }

    GWASsumStatsDEPICT=imap(sumStatsList, function(x, y) prepareDEPdata(x,y))

    # delete existing files
    system(paste0("rm ",HOME, "/programs/depict/SNPfile/*"))
    dir(paste0(HOME, "/programs/depict/SNPfile/"))

    # select SNP data and clump
    for(i in 1:length(phenoNames)){
      print(paste0("***START clumping:  ", phenoNames[i], "****"))
      sumStatsDepict=GWASsumStatsDEPICT[[phenoNames[i]]]

      print(paste0("Number of included SNPs for ", phenoNames[i], ": ", NROW(sumStatsDepict)))
      write.table(sumStatsDepict,
                  file= paste0(HOME, "/programs/depict/SNPfile/sumStatsDepict.", phenoNames[i]),
                  sep="\t",
                  row.names = FALSE,
                  col.names=TRUE,
                  quote=FALSE)

      print("Clump files")
      system(paste0(HOME, "/data/clump/plink --bfile ", HOME, "/data/clump/g1000_eur --clump ", HOME, "/programs/depict/SNPfile/sumStatsDepict.", phenoNames[i]," --clump-snp-field SNP --clump-p1 1 --clump-p2 1 --clump-field P --out ", HOME, "/programs/depict/SNPfile/sumStatsDepict_",phenoNames[i], " --clump-kb 500 --clump-r2 0.05"))

      SNPsClump=read.table(file= paste0(HOME, "/programs/depict/SNPfile/sumStatsDepict_",phenoNames[i], ".clumped"),
                           header = TRUE)

      SNPsClumpSel=subset(sumStatsDepict, SNP %in% unique(SNPsClump$SNP), select=c(SNP))

      print("Save clumped file formatted for DEPICT")
      write.table(SNPsClumpSel,
                  file= paste0(HOME, "/programs/depict/SNPfile/sumStatsDepictClumped.", phenoNames[i]),
                  sep="\t",
                  row.names = FALSE,
                  col.names=TRUE,
                  quote=FALSE)

    }

    # ============ RUN DEPICT
    # remove old files
    #system(paste0("rm ",HOME,  "/output/depict/*"))
    #system(paste0("rm ",HOME,  "/programs/depict/results/*"))

    for(i in 1:length(phenoNames)){
      SNPSepictClumped=read.table(file= paste0(HOME, "/programs/depict/SNPfile/sumStatsDepictClumped.", phenoNames[i]),
                                  header = TRUE)
      print(paste0("Number of included SNPs: ",   NROW(SNPSepictClumped)))

      write(paste0("Number of included SNPs for ",phenoNames[i],": ",   NROW(SNPSepictClumped)),
            file=paste0(HOME,  "/output/depict/DEPICT.log"),
            append=TRUE)

      print("Remove temporary files")
      setwd(paste0(HOME, "/programs/depict")) # has to be in the directory where the program is stored
      system(paste0("rm ",HOME, "/programs/depict/results/*")) # remove temp files
      system(paste0("chmod +x ",HOME, "/programs/depict/DepictGenSEM.py"))
      system(paste0("python2 ", HOME, "/programs/depict/DepictGenSEM.py ", phenoNames[i], " ", HOME))

      enrichDepict=read.table(file= paste0(HOME, "/programs/depict/results/sumStatsDepict_",phenoNames[i],"_tissueenrichment.txt"),
                              header = TRUE,
                              sep="\t")

      enrichDepict_sel=enrichDepict[,2:5]
      colnames(enrichDepict_sel)=c("PathName", "PathCat", "pvalDepict", "pvalDepict_fdr")


      system(paste0("rm ",HOME,  "/output/depict/", phenoNames[i], "_depictOut"))
      write.table(enrichDepict_sel,
                  file= paste0(HOME, "/output/depict/", phenoNames[i], "_depictOut"),
                  sep="\t",
                  row.names = FALSE,
                  col.names=TRUE,
                  quote=FALSE)

    }

    # ======= Process depict =======
    readinDepict=function(name){
      print(paste0("Read in ", name))
      dfOut=read.table(file=paste0(HOME, "/output/depict/", name, "_depictOut"),
                       header = TRUE,
                       sep="\t")
      dfOut$p_corr = p.adjust(dfOut$pvalDepict, method = "fdr", n =   NROW(dfOut) )
      dfOut$p_corr_dich = ifelse(dfOut$p_corr < 0.05 , "p_cor_sig", "p_cor_ns")
      return(dfOut)
    }
    depictOut=lapply(phenoNames, function(x) readinDepict(x))
    names(depictOut)=phenoNames

    # Select and bind together significant pathways
    bindTogetherDEPICT = mapply(`[<-`, depictOut, 'Label', value = phenoNames, SIMPLIFY = FALSE)
    bindTogetherDEPICT_names=subset(do.call(rbind, bindTogetherDEPICT), pvalDepict<0.05)$PathName

    # Select paths
    listDepictTop=lapply(depictOut, function(x) subset(x, PathName %in% bindTogetherDEPICT_names) )

    # Select names of n paths
    PathsSelection=unique(do.call("rbind",lapply(listDepictTop, '[', 'PathName'))$PathName)
    length(PathsSelection)
    # Select the same set of paths from each GWA
    listDepictSelected=list()
    for (i in 1:length(phenoNames)) {
      df1=depictOut[[i]]
      df2=data.frame(PathName=PathsSelection)
      listDepictSelected[[i]]=merge(df2, df1, by="PathName", all.x=T)
    }
    names(listDepictSelected)=phenoNames

    # Convert to dataframe
    library(plyr)
    listdepictSelected_DF=mergeLongFormat(listDepictSelected, phenoNames, cleanLabel = FALSE)
    listdepictSelected_DF$p_sig=NA
    listdepictSelected_DF$p_sig[listdepictSelected_DF$pvalDepict_fdr=="No"] <- ""
    listdepictSelected_DF$p_sig[listdepictSelected_DF$pvalDepict_fdr=="Yes"] <- "*"

    # Create vector for order (based on common liability p-values)
    df_order=depictOut[["commonLiability"]] %>% arrange(pvalDepict)
    df_order$order=seq(1:length(df_order$PathName))
    df_order=subset(df_order, select=c(PathName, order))
    listdepictSelected_DF_ordered=merge(listdepictSelected_DF, df_order, by="PathName", all.x = T)
    # Function heatmap
    listdepictSelected_DF_ordered$p_corr_log= -log10(listdepictSelected_DF_ordered$pvalDepict)

    # Revalue according to class
    listdepictSelected_DF_ordered$Label=as.factor(orderLabels(as.factor(listdepictSelected_DF_ordered$Label)))
    listdepictSelected_DF_ordered$Cat=revalue(listdepictSelected_DF_ordered$Label, c("commonLiability"="Common \n liability",
                                                                                     "DrinksPerWeek" = "Alcohol",
                                                                                     "AlcoholDependency" = "Alcohol",
                                                                                     "CigarettesPerDay" = "Cigarette",
                                                                                     "CigaretteDependency" = "Cigarette",
                                                                                     "CannabisUseFrequency" = "Cannabis",
                                                                                     "CannabisUseDisorder" = "Cannabis"))
    listdepictSelected_DF_ordered$Cat = factor(listdepictSelected_DF_ordered$Cat, levels=c('Common \n liability',
                                                                                           'Alcohol',
                                                                                           'Cigarette',
                                                                                           "Cannabis"))

    # Order factor
    listdepictSelected_DF_ordered$Label=recodeNameBreak(listdepictSelected_DF_ordered$Label)


    # Plot results
    DEPICT_heatmap=ggplot(listdepictSelected_DF_ordered,
                          aes(y=reorder(PathName, -order),
                              x= Label,
                              fill= p_corr_log,
                              label = p_sig)) +
      geom_tile() +    theme_minimal() + theme(axis.title.x=element_blank(),
                                               axis.title.y=element_blank(),
                                               axis.text.y = element_text(size=9)) +
      facet_grid(cols = vars(Cat), scales = "free", space = "free") +
      scale_fill_gradient(low="white", high="blue",  na.value = "grey50", name=expression(-log[10](italic(p))) ) +
      geom_text(size=5, color="black",  hjust = 0, nudge_x = 0, nudge_y =-0.4)
    print(DEPICT_heatmap)

</br></br></br>

# LD score regression analysis including other traits

[](#ld-score-regression-analysis-including-other-traits)

#### Step 1: Perform LD score regression analysis

-   LD score regression analysis can be conducted executing the
    following R script, which uses LD score regression as implemented in
    GenomicSEM

-   the analysis was performed on a HPC cluster

-   Note: the script
    [‘ldscoreregressionAnalysis.R’](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/CommonLiabAddiction.R)
    will read in all the summary statistic files from
    [GenSEM\_GWAS.xlsx](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/blob/master/analysis/)
    that are selected to be included in the multivariate genome-wide
    association analysis. Selected are those that those with
    `directory == "ldsc"`

-   it will also read in the summary statistic files generated for the
    common liability

<!-- -->

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

#### Step 2: Process LD score regression analysis results

``` r
################################################
# ======== Plot LD score regression results ====
################################################

# Import rds file
LdScoreResultsAll=readRDS(paste0(HOME, "/output/rds/ldscExternalOutDF.rds"))
NROW(LdScoreResultsAll)

gwasSumStast=NULL
gwasSumStast=read.xlsx(paste0(HOME,"/analysis/GenSEM_GWAS.xlsx"), 1)
gwasSumStast=subset(gwasSumStast, directory!="exclude", select=c(category, label, label_clean))
NROW(gwasSumStast)

# Remove traits to be excluded
LdScoreResults=merge(LdScoreResultsAll, gwasSumStast, by.x = "trait",by.y= "label", all.y=T)

LdScoreResults$Class=revalue(LdScoreResults$category, c("MentalHealth"="Mental \n health",
                                                        "PsychSocial" = "Social \n psychological",
                                                        "Anthropometric" = "Anthro- \n pometric",
                                                        "Cognition" = "Cognition \n SES",
                                                        "gensem" = "Substance use \n (indicators)",
                                                        "SubstanceUse" = "Consumption \n (other)"))

levels(as.factor(LdScoreResults$category))
# Create class label

table(LdScoreResults$Class)
LdScoreResults$Class <- factor(LdScoreResults$Class, levels = c("Anthro- \n pometric","Social \n psychological" ,"Cognition \n SES" ,"Personality" ,"Mental \n health", "Brain", "Consumption \n (other)", "Substance use \n (indicators)" ))

# Get p-value
LdScoreResults$rg_upper_CI=LdScoreResults$rg + 1.96 * LdScoreResults$rg_se
LdScoreResults$rg_lower_CI=LdScoreResults$rg - 1.96 * LdScoreResults$rg_se
LdScoreResults$zscore=LdScoreResults$rg/LdScoreResults$rg_se
LdScoreResults$pval=2*pnorm(-abs(LdScoreResults$zscore))

# Remove rows
LdScoreResults=subset(LdScoreResults, trait!= "commonLiability_filtered_corResConstr" & trait!= "commonLiability_corResConstr" )
LdScoreResults$commonGWAS=recodeName(LdScoreResults$commonGWAS)
# include only qsnp filtered results
LdScoreResults=subset(LdScoreResults, commonGWAS=="commonLiability_filtered_cancigalc_constr_corRes")

# Correct for multiple testing
nTestLDscoreRegression=NROW(LdScoreResults)
LdScoreResults$pval_fdr=p.adjust(LdScoreResults$pval, method = "fdr", n = nTestLDscoreRegression)
LdScoreResults$pval_fdr_dich[LdScoreResults$pval_fdr<=0.05]="*"
LdScoreResults$pval_fdr_dich[LdScoreResults$pval_fdr>0.05]=""

# Categorize according to class and fdr correction
LdScoreResults$Class_dich=ifelse(LdScoreResults$Class=="Substance use \n (indicators)", "SubstanceUse", "Other")


# Plot results
ldscPlot=function(df, colCom="navy", twoCols=F){

  ldscPlot <- ggplot(df, aes(x=label_clean,
                             y=rg,
                             fill=Class_dich)) +
    geom_bar(stat="identity", position = position_dodge()) +
    geom_errorbar(aes(ymin=rg_lower_CI, ymax=rg_upper_CI), width=.2,
                  position=position_dodge(.9)) +
    coord_flip()
  if(twoCols==T){
    ldscPlot = ldscPlot + facet_grid(rows = vars(Class), cols = vars(commonGWAS), scales = "free", space = "free", margins=FALSE)
  }else{
    ldscPlot = ldscPlot + facet_grid(rows = vars(Class), scales = "free", space = "free", margins=FALSE)
  }
  ldscPlot = ldscPlot +
    theme_pubr() +
    scale_fill_manual(name="",
                      values = c(colCom, "grey")) +
    theme_minimal() +
    theme(strip.text.y = element_text(angle = 90),
          #panel.border = element_rect(color = "grey", fill = NA, size = 2),
          #panel.grid.major = element_line(colour="grey", size = (0.2)),
          axis.title.y=element_blank(), legend.position = "none") +

    guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
    geom_text(mapping = aes(x=label_clean, y=rg, label = pval_fdr_dich),
              size = 8, colour = "black", nudge_y = ifelse(df$rg > 0, 0.01, -0.01)) +
    scale_y_continuous(name =expression((italic(r[g]))),
                       limits=c(round(min(df$rg_lower_CI, na.rm = T),1),
                                round(max(df$rg_upper_CI, na.rm = T)+0.1,1)))
  return(ldscPlot)
}


# PLot results
ldscPlot_filter=ldscPlot(LdScoreResults, twoCols=F)
print(ldscPlot_filter)

# Save plot
ldscPlot_width=16
ldscPlot_height=27
ggsave(paste0(HOME,"/results/figures/PlotLDScore.pdf"), ldscPlot_filter,  width = ldscPlot_width, height = ldscPlot_height, units = "cm", limitsize = FALSE)
ggsave(paste0(HOME,"/results/figures/PlotLDScore.svg"), ldscPlot_filter,  width = ldscPlot_width, height = ldscPlot_height, units = "cm", limitsize = FALSE)


# +++++++++ Create Tables for supplementary: Results LD score regression analysis
LdScoreResultsExport=subset(LdScoreResults, select = c(label_clean, rg, rg_se, pval, pval_fdr))
LdScoreResultsExport$pval=formatC(LdScoreResultsExport$pval, digits=2) # Format p-values
LdScoreResultsExport$pval_fdr=formatC(LdScoreResultsExport$pval_fdr, digits=2) # Format p-values
# Format decimals
LdScoreResultsExport=formatNum(vecNum=c("rg", "rg_se"), vecP=c("pval", "pval_fdr"), df=LdScoreResultsExport)
# Change column names
colnames(LdScoreResultsExport)=c("external trait", "rg", "se (rg)", "p (rg)", "p (rg, fdr corrected)")
```

</br></br></br>

# Mendelian Randomization Analysis

[](#mendelian-randomization-analysis)

``` r
################################################
# ======== Mendelian Randomization =============
################################################

# Load MAF data
# Download reference file with MAF data: https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v/file/576598996073
MAFdat<-fread(paste0(HOME,"/data/reference.1000G.maf.0.005.txt"))
str(MAFdat)


clumpData=function(clump_input, name){
  print(paste0("Iteration: for ",name))
  print("Save SNP file")
  print(paste0("Number of included SNPs: ", NROW(clump_input)))
  system(paste0("rm ", HOME, "/data/clump/", name, ".pvalsfile.clumped" ))
  write.table(clump_input,
              file= paste0(HOME, "/data/clump/",  name, ".pvalsfile"),
              sep="\t",
              row.names = FALSE,
              col.names=T,
              quote=F)
  system(paste0(HOME, "/data/clump/plink --bfile ", HOME, "/data/clump/g1000_eur --clump ", HOME, "/data/clump/",  name, ".pvalsfile --clump-snp-field SNP --clump-field P --out ", HOME, "/data/clump/", name, ".pvalsfile --clump-kb 250 --clump-r2 0.1 --clump-p1 1"))
  print("Remove SNP file")
  system(paste0("rm ", HOME, "/data/clump/", name, ".pvalsfile" ))
  clumped=NULL
  clumped=read.table(paste0(HOME, "/data/clump/", name, ".pvalsfile.clumped"), header = TRUE)
  clumpedSNPsMerged=subset(clump_input, SNP %in% unique(  clumped$SNP))
  print(paste0("Numbr of SNPs included after clumping ", NROW(clumpedSNPsMerged) ) )
  return(clumpedSNPsMerged)
}

# Filter according to Qsnp
postGWA_clumpFiltered=imap(postGWA_clump, function(x, y) filterQtest(x,y))
# Prepare sumstats for GWA without significant hits
sumStatsListFiltered=imap(sumStatsList, function(x, y) filterQtest(x,y))
clumpedMR_filtered=list()


# Extraxt exposure SNPs
for ( i in 1:length(phenoNames) ) {
  name=names(postGWA_clumpFiltered)[i]
  clump_input=postGWA_clumpFiltered[[i]]
  print(paste0("Iteration: for ",name))
  print(paste0("Number of significant SNPs ", NROW(clump_input) ) )

  nTopsnpsSelect=10
  if(NROW(clump_input)<=nTopsnpsSelect){
    print(paste0("Less than ", nTopsnpsSelect, " GWA significant LD independent hits"))
    clumpedSNPsMerged=clumpData(sumStatsListFiltered[[i]], name=name)
    clump_output=as.data.frame(clumpedSNPsMerged %>%
                                 top_n(nTopsnpsSelect, -P))
    print(paste0(NROW(clump_output), " top SNPs included after clumping"))
    clumpedMR_filtered[[i]]=subset(clump_output, select=c(SNP,BETA_STD, SE_STD, A1, A2, P ))
  } else {
    print(paste0("Clumping of ", NROW(clump_input), " GWA significant hits" ) )
    clumpedMR_filtered[[i]]=subset(clump_input, select=c(SNP,BETA_STD, SE_STD, A1, A2, P ))
  }

}
names(clumpedMR_filtered)=phenoNames
str(clumpedMR_filtered)

# Report on SNPs used as instruments
mrSNPdescription=function(df){
  P_rangeL=range(df$P, na.rm = FALSE)[1]
  P_rangeU=range(df$P, na.rm = FALSE)[2]

  nSNPsIncluded=NROW(df)
  nSNPsSig=NROW(subset(df, P < 5e-8))
  nSNPsNs=NROW(subset(df, P >= 5e-8))

  dfOut=data.frame(P_rangeL,
                   P_rangeU,
                   nSNPsIncluded,
                   nSNPsSig,
                   nSNPsNs)
  return(dfOut)
}

mrSNPdescriptionList=lapply(clumpedMR_filtered, function(x) mrSNPdescription(x))



# Add column containing Minor Allele Frequency
prepareExposure=function(df, name, SNPlabel="SNP"){
  # Extract N per GWAS
  df_comb=merge(df, MAFdat, by="SNP", all.x=T, suffixes=c("", ".mafDF") )

  print(paste0("Prepare exposure data for ",  name))
  df=df_comb
  if(NROW(df_comb)==0){
    print("No hits included")
    df_out_format=df_comb
  } else{
    df_out = data.frame(SNP = df$SNP,
                        beta = df$BETA_STD,
                        se = df$SE_STD,
                        effect_allele = df$A1,
                        other_allele =df$A2,
                        pval = df$P,
                        exposure = name,
                        eaf = df$MAF)

    df_out_format=format_data(df_out, type="exposure")
  }
  df_out_format$exposure=name
  return(df_out_format)
}


# Prepare exposure format
postGWA_clump_Qfilter_MAF_exp=imap(clumpedMR_filtered, function(x, y) prepareExposure(x,y))
names(postGWA_clump_Qfilter_MAF_exp)=phenoNames
str(postGWA_clump_Qfilter_MAF_exp)



# ===== Harmonize exposure data with outcome data and run MR

# Create lists for output
list_out=list()
MRout_out=list()

# Run MR in loop
for ( j in 1:length(phenoNames) ) {
  df_exp=postGWA_clump_Qfilter_MAF_exp[[j]]
  name=phenoNames[j]

  print(paste0("Start MR for ",name))

  for ( i in 1:length(phenoNames) ) {
    print(paste0("Read in ", phenoNames[i] ))
    df_out=sumStatsList[[i]]
    label_out=phenoNames[i]

    # Select instruments
    exposureSNP=df_exp$SNP
    df_out = subset(df_out, df_out$SNP %in% exposureSNP)
    df_out=merge(df_out, MAFdat, by="SNP", all.x=T, suffixes=c("", ".mafDF") )

    # Format data frame
    outcome = data.frame(SNP = df_out$SNP,
                         beta = df_out$BETA_STD,
                         se = df_out$SE_STD,
                         effect_allele = df_out$A1,
                         other_allele =df_out$A2,
                         pval = df_out$P,
                         eaf = df_out$MAF)

    # Convert to MR format
    outcome_df <- format_data(outcome, type="outcome")
    outcome_df$outcome=label_out

    # Harmonize
    dat <- harmonise_data(df_exp, outcome_df)
    # Perform MR
    MRout = mr(dat,  method_list=c( "mr_ivw"))

    # mr_egger_regression
    MRegger = mr(dat,  method_list=c( "mr_egger_regression")) # "mr_egger_regression" "mr_weighted_median", "mr_weighted_mode"
    colnames(MRegger)=paste0(colnames(MRegger), "_egger")
    MReggerSub=subset(MRegger,  select=c(b_egger, se_egger, pval_egger))

    # Check for pleiotropy
    pleiotropy=mr_pleiotropy_test(dat)
    pleiotropySub=subset(pleiotropy,  select=c(egger_intercept, pval))
    colnames(pleiotropySub)=c("intercept_egger" ,"intercept_pval_egger")
    # Heterogeneity test
    # Heterogeneity in causal effects amongst instruments is an indicator of potential violations of IV assumptions (Bowden et al., 2017a). Heterogeneity can be calculated for the IVW and Egger estimates, and this can be used to navigate between models of horizontal pleiotropy (Bowden et al., 2017a).
    HETmr=mr_heterogeneity(dat, method_list=c( "mr_ivw"))
    HETmrSub=subset(HETmr, select=c(Q_pval))

    if(NROW(pleiotropySub)==0){
      pleiotropySub[nrow(pleiotropySub)+1,] <- NA
    }
    if(NROW(MReggerSub)==0){
      MReggerSub[nrow(MReggerSub)+1,] <- NA
    }
    if(NROW(HETmrSub)==0){
      HETmrSub[nrow(HETmrSub)+1,] <- NA
    }
    list_out[[i]]=cbind(MRout, MReggerSub, pleiotropySub, HETmrSub)

  }

  names(list_out)=phenoNames
  print(paste0("Finished MR for ",name))
  outcomeDF=do.call(rbind, list_out)
  MRout_out[[j]]=outcomeDF

}

# Bind lists
MRout_bind=do.call(rbind, MRout_out)
# Remove estimates where exposure and outcome are the same phenotype
MRout_bind=subset(MRout_bind, MRout_bind$outcome!=MRout_bind$exposure)
MRout_bind=subset(MRout_bind, MRout_bind$outcome!=MRout_bind$exposure)

# Derive confidence interval
MRout_bind$ci.U=MRout_bind$b + 1.96 * MRout_bind$se
MRout_bind$ci.L=MRout_bind$b - 1.96 * MRout_bind$se

# ===== Plot MR results
MRout_plot=subset(MRout_bind, exposure %in% phenoNames & outcome %in% phenoNames)
# Remove estimates from MR including only a single SNP
MRout_plot$b[MRout_plot$nsnp==1]=NA
MRout_plot$ci.U[MRout_plot$nsnp==1]=NA
MRout_plot$bci.L[MRout_plot$nsnp==1]=NA
# Create label indicating direct causation verus reverse causation
MRout_exposure=subset(MRout_plot, MRout_plot$exposure=="commonLiability" & MRout_plot$outcome!="commonLiability")
MRout_exposure$model="exposure"
MRout_exposure$label=MRout_exposure$outcome
MRout_outcome=subset(MRout_plot, MRout_plot$outcome=="commonLiability" & MRout_plot$exposure!="commonLiability")
MRout_outcome$model="outcome"
MRout_outcome$label=MRout_outcome$exposure

# Combine
MRout_exposure_outcome=rbind(MRout_exposure, MRout_outcome)

# Relabel
MRout_exposure_outcome$label=recodeName(MRout_exposure_outcome$label)
MRout_exposure_outcome$model=revalue(MRout_exposure_outcome$model, c("outcome"="reverse causation"))
MRout_exposure_outcome$model=revalue(MRout_exposure_outcome$model, c("exposure"="direct causation"))

MRout_exposure_outcome$loadings=ifelse(grepl("depen", MRout_exposure_outcome$label)==TRUE, 0.79, 0.39)
head(MRout_exposure_outcome)
MRout_exposure_sel=subset(MRout_exposure_outcome, exposure=="commonLiability")

MRList=list(MRout_exposure_sel, MRout_exposure_outcome)
# Generate pLot
MR_plotList=list()
for ( i in 1:length(MRList) ) {
  MR_plotList[[i]] <- ggplot(MRList[[i]], aes(col=model)) +
    scale_color_manual(values = c("navy", "grey")) +
    geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
    geom_linerange(aes(x = label, ymin = ci.L, ymax = ci.U),
                   lwd = 1, position = position_dodge(width = 1/2)) +
    geom_pointrange(aes(x = label, y = b, ymin = ci.L,
                        ymax = ci.U),
                    lwd = 1/2, position = position_dodge(width = 1/2)) +
    coord_flip() +
    theme_classic() +
    theme(legend.position="none") +
    labs(title = "", x = "", y = "", color = c(""))  +   scale_shape_discrete(name  ="") +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10)) +
    scale_y_continuous(expression((italic(beta)[std])))

  if((i==1)==TRUE) {
    print("add dots")
    MR_plotList[[i]] =   MR_plotList[[i]] +
    geom_point(data=MRout_exposure_outcome, aes(y=loadings,x=label), colour="darkred", size = 3,   position = position_dodge(width = 1))
    }

  print(  MR_plotList[[i]])
}

MR_plot = MR_plotList[[1]]
MR_plot_supp = MR_plotList[[2]] + theme(legend.position="top")


# Save plot
MrPlot_width=12
MrPlot_height=15
ggsave(paste0(HOME,"/results/figures/MrPlot.pdf"), MR_plot,  width = MrPlot_width, height = MrPlot_height, units = "cm", limitsize = FALSE)
ggsave(paste0(HOME,"/results/figures/MrPlot.svg"), MR_plot,  width = MrPlot_width, height = MrPlot_height, units = "cm", limitsize = FALSE, bg="white")


# Save plot (supplement)
ggsave(paste0(HOME,"/results/figures/MrPlot_supp.pdf"), MR_plot_supp,  width = MrPlot_width, height = MrPlot_height, units = "cm", limitsize = FALSE)

# Create Table for supplement
# Add SNP stats
mrSNPdescriptionList_DF=mergeLongFormat(mrSNPdescriptionList, names(mrSNPdescriptionList), cleanLabel = FALSE)
# Format decimals
mrSNPdescriptionList_DF=formatNum(vecP=c("P_rangeL", "P_rangeU"), df=mrSNPdescriptionList_DF)
mrSNPdescriptionList_DF$description=paste0( mrSNPdescriptionList_DF$nSNPsIncluded, " genetic variants were selected as instruments, of which ", mrSNPdescriptionList_DF$nSNPsSig, " variants were genome-wide significant. The p-values of the selected variants range from p=", mrSNPdescriptionList_DF$P_rangeL, " to ", mrSNPdescriptionList_DF$P_rangeU)
MRout_export = merge(MRout_bind, mrSNPdescriptionList_DF, by.x="exposure", by.y="Label", all.x=T)

# Get vector indicating the number of instruments used to index the common liability
nInsturmentsComLiab=subset(MRout_export, exposure=="commonLiability")$nsnp[1]
# Recode names
MRout_export$exposure=recodeName(MRout_export$exposure)
MRout_export$outcome=recodeName(MRout_export$outcome)
# Format decimals
MRout_export=formatNum(vecNum=c("b", "se", "ci.L", "ci.U"),vecP=c("pval", "intercept_pval_egger"), df=MRout_export)
# Create label summarizing beta (95%CI)
MRout_export$est_CI=paste0(MRout_export$b, " (", MRout_export$ci.L, "; ", MRout_export$ci.U, ")" )
# Select subset
MRout_export=subset(MRout_export, select=c(exposure, outcome, est_CI, pval, intercept_pval_egger , nsnp, description))
# Rename column names
colnames(MRout_export)=c("Exposure", "Outcome", "beta (95% CI) (IVW)", "p-value (IVW)","p-value (MR-Egger intercept)","number of genetic instruments used in MR", "description of SNPs selected as instruments")
```

</br></br></br>

# Create Supplementary tables

[](#create-supplementary-tables)

``` r
################################################################
# ====================== Create output tables ==================
################################################################

ColStart=2
RowHeader=2
RowSubheaderStart=3
RowSubheaderEnds=6
RowTable=7

# Create info text
createInfo=function(dataInfoPath){
  datOut=read.csv(dataInfoPath,header=T)
  datOut$X=NULL
  datOut_merged=paste0(datOut[,1],": " ,datOut[,2])
  return(datOut_merged)
}

# define style
hs1 <- createStyle(halign = "CENTER", textDecoration = "Bold",
                   border = "Bottom", fontColour = "black", fgFill = "white")

h_info <- createStyle(halign = "left", textDecoration = "Bold",
                      border = "Bottom", fontColour = "black", fgFill = "white")

addTable=function(sheet, table){
  writeDataTable(wb, sheet, table, headerStyle=hs1, tableStyle = "TableStyleLight1",
                 startRow = RowTable, startCol = ColStart)
  setColWidths(wb, sheet, cols=2:10, widths = 15)
}

# HEADER
headerFunc=function(TITLE, sheet){
  writeData(wb, sheet = sheet, TITLE,
            colNames = FALSE, rowNames = FALSE,
            startCol = ColStart, startRow = RowHeader)
}
# INFO ROW
InfoFunc=function(TITLE, sheet){
  writeData(wb, sheet = sheet, TITLE,
            colNames = FALSE, rowNames = FALSE,
            startCol = ColStart, startRow = RowSubheaderStart)
}

# Create new workbook
setwd(paste0(HOME,"/results/tables/"))
wb <- openxlsx::createWorkbook()




# ================================ TABLE 1 ================================ #
# +++++++++  Read in included GWAS summary statistics
gwasSumStastAll=openxlsx::read.xlsx(paste0(HOME,"/analysis/GenSEM_GWAS.xlsx"))
gwasSumStast=subset(gwasSumStastAll, directory!="exclude")
# Add row with common liability results
dfCommon=as.data.frame(matrix(ncol=length(colnames(gwasSumStastAll)), nrow=1))
colnames(dfCommon)=colnames(gwasSumStastAll)
dfCommon$label="commonLiability"
dfCommon$label_clean = "Common liability"
dfCommon$N=sumStatsList[[1]]$N[1]
gwasSumStastRbind=rbind(gwasSumStast, dfCommon)

# Add heritability estimates
ldH2=readRDS(paste0(HOME, "/output/rds/h2outDF.rds"))
gwasSumStast=merge(gwasSumStastRbind, ldH2, by.x="label", by.y="trait", all.x=T)
# Format for supplementary table
gwasSumStast=subset(gwasSumStast, select=c(label_clean, N, h2, h2_se, intercept, link, Publication))
colnames(gwasSumStast)=c("Summary statistic file", "Sample size", "SNP-based heritability (h2 estimate)", "SNP-based heritability (h2 standard error)", "Intercept", "Link to summary statistics file", "Link to study")

# Create new sheet
LDscoreRegressionData="sTable 1"
addWorksheet(wb, LDscoreRegressionData)
# Add datatable
title_name=paste0("sTable 1. Overview of summary statistic files used in the multivariate genome-wide association study")
infoTextGWAoverview=paste0("List of summary statistic files used to derive the common hertiable liability to addiction, as well as all summary statistic files used in LD score regression analysis assessing the genetic correlations between the common liability and other traits. Heritability (h2) estimates and the intercepts were estimated using univariate LD score regression implemented in GenomicSEM")
Info_text=infoTextGWAoverview
sheet=LDscoreRegressionData
table=gwasSumStast
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT1=list(data=gwasSumStast,
     description=infoTextGWAoverview,
     title=title_name)


# ================================ TABLE 2 ================================ #
# Read in data
modelName="cancigalc"
GenomicSEM_ldscAll=readRDS(paste0(HOME,"/output/rds/LDSCoutput_", modelName, "_SUD.rds"))
data.frame(pheno=colnames(GenomicSEM_ldscAll$S),
           intercept = diag(GenomicSEM_ldscAll$I))

# Get mean correlations
GenomicSEM_ldsc_cor=GenomicSEM_ldscAll$S_Stand
GenomicSEM_ldsc_cor=ifelse(as.character(GenomicSEM_ldsc_cor)=="1", NA, as.numeric(GenomicSEM_ldsc_cor))
# Format decimals
GenomicSEM_ldscAll$S_Stand=round(GenomicSEM_ldscAll$S_Stand, 3)

# Replace cor=1 with h2 estimates
for (i in 1:length(colnames(GenomicSEM_ldscAll$S))) {
  GenomicSEM_ldscAll$S_Stand[i,i]=round(GenomicSEM_ldscAll$S[i,i],3)
}

colnames(GenomicSEM_ldscAll$S_Stand)=recodeName(colnames(GenomicSEM_ldscAll$S_Stand))
rownames(GenomicSEM_ldscAll$S_Stand)=colnames(GenomicSEM_ldscAll$S_Stand)

# Extract h2:
# => SNP-based heritability on the diagonal and genetic correlations off-diagonal
# => S consists of heritabilities on the diagonal and genetic covariances (co-heritabilities) on the off-diagonal

# Format table
GenomicSEM_ldscTable=cbind(Phenotype = colnames(GenomicSEM_ldscAll$S_Stand), GenomicSEM_ldscAll$S_Stand)

# Add new sheet
GeneticCorrelation_supplement="sTable 2"
addWorksheet(wb, GeneticCorrelation_supplement)
# Add parameters
title_name="sTable 2. Genetic correlations between the individual substance use phenotypes"
sheet=GeneticCorrelation_supplement
table=as.data.frame(GenomicSEM_ldscTable)
Info_text=paste0("Shown are the genetic correlations between each of the cigarette, alcohol and cannabis use phenotypes, with heritability estimates displayed down the diagonal. The mean genetic correlation is rg=", round(mean(GenomicSEM_ldsc_cor, na.rm=TRUE),2), " [sd=", round(sd(GenomicSEM_ldsc_cor, na.rm=TRUE),2), ", median=", round(median(GenomicSEM_ldsc_cor, na.rm=TRUE),2), " and range (", round(range(GenomicSEM_ldsc_cor, na.rm=TRUE),2)[1], "-", round(range(GenomicSEM_ldsc_cor, na.rm=TRUE),2)[2], ")]")

# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT2=list(data=as.data.frame(GenomicSEM_ldscTable),
            description=Info_text,
            title=title_name)


# ================================ TABLE 3 ================================ #
# Import correlation results from ldsc regression analysis
CommonFac=readRDS(paste0(HOME,"/output/rds/CommonFac_model.rds"))

# ++ Extract standardized estimates +++
CommonFacRes=subset(as.data.frame(CommonFac$results), select=c(lhs,rhs ,STD_All ), lhs=="F1" & rhs!= "F1" )
CommonFacRes$STD_All=round(CommonFacRes$STD_All,2)

# Extract mean fit indicies
loadingsSum=subset(CommonFacRes, lhs != rhs)
loadings_supplement_info=paste0("The mean loading of the indicators is m=", round(mean(loadingsSum$STD_All, na.rm=TRUE),2), " [sd=", round(sd(loadingsSum$STD_All, na.rm=TRUE),2), " , range (", round(range(loadingsSum$STD_All, na.rm=TRUE),2)[1], "-", round(range(loadingsSum$STD_All, na.rm=TRUE),2)[2], ")]")

# Derive shared non-sphared variance estimates
loadingsSum$var_shared=round(loadingsSum$STD_All^2,2)
loadingsSum$non_var_shared=round(1-loadingsSum$STD_All^2,2)
explained_variance_factor=round(mean(loadingsSum$STD_All^2, na.rm=TRUE),4)*100
explained_variance_sd=round(sd(loadingsSum$STD_All^2, na.rm=TRUE),4)*100
explained_variance_range=round(range(loadingsSum$STD_All^2, na.rm=TRUE),4)*100

# Select columns for table
loadingsSum$lhs=recodeName(loadingsSum$lhs)
loadingsSum$rhs=recodeName(loadingsSum$rhs)
loadingsSum$Indicator=ifelse(loadingsSum$lhs=="Common liability", paste0(loadingsSum$lhs, " ~ ", loadingsSum$rhs), paste0(loadingsSum$lhs, " ~~ ", loadingsSum$rhs))
CommonFacRes_selected=subset(loadingsSum, select = c(Indicator, STD_All, var_shared, non_var_shared ))
colnames(CommonFacRes_selected)=c("Indicator", "Estimate (standardized)", "Variance (common liability)", "Variance (specific)")

# Create new sheet
modelEstimates="sTable 3"
addWorksheet(wb, modelEstimates)
modelEstimates_info=paste0("All substance use phenotypes are scaled so that higher scores indicate more problematic substance use, such as higher frequency of use and the presence of substance use dependence. Estimate (standardized) = the standardized linear relationship between the factor and each of the individual substance use phenotypes. Variance (common liability) = variance per phenotype explained by the common liability. Variance (specific) = variance per phenotype not explained by the common liability. On average, the common liability factor accounted for ", explained_variance_factor, "% (range ", explained_variance_range[1], "%-", explained_variance_range[2], "%) of the genetic variance in the six substance use phenotypes.")
# Add parameters
title_name="sTable 3. Estimates of the genomic factor model"
sheet=modelEstimates
table=CommonFacRes_selected
Info_text=modelEstimates_info
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT3=list(data=as.data.frame(CommonFacRes_selected),
            description=modelEstimates_info,
            title=title_name)



# ================================ TABLE 4 ================================ #
# Format decimals
CommonFac$modelfit=formatNum(vecNum=c("chisq", "CFI", "SRMR"), vecP=c("p_chisq"),df=CommonFac$modelfit)
CommonFac$modelfit$AIC=NULL

# Create new sheet
modelFit="sTable 4"
addWorksheet(wb, modelFit)
modelFit_info="chisq = chi-square statistic; CFI = the comparative fit index; SRMR =  the standardized root mean square residual"
# Add parameters
title_name="sTable 4. Fit indicies of the factor model representing the common heritable liability to addiction"
sheet=modelFit
table=CommonFac$modelfit
Info_text=modelFit_info
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)


# Store as list
listT4=list(data=as.data.frame(CommonFac$modelfit),
            description=modelFit_info,
            title=title_name)


# ================================ TABLE 5 ================================ #
GWA_shortSum$GWAS=recodeName(phenoNames)
# Create new sheet
GWA_short="sTable 5"
addWorksheet(wb, GWA_short)
GWA_short_info=paste0("SNP = single nucleotide polymorphism. The effective sample size of the common liability genome-wide association study (GWAS) was calculated using the formula described in the sMethods (Supplement). SNPs (shared) represent SNPs that operate via the common liability (i.e., with Qsnp p>5×10−8) versus SNPs (non-shared) that show heterogeneous effects across the individual cigarette, alcohol and cannabis use phenotyes (i.e., Qsnp p<5×10−8)")
# Add parameters
title_name="sTable 5. Overview of the results from the multi- and univariate genome-wide association analyses"
sheet=GWA_short
table=GWA_shortSum
Info_text=GWA_short_info
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT5=list(data=as.data.frame(GWA_shortSum),
            description=GWA_short_info,
            title=title_name)


# ================================ TABLE 6 ================================ #
# Format as Table
postGWA_clumpOrdered=lapply(postGWA_clump, function(x) arrange(x, P)) # Arrange according to p-value
clumpResTable=mergeLongFormat(postGWA_clumpOrdered, recodeName(phenoNames), cleanLabel = TRUE)

# Format decimals
clumpResTable=formatNum(vecP=c("Q_chisq_pval"," P"), vecNum=c("BETA", "SE", "P_logp"), df=clumpResTable)
clumpResTable=subset(clumpResTable, select=c(Label, ID_GENE, BP, CHR, A1, A2, BETA, SE, P, Q_chisq_pval, chi_dich, consequence, description))
colnames(clumpResTable)=c("GWA data","SNP (annotated gene)", "Position", "Chromosome", "A1", "A2", "BETA", "SE", "P (GWA)", "P (chi-square)", "P (chi-square, dichotomized)","consequence", "description")

# Create new sheet
GWA_clump="sTable 6"
addWorksheet(wb, GWA_clump)
GWA_clump_info="Shown are only lead single nucleotide polymorphisms (SNPs), defined as LD-independent genome-wide significant variants. The columns 'P (chi-square, dichotomized)' indicates whether the genetic variant is likely to operate via the common liability (Q_ns = Qsnp p>5×10−8) or shows heterogeneous effects across the individual substance use phenotyes (Q_sig = Qsnp p<5×10−8)"
# Add parameters
title_name="sTable 6. Summary of lead genetic variants associated with the commmon liability to addiction"
sheet=GWA_clump
table=clumpResTable
Info_text=GWA_clump_info
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT6=list(data=as.data.frame(clumpResTable),
            description=GWA_clump_info,
            title=title_name)

# ================================ TABLE 7 ================================ #
postGWA_clump_genesDF$Label=recodeName(postGWA_clump_genesDF$Label)
postGWA_clump_genes_select=subset(postGWA_clump_genesDF, select = c(SNP, BETA, P, Q_chisq_pval, ID_GENE, Label))

# Create new sheet
GWA_clump_long="sTable 7"
addWorksheet(wb, GWA_clump_long)
GWA_clump_long_info="Table of the results displayed in Figure 2 of the main manuscript"
# Add parameters
title_name="sTable 7. Associations of genetic variants with the common liability to addiction and the individual substance use phenotypes"
sheet=GWA_clump_long
table=postGWA_clump_genes_select
Info_text=GWA_clump_long_info
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT7=list(data=as.data.frame(postGWA_clump_genes_select),
            description=GWA_clump_long_info,
            title=title_name)


# ================================ TABLE 8 ================================ #
# Remove SNP based on Qsnp
mergeQTLTableList=imap(postGWA_pos_eqtl, function(x, y) filterQtest(x,y))
mergeQTLTableListOrdered=lapply(mergeQTLTableList, function(x) arrange(x, P_GWAS))
# Remove lists with empty rows
mergeQTLTableClean=mergeQTLTableListOrdered[sapply(mergeQTLTableListOrdered, nrow)>0]
# Format dataframe
mergeQTLTableDF=mergeLongFormat(mergeQTLTableClean,
                                recodeName(names(mergeQTLTableClean)), cleanLabel = TRUE)
# Add gene description
enseQTL_select =gconvert(query = na.omit(mergeQTLTableDF$ens_eqtl), organism = "hsapiens",
                         target="ENSG", mthreshold = Inf, filter_na = TRUE)[c("input", "description")]
enseQTL_selectClean=as.data.frame(subset(enseQTL_select, input!="None", select=c(input, description))) %>% distinct(input, .keep_all= TRUE)
mergeQTLTableMappedDF=merge(mergeQTLTableDF, enseQTL_selectClean, by.x="ens_eqtl", by.y="input", all.x=T, sort = F)

# Select columns
mergeQTLTable=subset(mergeQTLTableMappedDF, select=c(SNP, BP_pos, gene_pos, gene_eqtl, distance, description, Label))
colnames(mergeQTLTable)=c("SNP", "Position", "Gene name (positional mapping)","Gene name (eQTL mapping)", "Distance (in kilobase)", "description", "GWA data")
# is_best:        if it has the overall lowest P-value in a specific tissue of a QTL study
# n_qtls:         the number of putative QTLs for each gene in a specific tissue of a QTL study

# Create new sheet
eQTLData="sTable 8"
addWorksheet(wb, eQTLData)
# Add datatable
title_name="sTable 8. Summary of eQTL annotation of lead SNPs associated with the common liability and the individual substance use phenotypes"
Info_text=paste0("eQTL = expression quantitative trait loci; SNP = single nucleotide polymorphism. Position = variant position based on GRCh37. Distance = the distance of the eQTL SNP from the corresponding eQTL gene in kilobases (according to GRCh37). eQTL mapping was done for lead SNPs associated with the common liability to addiction and the individual substance use phenotypes, defined as LD-independent genome-wide significant variants. SNPs used for eQTL mapping were selected based on the Qsnp statistic: From the common liability GWA, only SNPs acting via the common liability were selected (Qsnp p>5×10−8); for the individual substance use phenotypes, only SNPs with heterogeneous effects across the individual substance use phenotypes were selected (Qsnp p<5×10−8). eQTL mapping was done using Qtilizer. Mapped to genes were only variants acting as cis-eQTL (variants in a +/- 1 megabase window around the transcription start site of a given gene) and that remained significant (p<0.05) after false discovery rate (FDR)/family-wise error rate (FWER) correction as implemented in Qtlizer.")
sheet=eQTLData
table=mergeQTLTable
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT8=list(data=as.data.frame(mergeQTLTable),
            description=Info_text,
            title=title_name)


# ================================ TABLE 9 ================================ #
# Format decimals
postGWA_clumpPhenoDF_suppl=formatNum(vecNum=c("BETA"), vecP=c("P", "P_pheno"), df=postGWA_clumpPhenoDF)
postGWA_clumpPhenoDF_suppl = subset(postGWA_clumpPhenoDF_suppl, select=c(Label, ID_GENE,  P, trait, P_pheno))
colnames(postGWA_clumpPhenoDF_suppl)=c("GWA data","SNP (mapped gene)", "p-value (current GWA)", "Trait (phenoscanner)", "p-value (phenoscanner)")
head(postGWA_clumpPhenoDF_suppl)
# Create new sheet
PhenoScanData="sTable 9"
addWorksheet(wb, PhenoScanData)
# Add datatable
title_name="sTable 9. Summary of PhenoScanner results for lead SNPs"
Info_text=paste0("Phenotypic associations for lead single nucleotide polymorphisms (SNPs) associated with each of the GWA analysis, defined as the top n=",topNGWA, " LD-independent SNPs associated with the common liability to addiction and the individual substance use phenotypes. SNPs included in the phenoscanner search were selected based on the Qsnp statistic: From the common liability GWA, only SNPs acting via the common liability were selected (Qsnp p<5×10−8); for the individual substance use phenotypes, only SNPs with heterogeneous effects across the individual substance use phenotypes were selected (Qsnp p<5×10−8). Listed are the ",nTopSNPsPheno, " most significant genotype-phenotype associations from the PhenoScanner database")
sheet=PhenoScanData
table=postGWA_clumpPhenoDF_suppl
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT9=list(data=as.data.frame(postGWA_clumpPhenoDF_suppl),
            description=Info_text,
            title=title_name)

# ================================ TABLE 10 ================================ #
depictOutTable=mergeLongFormat(depictOut, phenoNames, cleanLabel = TRUE)
depictOutTable$p_corr=NULL
depictOutTable$p_corr_dich=NULL
depictOutTable$Label=recodeName(depictOutTable$Label)

# Create new sheet
DEPICTData="sTable 10"
addWorksheet(wb, DEPICTData)
# Add datatable
title_name="sTable 10. Results obtained from DEPICT tissue/cell type enrichment analysis"
Info_text="Tissue/cell type enrichment analysis was conducted for the common liability GWA (including only SNPs with Qsnp p>5×10−8). For all summary statistics analysed in DEPICT, included were a set of LD-independent SNPs (r2<0.05 within 500 kb) outside genome-wide significance (p<1×10-5)."
sheet=DEPICTData
table=depictOutTable
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT10=list(data=as.data.frame(depictOutTable),
            description=Info_text,
            title=title_name)



# ================================ TABLE 11 ================================ #
# Create Table for supplement
PascalSuppAllDF=listPascalSelected_DF_ordered
PascalSuppAllDF$label_fac= orderLabels(PascalSuppAllDF$label)
PascalSuppAllDF=PascalSuppAllDF[with(PascalSuppAllDF, order(label_fac, chi2Pvalue)),]
PascalSuppAllDF$label_rec=recodeName(PascalSuppAllDF$label)
PascalSuppSel=subset(PascalSuppAllDF, select=c(label_rec, Name, pathAnno, chi2Pvalue, p_corr, nTests, PathData))
colnames(PascalSuppSel)=c("GWA data", "Pathway name","Pathway name (annotated)" , "p-value", "p-value (FDR corrected)", "Number of tests per gene-set", "Gene-set input")

# Create new sheet
PascalSupp="sTable 11"
addWorksheet(wb, PascalSupp)
# Add datatable
title_name="sTable 11. Results from Pascal pathway analysis on the common liability and substance use GWAs"
Info_text=paste0("From the common liability GWA results, included were SNPs that likely operated through the common liability (filtered according to Qsnp p>5×10−8). Listed in the table are all the pathways that were significant after FDR correction for at least one of the included phenotypes, including the common liability to addiction and the individual substance use phenotypes")
sheet=PascalSupp
table=PascalSuppSel
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT11=list(data=as.data.frame(PascalSuppSel),
             description=Info_text,
             title=title_name)

# ================================ TABLE 12 ================================ #
# Create new sheet
LDscoreRegression="sTable 12"
addWorksheet(wb, LDscoreRegression)
# Add datatable
title_name="sTable 12. Genetic correlations between the common liability to addiction and other traits"
Info_text=""
sheet=LDscoreRegression
table=LdScoreResultsExport
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT12=list(data=as.data.frame(LdScoreResultsExport),
             description=Info_text,
             title=title_name)


# ================================ TABLE 13 ================================ #
# Create new sheet
MR_res="sTable 13"
addWorksheet(wb, MR_res)
# Add datatable
title_name="sTable 13. Results from Mendelian Randomization analysis"
Info_text=paste0("Reported are the standardized beta coefficients obtained from Mendelian Randomization (MR) analysis assessing the bi-directional relations between the common heritable liability to addiction and the six substance use phenotypes. Direct causation was estimated as the effects of the common liability on the substance use phenotypes, including n=", nInsturmentsComLiab, " genome-wide significant genetic variants (p<5×10-8) operating through the common liability (Qsnp p>5×10-8) as instruments. Reverse causation was assessed by estimating the effects of the individual substance use phenotypes on the common liability, including only genome-wide significant genetic variants that were unlikely to operate through the common liability (Qsnp p<5×10-8) as instruments. In instances where the GWA data did not contain genetic variants reaching genome-wide significance, we selected the top ", nTopsnpsSelect, " LD-independent SNPs as instruments instead")
sheet=MR_res
table=MRout_export
# Run functions
addTable(sheet, table)
headerFunc(title_name, sheet)
InfoFunc(Info_text, sheet)

# Store as list
listT13=list(data=as.data.frame(MRout_export),
             description=Info_text,
             title=title_name)




################################################################
# ====================== Export table ==================
################################################################

listAll=list(listT1, listT2, listT3, listT4, listT5, listT6, listT7, listT8, listT9, listT10, listT11, listT12, listT13)
saveRDS(listAll, paste0(HOME,"/results/tables/TableSum.rds"))


# Create new styles
s <- createStyle(fgFill = "#FFFFFF")
h_info <- createStyle(halign = "left",
                      border = "BOTTOM", fontColour = "black", fgFill = "white", fontSize=16, textDecoration = "Bold", numFmt="TEXT", borderColour = "black")
info_info <- createStyle(halign = "left",
                         border = NULL, fontColour = "black", fgFill = "white", fontSize=14, textDecoration = NULL, numFmt="TEXT", wrapText=TRUE)
# Run loop
for(curr_sheet in names(wb)){
  addStyle(wb,sheet = curr_sheet, s, cols=1:40, rows=1:2000, gridExpand = TRUE)
  setColWidths(wb, sheet = curr_sheet, cols=1:40, widths = 20)
  addStyle(wb,sheet = curr_sheet, h_info, cols=ColStart:20, rows=RowHeader, gridExpand = TRUE)
  addStyle(wb,sheet = curr_sheet, info_info, cols=ColStart:5, rows=RowSubheaderStart, gridExpand = TRUE)
  mergeCells(wb,sheet = curr_sheet, cols = 2:8, rows = RowSubheaderStart:RowSubheaderEnds)
}

library( openxlsx)
openxlsx::saveWorkbook(wb, paste0(HOME,"/results/tables/CommomLiabAddiction_Tables.xlsx"), overwrite = TRUE)
# Open File
openXL(wb)
```
