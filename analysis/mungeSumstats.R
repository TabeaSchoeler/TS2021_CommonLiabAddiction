#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
print(args[1])
HOME=args[1]
print(HOME)
source(paste0(HOME, "/analysis/input.R"))

# Import data from dropbox 
drop_auth(rdstoken = paste0(HOME, "/token.rds"))  # authentication for dropbox
drop_acc()
drop_download(
  local_path = paste0(HOME,"/analysis/GenSEM_GWAS.xlsx"),
  path = paste0(LOCAL,"/analysis/GenSEM_GWAS.xlsx"), overwrite=TRUE)
gwasSumStastAll=read.xlsx(paste0(HOME,"/analysis/GenSEM_GWAS.xlsx"), na.strings="NA")

# Create path name
gwasSumStastAll$DirPath=paste0(gwasDir,gwasSumStastAll$fileName)
# Get name for zipped / munged files
gwasSumStastAll$DirPathgz=paste0(HOME,"/data/processed/", gwasSumStastAll$label, ".sumstats.gz") # Include original GWAS as well as munged GWAS files
# Select only substance use traits
gwasSumStast = subset(gwasSumStastAll, directory == "gensem" )

# Get sample prevalence
# Note: for binary traits - if this reflect the sum of effective sample sizes across contributing cohorts for case/control designs, the sample prevalence should then be entered as 0.5 when running ldsc to reflect the fact that effective sample size already corrects for sample ascertainment. 
# If the input contains SNP-specific sample sizes for each row, then the munge column will only use this N if the user does not provide their own.
samplePrevalence=ifelse(is.na(gwasSumStast$N_eff)==F, 0.5, gwasSumStast$sample.prev)


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

# ======================================== Munge files in loop
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

## CHECK IF LOGIT SCALE
# If the standard error column is the standard error of the odds ratio you would set se.logit = FALSE, and if it is the standard error of a logsitic beta you would set se.logit = TRUE. Off the top of my head I'm not positive for METAL output what the scale of the standard error column is, but my guess is that it is the standard error of a logistic beta in which case as you say you would set se.logit = TRUE

# Note GENSEM
# You need to know whether the GWAS was a logistic regression or a linear regression
# Note that not all case/control studies use logistic regression. 
#     =>  e.g., binary outcomes for the HAIL GWAS), where dichotomous outcome (e.g. a case/control trait) is analyzed using a linear regression
#     =>  this is called a "linear probability model" 
#     =>  The function sumstats does know how to deal with this scenario using the linprob argument
#     =>  The package can also deal with a GWAS of a continuous trait being analyzed using linear regression (use the OLS flag in sumstats to indicate which GWAS are of continuous traits), or a case/control traits analyzed using logistic regression (the default in sumstats). The decision tree directly below can be used to decide what the correct arguments are for the sumstats function based on the scale of the GWAS outcome and how that outcome was analyzed. In the decision tree, we list the necessary arguments for each individual summary statistics file. In practice, you would combine the appropriate elements passed to each argument into a single vector to match the order the summary statistics are listed.


# On logit scale?
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
