#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
print(args[1])
HOME=args[1]
print(HOME)

modelName=args[2]
print(modelName)

source(paste0(HOME, "/analysis/input.R"))



print("Read in common liability GWA")
namesComldsc=c("commonLiability")
listUserGWA=readRDS(paste0(HOME,"/output/rds/",modelName, "_commonFacSumstats.rds"))

print("Prepare common liability GWAS data")
commonList=list(listUserGWA[["commonLiability"]],
                subset(listUserGWA[["commonLiability"]],   Q_chisq_pval >= 5e-8))

# Select relevant columns
commonList_sel=lapply(commonList, function(x) subset(x, select=c("SNP", "A1", "A2", "BETA", "SE", "P", "N")))
namesComldsc=c("commonLiability", "commonLiability_filtered")
names(commonList_sel)=namesComldsc

print("Save output on cluster")
sapply(names(commonList_sel), function(x)
  write.table(commonList_sel[[x]], file=paste0(HOME, "/output/ldsc/", x, "_", modelName, "_clean"),
              col.names=T, row.names=F, quote=F, sep="\t") )


print("Prepare summary statistic files for external traits")
# Import data from dropbox
drop_auth(rdstoken = paste0(HOME, "/token.rds"))  # authentication for dropbox
drop_acc()
drop_download(
  local_path = paste0(HOME,"/analysis/GenSEM_GWAS.xlsx"),
  path = paste0(LOCAL,"/analysis/GenSEM_GWAS.xlsx"), overwrite=TRUE)
gwasSumStastAll=read.xlsx(paste0(HOME,"/analysis/GenSEM_GWAS.xlsx"), na.strings="NA")
gwasSumStast=subset(gwasSumStastAll, directory!="exclude")

# Create path name
gwasSumStast$DirPath=paste0(gwasDir,gwasSumStast$fileName)

print("Prepare clean datasets for munging")
for ( i in 1:length(gwasSumStast$DirPath) ) {
  setwd(paste0(HOME,"/output/ldsc"))
  print(gwasSumStast$DirPath[i])
  gwa_in=fread(gwasSumStast$DirPath[i], fill=TRUE )
  colnames(gwa_in)
  gwa_select=subset(gwa_in , select=c(gwasSumStast$SNP[i],
                                      gwasSumStast$A1[i],
                                      gwasSumStast$A2[i],
                                      gwasSumStast$effect[i],
                                      gwasSumStast$se[i],
                                      gwasSumStast$P[i]))
  gwa_select$N=gwasSumStast$N[i]
  colnames(gwa_select)=c("SNP", "A1", "A2", "BETA", "SE", "P", "N")
  gwa_select$BETA=as.numeric(gwa_select$BETA)
  gwa_select$SE=as.numeric(gwa_select$SE)
  gwa_select$P=as.numeric(gwa_select$P)

  gwa_unique=gwa_select[!duplicated(gwa_select$SNP)]
  print(paste0("Remove ", NROW(gwa_select)-NROW(gwa_unique), " duplicate SNPs"))
  # Save output
  print("Save output")
  write.table(gwa_unique,
              file=paste0(HOME, "/output/ldsc/", gwasSumStast$label[i], "_clean"),
              sep="\t",
              row.names = FALSE,
              col.names = TRUE,
              quote=F)
}

system(paste0("rm ", HOME, "/output/ldsc/*.sumstats.gz"))
system(paste0("rm ", HOME, "/output/ldsc/*.log"))
system(paste0("rm ", HOME, "/output/ldsc/*.txt"))
system(paste0("rm ", HOME, "/output/ldsc/ldsc_externa*"))

print("Add common liability gwas for munging")
mungeFiles=c(paste0(gwasSumStast[["label"]], "_clean"),
             paste0(namesComldsc, "_",modelName, "_clean"))

print("Munge files")
mungeData=function(name){
  paste0("start munging for: ", name)
  setwd(paste0(HOME,"/output/ldsc"))

  munge(files = paste0(HOME, "/output/ldsc/", name),
        hm3 = paste0(HOME,"/data/processed/w_hm3.snplist"),
        trait.names=name,
        info.filter = 0.9,
        maf.filter = 0.01)
  paste0("finished munging for: ", name)
}
lapply(mungeFiles, function(x) mungeData(x))


print("Get heritability estimates and intercept")
dfldsc1=subset(gwasSumStast, select=c(label, sample.prev, population.prev))
dfldsc2=data.frame(label=paste0(namesComldsc, "_",modelName), sample.prev=NA , population.prev=NA)
ldscFiles=rbind(dfldsc1, dfldsc2)



h2outList=list()
# function for ldscore
for ( i in 1:NROW(ldscFiles) ) {

  print(paste0("Start ", ldscFiles[["label"]][i]))

  print("Estimate heritability")

  ldsch2=NULL
  ldsch2=ldsc(traits=rep(paste0(HOME,"/output/ldsc/", ldscFiles[["label"]][i], "_clean.sumstats.gz"),2),
              sample.prev= rep(ldscFiles[["sample.prev"]][i],2),
              population.prev=rep(ldscFiles[["population.prev"]][i],2) ,
              ld=paste0(HOME,"/data/processed/eur_w_ld_chr/"),
              wld=paste0(HOME,"/data/processed/eur_w_ld_chr/"),
              trait.names = rep(ldscFiles[["label"]][i],2),
              stand=T)

  print("Get standard error")
  h2_t1=ldsch2[["S"]][1]
  SE=NULL
  SE<-matrix(0, nrow(ldsch2[["S"]]), nrow(ldsch2[["S"]]))
  SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(ldsch2[["V"]]))
  h2_t1_se=SE[1]

  h2outList[[i]]=data.frame(trait=ldscFiles[["label"]][i],
                            h2=h2_t1,
                            h2_se=h2_t1_se,
                            intercept=ldsch2[["I"]][1])

}

h2outDF=as.data.frame(do.call("rbind", h2outList))
saveRDS(h2outDF, paste0(HOME, "/output/rds/h2outDF.rds"))
print("Export results to dropbox")
drop_auth(rdstoken = paste0(HOME, "/token.rds"))  # authentication for dropbox
drop_acc()
drop_upload(paste0(HOME,"/output/rds/h2outDF.rds"), path = paste0(LOCAL, "/output/rds"))




# =========== LD score regression
print("Start ld score regression analysis")
ldscCrossTrais=function(name, listLDSCallsub=list()){

  writeOut=function(line){
    write(line,
          file=paste0(HOME, "/output/ldsc/", name, "_ldscOut.txt"),
          append=TRUE)
  }

  setwd(paste0(HOME,"/output/ldsc"))

  COMMON_trait=paste0(HOME,"/output/ldsc/", name, "_clean.sumstats.gz")
  COMMON_trait_name=name

  for ( i in 1:length(gwasSumStast[["label"]]) ) {
    print(paste0("Start ", gwasSumStast[["label"]][i]))

    writeOut(paste0("Start ", gwasSumStast[["label"]][i]))

    LDSCoutputLoop=ldsc(traits=c(COMMON_trait,
                                 paste0(HOME,"/output/ldsc/", gwasSumStast[["label"]][i], "_clean.sumstats.gz")),
                        sample.prev= c(NA, gwasSumStast[["sample.prev"]][i]),
                        population.prev=c(NA,gwasSumStast[["population.prev"]][i]),
                        ld=paste0(HOME,"/data/processed/eur_w_ld_chr/"),
                        wld=paste0(HOME,"/data/processed/eur_w_ld_chr/"),
                        trait.names = c(COMMON_trait_name, gwasSumStast[["label"]][i]),
                        stand=T)

    print("Get standard errors of ldscore regression")
    k<-nrow(LDSCoutputLoop[["S_Stand"]])
    SE<-matrix(0, k, k)
    SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutputLoop[["V_Stand"]]  ))
    rg_se=SE[2,1]
    print(paste0("Finished ", gwasSumStast[["label"]][i]))

    listLDSCallsub[[i]]=data.frame(commonGWAS=name,
                                   trait=gwasSumStast[["label"]][i],
                                   rg=LDSCoutputLoop[["S_Stand"]][1,2],
                                   rg_se=rg_se)

    writeOut(paste0("Finished ", gwasSumStast[["label"]][i]))

  }
  return(listLDSCallsub)
}

int <- 2
cl <- makeCluster(int)
clusterExport(cl, varlist=c("ldscCrossTrais", "R_libPaths", "HOME", "gwasSumStast"))
clusterEvalQ(cl, R_libPaths )
clusterEvalQ(cl, .libPaths(R_libPaths) )
clusterEvalQ(cl, library(GenomicSEM) )
clusterEvalQ(cl, c("ldscCrossTrais", "R_libPaths", "HOME", "gwasSumStast", library(GenomicSEM) ))

# Run in loop for Qsnp filtered and unfiltered GWAS
listLDname=paste0(namesComldsc, "_",modelName)
ldscExternalOut=parLapply(cl, listLDname, function(x) ldscCrossTrais(name=x, listLDSCallsub=list()))
ldscExternalOutDF1=as.data.frame(do.call("rbind", ldscExternalOut[[1]]))
ldscExternalOutDF2=as.data.frame(do.call("rbind", ldscExternalOut[[2]]))
ldscExternalOutDF=as.data.frame(rbind(ldscExternalOutDF1, ldscExternalOutDF2))
print("Save results on cluster")
saveRDS(ldscExternalOutDF, paste0(HOME, "/output/rds/ldscExternalOutDF.rds"))

print("Export results to dropbox")
drop_auth(rdstoken = paste0(HOME, "/token.rds"))  # authentication for dropbox
drop_acc()
drop_upload(paste0(HOME,"/output/rds/ldscExternalOutDF.rds"), path = paste0(LOCAL, "/output/rds"))




