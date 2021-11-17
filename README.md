Genome-wide association study on the common liability to addiction
================

## Table of Contents

1.  [Description](#description)
2.  [Model without SNP effects](#model)
3.  [Model with SNP effects](#modelSNP)
4.  [Gene annotation](#functional)
5.  [Gene-set enrichment](#enrichment)
6.  [LD score regression](#ldsc)

# [Description](#description)

</br></br>

The following documentation provides a description for all analytical
steps involved in the multivariate genome-wide association study on the
common liability to addiction.

1.  Getting started with GenomicSEM

-   Read the [wiki and
    tutorial](https://github.com/GenomicSEM/GenomicSEM/wiki) for the
    R-package `GenomicSEM`

1.  Install R packages

2.  Download all required software tool, including

-   plink
-   DEPICT
-   PASCAL

</br></br>

1.  Download summary statistic files

# [Model without SNP effects](#model)

</br></br></br>

## Estimate genetic correlations

</br></br>

![](results/figures/CorrGWA.svg) </br>

To prepare the data and run all analytical steps, follow the script
published
[here](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/analysis)

</br></br>

## Specify the structural model

![](results/figures/strucModel.png)

</br></br></br></br>

# [Run the multivariate genome-wide association study](#modelSNP)

![](results/figures/ManHplot_commonLiability.jpeg)

</br></br>

# [Functional follow-up](#functional)

# [LD score regression](#ldsc)
