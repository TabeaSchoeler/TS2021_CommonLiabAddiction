Genome-wide association study on the common liability to addiction
================

</br></br></br></br></br>&lt;

## Table of Contents

1.  [Description](#description)
2.  [Model without SNP effects](#model)
3.  [Model with SNP effects](#modelSNP)
4.  [Gene annotation](#functional)
5.  [Gene-set enrichment](#enrichment)
6.  [LD score regression](#ldsc)

</br></br>

#### *The documentation below provides a brief overview of all steps involved in the analysis. All scripts and detailed documentations can be found in [ANALYSIS](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis)*

</br></br>

# [Getting started](#description)

</br></br>

**Getting started with GenomicSEM**

-   Read the [wiki and
    tutorial](https://github.com/GenomicSEM/GenomicSEM/wiki) for the
    R-package `GenomicSEM`

**Install R packages**

-   See
    [here](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis)
    for the list packages to be installed

**Download the required software tools, including**

-   [PLINK](https://www.cog-genomics.org/plink/)
-   [DEPICT](https://data.broadinstitute.org/mpg/depict/)
-   [PASCAL](https://www2.unil.ch/cbg/index.php?title=Pascal)

**Download summary statistic files**

-   See
    [here](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis)
    for the list of summary statistic files included in this analysis

</br></br></br>

# [Model without SNP effects](#model)

</br>

## Estimate genetic correlations

</br></br>

![](results/figures/CorrGWA.svg) </br>

To prepare the data and run all analytical steps involved in deriving
the genetic correlations, follow the script published
[here](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis)

</br></br>

## Specify the structural model

![](results/figures/strucModel.png)

</br></br></br></br>

# [Run the multivariate genome-wide association study](#modelSNP)

![](results/figures/ManHplot_commonLiability.jpeg)

</br></br>

# [Functional follow-up](#functional)

# [LD score regression](#ldsc)
