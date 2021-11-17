Genome-wide association study on the common liability to addiction
================

</br></br></br></br></br>

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

</br>

**Getting started with GenomicSEM**

-   Read the [wiki and
    tutorial](https://github.com/GenomicSEM/GenomicSEM/wiki) for the
    R-package `GenomicSEM`

**Install R packages**

-   See
    [here](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis#install-r-packages)
    for the list packages to be installed

**Download the required software tools, including**

-   [PLINK](https://www.cog-genomics.org/plink/)
-   [DEPICT](https://data.broadinstitute.org/mpg/depict/)
-   [PASCAL](https://www2.unil.ch/cbg/index.php?title=Pascal)

**Download summary statistic files**

-   See
    [here](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis#download-summary-statistic-files)
    for the list of summary statistic files included in this analysis

</br></br></br>

# [Pre-processing of the summary statistic files](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis#pre-processing-of-the-summary-statistic-files)

-   includes formatting of the summary statistic files
-   running the `munge()` function to prepare the data for LD score
    regression analysis

</br></br></br>

# [Estimate genetic correlations](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis#estimate-genetic-correlations)

![](results/figures/CorrGWA.svg) </br></br>

# [Specify the structural model](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis#specify-the-structural-model)

![](results/figures/strucModel.png)

</br></br>

## [Run the multivariate genome-wide association study](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis#run-the-multivariate-genome-wide-association-study)

![](results/figures/ManHplot_commonLiability.jpeg)

</br></br>

## [Pasthway analysis in PASCAL and DEPICT](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis)

![](results/figures/pascalPlot_comb.svg)

</br></br>

## [LD score regression analysis including other traits](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis#ld-score-regression-analysis-including-other-traits)

![](results/figures/PlotLDScore.svg)

</br></br>

## [Mendelian Randomization analysis](https://github.com/TabeaSchoeler/TS2021_CommonLiabAddiction/tree/master/analysis#mendelian-randomization-analysis)

</br></br>
