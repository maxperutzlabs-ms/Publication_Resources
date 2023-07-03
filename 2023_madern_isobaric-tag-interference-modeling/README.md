# Publication: "A causal model of ion interference enables assessment and correction of ratio compression in multiplex proteomics"

This README gives an overview and description of all R code files used in the study. The individual files are located in the various subfolders of this repository, which represent independent experiments and corresponding analyses. Additionally, the table **Metadata_table.xlsx** gives an overview of the Thermo raw files and associated database search result files used in each experiment. Both are available for download on PRIDE (identifier PXD040449).

<br>
<br>

## Code Overview

### Interference Modeling of Yeast-Human Dataset (Main Experiment)
1)  **IM_v1.2__MS2_YeastHuman.Rmd**: Code to perform interference modeling and interference-correction for PSM level data of yeast-human mixture samples measured via MS2-based quantification.
2)  **IM_v1.2__MS2_YeastHuman_PPF-basedCorrection.Rmd**: Same as 1) but 1-PPF was used instead of EIL in the interference-correction algorithm.
3)  **functions_IM.R**: Contains larger functions called in 1) and 2).
4)  **Evaluation_InterferenceCorrectedProteinIntensities.Rmd**: Code for evaluating interference-correction results obtained from 1) and 2) via density plots, ROC-curves and more.
5)  **functions_ROC&Volcano.R**: Contains larger functions called in 4).
6)  **Labeling_Efficiency_Calculation_From_MIX.Rmd** Script to estimate channel-wise TMT labeling efficiencies from mixture of TMT labeling reactions (i.e. before quenching) on the basis of partially labeled peptide intensities.


### Comparative Measurements 
1)  **Interference_Characterization_in_Comparative_Measurements.Rmd**: Code to analyze variations in sample- and measurement-specific parameters in their effect on ion interference and more.
2)  **parameters_sampleComplexity.txt**: Parameters used in 1). Specific to comparing variations in sample complexity.
3)  **parameters_isolationWindowWidth.txt**: Parameters used in 1). Specific to comparing variations in isolation window width.
4)  **parameters_quantificationStrategy.txt**: Parameters used in 1). Specific to comparing variations in quantification strategy (i.e. MS2, FAIMS-MS2, MS3).
5)  **parameters_injectionTime.txt**: Parameters used in 1). Specific to comparing variations in maximum injection time.
6)  **parameters_injectionAmount.txt**: Parameters used in 1). Specific to comparing variations in injection amount.
7)  **parameters_gradientLengths.txt**: Parameters used in 1). Specific to comparing variations in gradient length.


### Interference Modeling of TKO9 Dataset
1)  **IM_v1.2__TKO9.Rmd**: Code to perform interference modeling etc. for PSM level data of TKO9 standard measured via MS2-based quantification.
2)  **functions_IM.R**: Contains larger functions called in 1).


### Interference Modeling of TKO11 Dataset
1)  **IM_v1.2__TKO11.Rmd**: Code to perform interference modeling etc. for PSM level data of TKO9 standard measured via MS2-based quantification.
2)  **functions_IM.R**: Contains larger functions called in 1).


### Interference Modeling of Yeast-Human Dataset (FAIMS-MS2 quantified)
1)  **IM_v1.2__FAIMS_MS2_YeastHuman.Rmd**: Code to perform interference modeling etc. for PSM level data of yeast-human mixture samples measured via FAIMS MS2-based quantification.
2)  **functions_IM.R**: Contains larger functions called in 1).


### Comparison of Isolation Window Purity Metrics 
1)  **IM_Fragpipe.Rmd**: Code to perform interference modeling etc. for PSM level data of yeast-human mixture samples which was database-searched via FrapPipe.
2)  **IM_PD_MSAmanda.Rmd**: Code to perform interference modeling etc. for PSM level data of yeast-human mixture samples which was database-searched via Proteome Discoverer software MSAmanda.
3)  **functions_IM.R**: Contains larger functions called in 1) and 2).


### Site-To-Protein Normalization 
1)  **IM_v1.2__AcetylData_YeastHuman.Rmd**: Code to perform interference modeling for PSM level data of acetyl (K) peptide-enriched yeast-human mixture samples measured via MS2-based quantification.
2)  **IM_v1.2__PhosphoData_YeastHuman.Rmd**: Code to perform interference modeling for PSM level data of phospho (STY) peptide-enriched yeast-human mixture samples measured via MS2-based quantification.
3)  **IM_v1.2__MS2_YeastHuman_Pool2Pool5.Rmd**: Code to perform interference modeling for PSM level data of unmodified peptides (i.e. "proteome") in yeast-human mixture samples measured via MS2-based quantification.
4)  **functions_IM.R**: Contains larger functions called in 1), 2) and 3).
5)  **Process_proteinGroups_v0.7__MS3Proteome.Rmd**: Code to perform filtering, isotopic impurity correction, between-sample normalization etc. for protein level data (in MaxQuant's "proteinGroups.txt") measured via MS3-based quantification.
6)  **Aggregate_PSMs_to_Sites_v0.8__AcetylData_YeastHuman.Rmd**: Code to perform aggregation of PSM level output from 1) to acetyl-site features listed in MaxQuant's "Acetyl (K)Sites.txt".
7)  **Aggregate_PSMs_to_Sites_v0.8__PhosphoData_YeastHuman.Rmd** : Code to perform aggregation of PSM level output from 2) to phospho-site features listed in MaxQuant's "Phospho (STY)Sites.txt".
8)  **Aggregate_PSMs_to_Proteins_v0.6__Pool2Pool5.Rmd**: Code to perform aggregation of PSM level output from 3) to MS2-quantified protein features as listed in MaxQuant's "proteinGroups.txt".
9)  **Normalize_AcetylSiteMS2ToProteinMS3_v0.16.Rmd**: Code to perform site-to-protein normalization of MS2-quantified acetyl-site features (output from 6)) to MS3-quantified protein features (output from 5)).
10) **Normalize_PhosphoSiteMS2ToProteinMS3_v0.16.Rmd**: Code to perform site-to-protein normalization of MS2-quantified phospho-site features (output from 7)) to MS3-quantified protein features (output from 5)).
11) **Normalize_AcetylSiteMS2ToProteinMS2_v0.16.Rmd**: Code to perform site-to-protein normalization of MS2-quantified acetyl-site features (output from 6)) to MS2-quantified protein features (output from 8)).
12) **Normalize_PhosphoSiteMS2ToProteinMS2_v0.16.Rmd**: Code to perform site-to-protein normalization of MS2-quantified phospho-site features (output from 7)) to MS2-quantified protein features (output from 8)).
13) **functions_Site_To_Protein.R**: Contains larger functions called in 9), 10), 11) and 12).

<br>
<br>

## Session Info

```
R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] msqrob2_1.2.0               QFeatures_1.4.0             MultiAssayExperiment_1.20.0 DESeq2_1.34.0              
 [5] SummarizedExperiment_1.24.0 MatrixGenerics_1.6.0        matrixStats_0.61.0          GenomicRanges_1.46.1       
 [9] GenomeInfoDb_1.30.0         IRanges_2.28.0              limma_3.50.0                MSnbase_2.20.0             
[13] ProtGenerics_1.26.0         S4Vectors_0.32.3            mzR_2.28.0                  Rcpp_1.0.7                 
[17] Biobase_2.54.0              BiocGenerics_0.40.0         cowplot_1.1.1               fields_13.3                
[21] viridis_0.6.2               viridisLite_0.4.0           spam_2.8-0                  doParallel_1.0.17          
[25] iterators_1.0.13            foreach_1.5.2               rlist_0.4.6.2               gridExtra_2.3              
[29] MASS_7.3-54                 plot3D_1.4                  pracma_2.3.8                forcats_0.5.1              
[33] stringr_1.4.0               dplyr_1.0.7                 purrr_0.3.4                 readr_2.1.1                
[37] tidyr_1.1.4                 tibble_3.1.6                ggplot2_3.3.6               tidyverse_1.3.1 
[41] sva_3.42.0 				 plotly_4.10.0  			 ROCR_1.0-11           		 IHW_1.22.0   
[45] ggVennDiagram_1.2.0
```

<br>
<br>

## Used Libraries and other Resources

- sva: Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Zhang Y, Storey JD, Torres LC (2022). sva: Surrogate Variable Analysis. R package version 3.46.0.

- Msnbase: Gatto, L. & Lilley, K. S. Msnbase-an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics 28, 288–289 (2012).

- fields: Douglas Nychka, Reinhard Furrer, John Paige, S. S. (2021). “fields: Tools for spatial data.”

- limma: Ritchie, M. E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 43, e47 (2015).

- DESeq2: Anders, S. & Huber, W. Differential expression analysis for sequence count data. Genome Biol. 11, R106 (2010).

- msqrob2: Goeminne, L. J. E., Gevaert, K. & Clement, L. Peptide-level robust ridge regression improves estimation, sensitivity, and specificity in data-dependent quantitative label-free shotgun proteomics. Mol. Cell. Proteomics 15, 657–668 (2016).

- plot3D: Soetaert, K. plot3D: Plotting Multi-Dimensional Data.

- ggplot2: Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. (Springer-Verlag New York).

- MaxQuant: Tyanova, S., Temu, T. & Cox, J. The MaxQuant computational platform for mass spectrometry-based shotgun proteomics. Nat. Protoc. 11, 2301–2319 (2016).

- MS Amanda: Dorfer, V. et al. MS Amanda, a universal identification algorithm optimized for high accuracy tandem mass spectra. J. Proteome Res. 13, 3679–3684 (2014).

- FragPipe: Teo, G. C., Polasky, D. A., Yu, F. & Nesvizhskii, A. I. Fast Deisotoping Algorithm and Its Implementation in the MSFragger Search Engine. J. Proteome Res. 20, 498–505 (2021).

- msConvert: Chambers, M. C. et al. A cross-platform toolkit for mass spectrometry and proteomics. Nat. Biotechnol. 30, 918–20 (2012).

- PRIDE: Perez-Riverol, Y. et al. The PRIDE database resources in 2022: A hub for mass spectrometry-based proteomics evidences. Nucleic Acids Res. 50, D543–D552 (2022).




