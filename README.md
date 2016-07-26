# Project Epsilon
Biostatistics research project in Dr.Mar's lab at the Albert Einstein College of Medicine.
Benjamin Church and Henry Williams

##Methods:

####RNA Expression Datasets:
Datasets were combined into two analysis groups, one group collected via microarray chips and the other via RNA sequencing (RNA-Seq).

#####Microarray Datasets:
Six microarray datasets were collected, five cancer groups and one control group. An ovarian serious cystadenocarcinoma (OV) group (n = 568), a glioblastoma multiforme (GBM) group (n = 548), and a Luminal A (LumA) breast cancer subgroup (n = 284) were taken from The Cancer Genome Atlas (TCGA)[1]. Samples were collected from primary solid tumor and recurrent solid tumor and are log2 transformed and normalized via RPKM. These data can be accessed from the TCGA level 3 database. All three datasets were collected at University of North Carolina at Chapel Hill (UNC) using an Agilent 244K Custom Gene Expression Microarray (G4502A-07-3).

Three microarray datasets were taken from the NCBI Gene Expression Omnibus (GEO). We used two groups with acute myeloid leukemia (AML), one of exclusively individuals over the age of 60 (GSE6891) [2] (n = 461) with samples collected from both blood and bone marrow and the second (GSE15434) of exclusively normal karyotype (NK) AML (n = 251) [3, 4] with samples collected from mononuclear cells. These two datasets were collected using an Affymetrix Human Genome U133 Plus 2.0 Array. Our control for the microarray group is a HapMap expression profiling (GSE6536) [5] collected via a Sentrix Human-6 Expression BeadChip. The HapMap expression data were log2 transformed and normalized using quantile normalization [5].

####RNA-Seq Datasets:
We analyzed six RNA-Seq datasets, five cancer groups and one control group. Five cancer groups were taken from TCGA level 3 data. These were: skin cutaneous melanoma (SKCM) (n = 470), head and neck squamous cell carcinoma (HNSC) (n = 519), lower grade glioma (LGG) (n = 514), 
lung squamous cell carcinoma (LUSC) (n = 495), and kidney renal clear cell carcinoma (KI C) (n = 531). All five datasets contain samples from primary solid tumor and recurrent solid tumor collected at UNC using Illumina HiSeq 2000 RNA Sequencing Version 2 Analysis and log2 transformed. [6]

Our RNA-Seq control is RNA-sequencing of 465 lymphoblastoid cell lines from the 1000 Genomes project created by the Geuvadis consortium. [7] These data were collected using an Illumina HiSeq 2000, processed with GEM mapper 1.349, and log2 transformed.
   
 
####Batch Effect Correction:	
Microarray TCGA datasets were corrected for batch effects using median/standard deviation correction [8]. However, accurate batch correction could not be applied to the GEO datasets due to unmarked batches. Than being said, 76 extra samples added to GSE6891 at a later date (total number of samples: 537) were removed from the analysis after showing significant batch effects. The RNA-Seq TCGA datasets were batch corrected using a linear model from the R/Bioconductor package limma (version 3.24.15 run on R version 3.1.2).  

####Skewness Statistic:

Our chosen metric of skewness is calculated by dividing the cube root of the third moment of a distribution by its standard deviation. Explicitly, the relative skewness of a gene's transcription expression (g) over a population X, is calculated as:

This statistic was chosen to differentiate between wide slightly asymmetric distributions and narrow highly asymmetric distributions. Alone, the third moment cannot distinguish between these two qualitatively different distributions (Fig. 1). For each dataset, the distribution of this statistic applied over that sample population is calculated for the available genome. The genes are divided into two groups, the genes with positive skew and those with negative skew. The percentage of genes in each of these groups is referred to as the gene splitting of that particular dataset. The empirical results show that very few genes have zero skew in any dataset i.e. few distributions are symmetric.

####Hypothesis Tests:
A Shapiro-Wilk test was performed on the gene splitting data to assess the normality of the variation between datasets obtained using the same amplification method. A Welch two sample t-test was used to test the hypothesis that gene splitting is equally distributed in microarray and RNA-Seq datasets.

####Sample Size Analysis:
We analyzed the effect of sample size on the gene splitting as a surrogate for the sensitivity of skewness results to sample size. The gene splitting statistic was computed on random subsamples of the various datasets. The gene splitting is organized into a sequence indexed by increasing sample size. We used two metrics to assess the convergence of each sequence, the rate of convergence and the error to the limit. For an infinite sequence sn with limit L, the rate of convergence (C) is a number between zero and one defined as:	
And the nth error (E) as:    

However, only a finite section of the sampling sequence is known so limits cannot be directly computed. We assume that the convergence can be approximated by an exponential decrease of the form: with constants a, b, and L. The constants were determined using least-squares regression. Once these constants are determined, the limit of the sequence is approximately L and the rate of convergence is: .
