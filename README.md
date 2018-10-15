# Admixture Estimator

The code presented here calculates the admixture estimates given in:

Cox MP, TM Karafet, JS Lansing, H Sudoyo and MF Hammer. 2010. [https://doi.org/10.1098/rspb.2009.2041](Autosomal and X-linked single nucleotide polymorphisms reveal a steep Asian-Melanesian ancestry cline in eastern Indonesia and a sex bias in admixture rates). *Proceedings of the Royal Society B* 277:1589-1596.

The original admixture estimator was developed by:

Chakraborty R, MI Kamboh, M Nwankwo and RE Ferrell. 1992. [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1682537/](Caucasian genes in American Blacks: New data). *American Journal of Human Genetics* 50:145–155.

Specifically equation 6 (see Chakraborty *et al* (1992) for details):

![Chakraborty Equation 6](Chakraborty_Equation6.jpg  | width=100)

The [example code](admixture_estimator.R), written in base R, shows an analysis of Asian-Papuan admixture in the Rindi population of Sumba, eastern Indonesia.

While any data can be used with this code, it is intended to be used with Ancestry Informative Markers (AIMs) that carry substantial information about genetic ancestry.  The basic data structures required by the code are vectors of markers, the number of genotypes screened for each marker and the marker allele frequency relative to the reference population (here, the Asian variant):


|       | Marker 1 | Marker 2 | Marker 3 | Marker 4 | ... | Marker *n* |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Marker | SNP1_chr2_3996214 | SNP2_chr7_152691894 | SNP3_chr11_20238944 | SNP4_chr13_64798694 | ... | SNP39_chr22_1984584 |
| Number of Genotypes | 20 | 20 | 20 | 20 | ... | 20 |
| Reference Allele Frequency | 0.85 | 0.85 | 0.65 | 0.60 | ... | 0.60 |

The same information is needed for proxy parental populations.  In the example code, these are Han Chinese (Asian) and Papua New Guinea Highlanders (Papuan).

The code calculates the admixture estimator of Chakraborty *et al* (1992). However, exact admixture estimates can be extremely dependent on the dataset, especially for small sample sizes. To place statistical bounds on this uncertainty, the code generates random allele frequency probability densities by pulling samples of observed size *n* from the observed frequency profile. These simulated densities are generated for the admixture population, as well as both proxy parental populations.

The output are mean and median estimates of the admixture proportion, as well as assymetric 95% confidence intervals.

> Note: This code was designed for small numbers of ancestry informative markers and will not scale well. For very large numbers of markers (e.g., from SNP genotyping arrays), alternative approaches might be preferable.
