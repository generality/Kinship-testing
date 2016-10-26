Kinship-testing
===
PERL scripts for kinship testing by SNPs and/or STRs based on likelihood ratio (LR) model.
---

#LR calculation

Based on genotypes of autosomal SNPs and/or STRs from real cases or simulated pedigrees, LR for kinship testing was calculated by Elston-Stewart algorithm [1]. For a person within a given pedigree, the model compared likelihood values (L) based on genotypes of autosomal markers (G) of two alternative hypothesis: H0, the test person was the specific member in the relationship pedigree, and H1, the test person was unrelated. Namely, LR = L(G|H0)/L(G|H1), and log10(LR) at each loci can be added for unlinked markers. For STRs and SNPs, two distinct scripts were used.

##For STRs
Allele frequencies were collected from Bingbing Xie's data for Chinese population [2]. The stepwise STR mutation model of Brenner’s method was used [3], and locus-dependent one-step paternal mutation rates of STRs except D19S433 were obtained from the STRBase [4] and the American Association of Blood Banks [5]. For D19S433, a median rate of 0.0017 per trio was applied.

##For SNPs
Allele frequencies were collected from variant data of the East Asia (EAS) continental population in the 1000 Genome Project [6], and a mutation rate of 1.2×10-8 per bp per generation was used [7].

##Datasets
Datasets of SNP and STR genotypes from real cases were also included.

#Pedigree simulation

To simulate the SNP and STR genotypes of individuals in a given pedigree scenario, first the putative pedigree was generated from ancestors , whose genotypes were randomly generated according on population allele frequencies and mutation rates used above. Then, following the rules of genetic inheritance and considering a mutation rate, genotypes of pedigree members were derived, and pairs in relationships of interest were selected for kinship testing. For unrelated controls, the genotypes were directly generated based on allele frequencies.

#Fixation index

The fixation index, F_RT, measures the variance of allele frequencies between continental populations. For each SNP, F_RT was calculated based on genotypes of the 853 unrelated individuals from the 1000 Genome Project [6]. SNPs with a <= 0.03 F_RT are retained.

#Reference
[1] Elston R C, Stewart J F. A General Model for the Genetic Analysis of Pedigree Data[J]. Human Heredity, 1971, 21(6): 523-542.

[2] Xie B, Chen L, Yang Y, et al. Genetic distribution of 39 STR loci in 1027 unrelated Han individuals from Northern China[J]. Forensic Science International-genetics, 2015: 205-206.

[3] http://dna-view.com/mufeatur.htm

[4] http://www.cstl.nist.gov/strbase/

[5] http://www.aabb.org

[6] The 1000 Genomes Project Consortium, A global reference for human genetic variation, Nature. 526 (2015) 68–74.

[7] Scally A, Durbin R. Revising the human mutation rate: implications for understanding human evolution[J]. Nature Reviews Genetics, 2012, 13(10): 745-753.

