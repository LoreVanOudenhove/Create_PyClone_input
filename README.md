# Clonality Assessment

The mutational profile of a cancer genome is characterized by a series of alterations such as single nucleotide substitutions (SNVs), small insertions and deletions (indels), genomic rearrangements, and copy-number variants (CNVs). These genomic changes follow a Darwinian evolution process and can be identified by high-throughput sequencing technologies. Nevertheless, the identification of genomic alterations is complicated by the tumour heterogeneity comprising normal cell contamination and cancer cell subpopultations.  

The identification of these subclonal populations, however, requires a detailed and accurate reconstruction of the cancer genome structure. Mainly the VAF (i.e., the fraction of the reads mapping to a region that shares a specific genotype) is of importance here. Nevertheless, it should be noted that VAF is highly dependent on the local copy number state and tumour purity. By correcting for copy numbers, observed VAFs are typically transformed into cancer-cell fractions (CCF, i.e., the fraction of tumour cells carrying the mutation). To finally identify subclonal populations, the distribution of CCFs is evaluated for distinct clusters that represent individual subpopulations.  

These guidelines provides a workflow to infer the prevalence of somatic mutations in heterogeneous cancer samples from paired tumour-normal NGS data (Figure 2). To do so, copy number changes, tumour purity and ploidy is estimated using Sequenza and somatic mutations are clustered using PyClone. 

## Copy Number, Tumor Purity and Ploidy Estimation using Sequenza. 

[Sequenza](https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html) is an R package that enables the efficient estimation of tumour cellularity and ploidy, and generation of copy number, loss-of-heterozygosity, and mutation frequency profiles. More detailed information on how to run Sequenza can be found in the [Sequenza User Guide](https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html).

### Running Sequenza: Preprocessing of Input Files
- Process a FASTA file to produce a GC Wiggle track file
'''
sequenza−utils gc_wiggle −w 50 --fasta hg38.fa -o hg38.gc50Base.wig.gz
'''   


---
## References
Favero, F., Joshi, T., Marquard, A. M., Birkbak, N. J., Krzystanek, M., Li, Q., … Eklund, A. C. (2015). Sequenza: Allele-specific copy number and mutation profiles from tumor sequencing data. Annals of Oncology, 26(1), 64–70. https://doi.org/10.1093/annonc/mdu479
