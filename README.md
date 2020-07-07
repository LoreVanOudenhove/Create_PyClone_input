# Clonality Assessment

The mutational profile of a cancer genome is characterized by a series of alterations such as single nucleotide substitutions (SNVs), small insertions and deletions (indels), genomic rearrangements, and copy-number variants (CNVs). These genomic changes follow a Darwinian evolution process and can be identified by high-throughput sequencing technologies. Nevertheless, the identification of genomic alterations is complicated by the tumour heterogeneity comprising normal cell contamination and cancer cell subpopultations.  

The identification of these subclonal populations, however, requires a detailed and accurate reconstruction of the cancer genome structure. Mainly the VAF (i.e., the fraction of the reads mapping to a region that shares a specific genotype) is of importance here. Nevertheless, it should be noted that VAF is highly dependent on the local copy number state and tumour purity. By correcting for copy numbers, observed VAFs are typically transformed into cancer-cell fractions (CCF, i.e., the fraction of tumour cells carrying the mutation). To finally identify subclonal populations, the distribution of CCFs is evaluated for distinct clusters that represent individual subpopulations.  

These guidelines provides a workflow to infer the prevalence of somatic mutations in heterogeneous cancer samples from paired tumour-normal NGS data. To do so, copy number changes and tumour purity are estimated using **Sequenza** and somatic mutations are clustered using **PyClone**. To transform the results generetad by Sequenza into a valid input file for PyClone, the *Sequenza_to_PyClone.py* function can be used.

## Copy Number and Tumor Purity Estimation using Sequenza.

[Sequenza](https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html) is an R package that enables the efficient estimation of tumour cellularity and ploidy, and generation of copy number, loss-of-heterozygosity, and mutation frequency profiles. More detailed information on how to run Sequenza can be found in the [Sequenza User Guide](https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html).

### Preprocessing of Input Files
* Process a FASTA file to produce a GC Wiggle track file:
```
sequenza−utils gc_wiggle −w 50 --fasta hg38.fa -o hg38.gc50Base.wig.gz
```
* Process BAM and Wiggle files to produce a seqz file:
```
sequenza−utils bam2seqz -n normal.bam -t tumor.bam --fasta hg38.fa \
    -gc hg38.gc50Base.wig.gz -o out.seqz.gz
```
* Post-process by binning the original seqz file:
```
sequenza−utils seqz_binning --seqz out.seqz.gz -w 50 -o out small.seqz.gz
```

### Sequenza Analysis (in R)

```
library(sequenza)

data.file <- "small.seqz.gz"
analysis <- sequenza.extract(data.file, verbose = FALSE)
CP <- sequenza.fit(analysis)

# Write files and plots
sequenza.results(sequenza.extract = analysis, cp.table = CP, sample.id = "run1")

# Get confidence intervals and extract purity and ploidy
confint <- get.ci(CP)
purity <- confint$max.cellularity
```
This Sequenza analysis will generate a number of files and plots, among which the Run1_segments.txt file that contains detected segments, with estimated copy number state at each segment. This Run1_segments.tsv file, along with the tumour purity esitmate, will be used an input for PyClone.

# Create input for PyClone

(PyClone)[https://bitbucket.org/aroth85/pyclone/wiki/Home] is a statistical model for inference of clonal population structures in cancers. Its algorithm relies on a Bayesian clustering method for grouping sets of deeply sequenced somatic mutations into putative clonal clusters while estimating their cellular prevalences and accounting for allelic imbalances introduced by segmental copy-number changes and normal-cell contamination.

To run the PyClone analysis we can use the ```run_analysis_pipeline``` command. This command requires a tab delimited input file with the subsequent fields: 

* mutation_id - A unique ID to identify the mutation. 
* ref_counts - The number of reads covering the mutation which contain the reference allele. 
* var_counts - The number of reads covering the mutation which contain the variant allele. 
* normal_cn - The copy number at the locus in normal cells. For autosomal chromosomes this will be 2 and for sex chromosomes it could be either 1 or 2. 
* minor_cn - The minor copy number at the locus in cancer cells. 
* major_cn - The major copy number at the locus in cancer cells.

To obtain this tab delimited file, the Sequenza_to_PyClone.py function, provided in this repository can be used.
```
Sequenza_to_PyClone.py -i Run1_segments.txt -v variants.vcf -o PyClone_input.tsv
```
# Clonal Reconstruction using PyClone

Finally, the PyClone_input.tsv file and the estimated tumour purity can be used as an input for PyClone.

```
PyClone run_analysis_pipeline --in_files PyClone_input.tsv --tumour_contents $purity --prior major_copy_number
```

PyClone generates a number of tables and plots, summarising the identified clusters and the associated cellular prevalences. 
![PyClone output](https://github.com/LoreVanOudenhove/Create_PyClone_input/blob/master/images/density.png)

---
## References
Favero, F., Joshi, T., Marquard, A. M., Birkbak, N. J., Krzystanek, M., Li, Q., … Eklund, A. C. (2015). Sequenza: Allele-specific copy number and mutation profiles from tumor sequencing data. Annals of Oncology, 26(1), 64–70. https://doi.org/10.1093/annonc/mdu479

Roth, A., Khattra, J., Yap, D., Wan, A., Laks, E., Biele, J., Ha, G., Aparicio, S., Bouchard-Côté, A., & Shah, S. P. (2014). PyClone: statistical inference of clonal population structure in cancer. Nature methods, 11(4), 396–398. https://doi.org/10.1038/nmeth.2883
