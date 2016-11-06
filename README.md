## Frequency Spectrum Composite Likelihood method to detect signatures of selection

Detect selective sweeps in genome-wide genotype or sequencing data, with or without ascertainment bias. 

This program implements the method of Nielsen et al. 2005 for detecting recent strong selection from genome-wide genotype data in a population sample, similar to, and modeled after, Melissa Hubisz program SweepFinder. Shared between the programs and methods is a composite-likelihood based approach to detecting selection by a localized perturbation of the site frequency spectrum (SFS) relative to the genome-wide background SFS, with an increase in rare alleles and high frequency derived alleles. This method is explained fully in

Nielsen R, Williamson S, Kim Y, Hubisz MJ, Clark AG, and Bustamante CD (2005). [Genomic scans for selection sweeps using SNP data.](http://genome.cshlp.org/content/15/11/1566.long) (2005) Genome Research 15(11): 1566-75. [Pubmed 16251466](https://www.ncbi.nlm.nih.gov/pubmed/16251466)

### Permutation test for selection

This program differs from the method described there and from SweepFinder in that it performs a permutation test to assess the statistical significance of the composite likelihood ratio (CLR) of selective sweep models compared to the background frequency spectrum. Permutations are performed using contiguous blocks of the genome which preserve characteristics of background linkage-disequilbrium over short ranges, but break up the much longer extents of linkage-disequilbrium that is created by strong selection. This erases the "footprint" of a selective sweep but retains the background frequency distribution correlation characteristics that require a composite likelihood approach. The permutation is performed repeatedly up to the maximum number requested by the user, but for computational efficiency sites assessed by the program are automatically pruned from the permutation scans once at least 20 permutations have produced higher CLRs at the locus than the original unpermuted data. P-values are extrapolated from the number of permutations that produced CLRs higher than the original data vs. the number of permutations assessed at the locus by the negative binomial distribution.

While the Nielsen et al. paper shows that their original method can produce a very similar null distribution to a simulated null distribution over the major portion of the null distribution mass when various demographic scenarios have affected the population (see Figure 6), it is the tail of the distribution that is most important when assessing selection in a genome wide scan, due to multiplicity of tests. While the method as described performs well for p-values in the range 1 to 1e-2, it diverges sharply from the null distribution expectation under various non-selective population genetic forces further out to the tail past 1e-2. The permutation test produces more accurate calculation of p-values in the presence of population size changes, population structure, and recombination rate variation across the genome when no selection is present than the original method implemented in SweepFinder.


### Genome-wide SNP array data and Ascertainment Bias

As described in Nielsen et al., the program implements ascertainment bias correction when the ascertainment bias is a simple condition such as "variant allele must be seen in at least K genomes from a sample of size M". Thus, the program may be used with genome-wide SNP array data, or whole genome sequencing data. Even if ascertainment conditions for building the array are unknown or are more complicated, the simple condition can be created by choosing M samples and requiring at least K minor alleles to be observed to include the site in the input of the program. As long as the ascertainment conditions to design the array were more permissive than the simple condition for choice of K and M, the effective ascertainment conditions will be the simple condition you choose. The program is then executed with the ascertainment sample depth and minimum minor allele requirement options set accordingly.

### Whole genome sequencing data and rare allele calling bias

Note that an ascertainment condition such as above can be applied to a population sample for whole-genome sequencing data to minimize the potential problem that rare alleles, especially singletons, may be over- or under-called. Choose K>=2 and M>=K, apply the condition, and reduce the dataset to only SNPs with minor alleles occuring at least K times in the M samples you chose. Set the ascertainment bias correction options accordingly. This may be important because the selective sweep model expects an increase in rare alleles around selected sites and thus accuracy in calling rare alleles may strongly influence the model.

### How to use the program

You only need the frequency of the derived allele. Missing data or varying sample depths at each site is automatically handled by the program. Sites where which allele is the derived allele can not be determined can be specified as a "folded" site. The program will understand and treat that site as coming from the "folded site-frequency spectrum". In this case provide the *minor* allele count. Genotype data is not supplied to the program, just the derived allele or minor allele count and the number of samples with observed data at that site, along with chromosome and position of the site. 

Input file format, command line options, and examples coming soon.

### LICENSE

FSCL is licensed under the GNU Public License version 3 (GPLv3).

