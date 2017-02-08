#C - Stephanocoenia intersepta
*Contents:* Details the processing and analyses of the *Stephanocoenia intersepta* dataset, from raw sequences to eventual population genomic analyses. *Notebook author:* Pim Bongaerts. 

*[Click here to go back to the overview](https://github.com/pimbongaerts/bermuda-rad/)*

* **[C1 - PyRAD clustering] (#c1---pyrad-clustering)**
	* [C1a - Raw sequence data] (#c1a---raw-sequence-data)
	* [C1b - PyRAD clustering] (#c1b---pyrad-clustering)
	* [C1c - FASTA with reference loci] (#c1c---fasta-with-reference-loci)
* **[C2 - Filtering and QC of loci and SNPs] (#c2---filtering-and-qc-of-loci-and-snps)**
	* [C2a - Remove *Symbiodinium* and other contamination] (#c2a---remove-symbiodinium-and-other-contamination)
	* [C2b - Basic SNP QC and filtering] (#c2b---basic-snp-qc-and-filtering)
	* [C2c - Remove low-performance individuals and check for clones] (#c2c---remove-low-performance-individuals-and-check-for-clones)
	* [C2d - Minimum representation filter] (#c2d---minimum-representation-filter)
	* [C2e - Check for deviations from HWE] (#c2e---check-for-deviations-from-hwe)
	* [C2f - Final datasets] (#c2f---final-datasets)
	* [C2g - Sample performance] (#c2g---sample-performance)
* **[C3 - Outlier analyses] (#c3---outlier-analyses)**
	* [C3a - Lositan] (#c3a---lositan)
	* [C3b - BayeScan] (#c3b---bayescan)
	* [C3c - Summary] (#c3c---summary)
	* [C3d - Non-outlier dataset] (#c3d---non-outlier-dataset)
* **[C4 - Genetic structure] (#c4---genetic-structure)**
	* [C4a - STRUCTURE] (#c4a---structure)
	* [C4b - PCA] (#c4b---pca)
	* [C4c - GD matrix] (#c4c---gd-matrix)
	* [C4d - DAPC] (#c4d---dapc)
* **[C5 - Additional analyses] (#c5---additional-analyses)**
	* [C5a - Differentiation and diversity stats] (#c5a---differentiation-and-diversity-stats)
	* [C5b - Symbiont typing] (#c5c---symbiont-typing)
	* [C5c - Morphometrics] (#c5d---morphometrics)


##C1 - PyRAD clustering
Note: location of raw sequence data (hosted separately on the [NBCI SRA](https://www.ncbi.nlm.nih.gov/bioproject/361144)) is referred to as `storage_server_path/...`.

###C1a - Raw sequence data
Replace barcodes of `fastq` files by actual sample names [using [fastq_barcodes2samplenames.py](https://github.com/pimbongaerts/radseq/blob/master/fastq_barcodes2samplenames.py) script]:
	
	$ fastq_barcodes2samplenames.py storage_server_path/sint storage_server_path/sint_barcodes.txt

Assess number of reads for all 109 samples and calculate mean [using [fastq_seqcount.py](https://github.com/pimbongaerts/radseq/blob/master/fastq_seqcount.py) script]:

	$ fastq_seqcount.py storage_server_path/sint sint_seqcounts.txt
    mean: 1,383,902 reads
    min: 375,822 reads
    max: 3,791,658 reads

Output number of reads for 4 failed samples (excluded from analyses) [using [fastq_seqcount.py](https://github.com/pimbongaerts/radseq/blob/master/fastq_seqcount.py) script]:

	$ fastq_seqcount.py storage_server_path/sint_failed sint_failed_seqcounts.txt
	mean: 33,283 reads
	min: 3,891 reads
	max: 119,110 reads
    
	$ cat sint_failed_seqcounts.txt
	SIMGD6845H	5917 reads
	SIMJS6731H	3891 reads
	SIMWD6946H	119110 reads
	SIMWD6973H	4217 reads
	
	
###C1b - PyRAD clustering
Run PyRAD clustering pipeline (includes QC, clustering and variant calling):

	$ pyrad -p params.txt -s1234567

Params file (only modified parameters shown; MaxSH "disabled" by giving it a high value; parameter 6 can be left empty for nextRAD data):

	$ head -n 15 params.txt
	==** parameter inputs for pyRAD version 3.0.6  **======================== affected step ==
	./                        ## 1. Working directory                                 (all)
	                          ## 2. Loc. of non-demultiplexed files (if not line 16)  (s1)
	                          ## 3. Loc. of barcode file (if not line 16)             (s1)
	vsearch-1.1.3-linux-x86_64 ## 4. command (or path) to call vsearch (or usearch)    (s3,s6)
	muscle                    ## 5. command (or path) to call muscle                  (s3,s7)
	GTGTAGAGG                 ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)
	28                        ## 7. N processors (parallel)                           (all)
	6                         ## 8. Mindepth: min coverage for a cluster              (s4,s5)
	3                         ## 9. NQual: max # sites with qual < 20 (or see line 20)(s2)
	.90                       ## 10. Wclust: clustering threshold as a decimal        (s3,s6)
	rad	                      ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)
	50                        ## 12. MinCov: min samples in a final locus             (s7)
	200                       ## 13. MaxSH: max inds with shared hetero site          (s7)
	sint                      ## 14. Prefix name for final output (no spaces)         (s7)

Output number of RAD loci in `.loci` file (included as compressed file `sint.loci.gz`):

	$ tail -n 1 storage_server_path/sint/outfiles/sint.loci | cut -d '|' -f 2 # number of loci
    7591
	
Output number of SNPs (given that header of PyRAD vcf is 12 lines) in `.vcf` file:
	
	$ wc -l sint.vcf | awk '{print $1-12}' # number of SNPs
	57869
	
###C1c - FASTA with reference loci
Extract a single reference sequence (using first sample) for each RAD locus [using [pyrad2fasta.py](https://github.com/pimbongaerts/radseq/blob/master/pyrad2fasta.py) script]:

	$ pyrad2fasta.py storage_server_path/sint/outfiles/sint.loci > sint_1c.fa

##C2 - Filtering and QC of loci and SNPs
###C2a - Remove *Symbiodinium* contamination
#####Symbiodinium contamination
Identify *Symbiodinium* contamination through a `blastn` comparison of RAD loci against the three references below (used rather than BWA mapping for increased sensitivity).

* *Symbiodinium* subtraction reference (isolated from *Agaricia fragilis*, *Stephanocoenia intersepta* and other species; as detailed in notebook "[A - *Symbiodinium*](https://github.com/pimbongaerts/bermuda-rad/tree/master/A%20-%20Symbiodinium)")
*  *Symbiodinium* C1 draft genome (in preparation by [CX Chan](http://researchers.uq.edu.au/researcher/1283) et al. as part of the [ReFuGe2020](http://refuge2020.com) sequencing initiative)
*  *Symbiodinium minutum* B1 genome (as published in [Shoguchi et al. (2013) *Current Biology*](http://dx.doi.org/10.1016/j.cub.2013.05.062); available through: <http://marinegenomics.oist.jp>)

Create output directory:

	$ mkdir local_blastn_results

Comparison (`blastn`) against the three *Symbiodinium* references: 
	
	$ OUTPUT_FORMAT="7 qseqid sseqid length nident pident evalue bitscore"
	$ blastn -query sint_1c.fa -db SymRAD16 -task blastn -outfmt $OUTPUT_FORMAT -out local_blastn_results/sint_2a_blastn_SymRAD16.txt
	$ blastn -query sint_1c.fa -db SymC -task blastn -outfmt $OUTPUT_FORMAT -out local_blastn_results/sint_2a_blastn_SymC.txt
	$ blastn -query sint_1c.fa -db SymB -task blastn -outfmt $OUTPUT_FORMAT -out local_blastn_results/sint_2a_blastn_SymB.txt

Extract `blastn` matches with an E-value lower than 10<sup>-15</sup> [using [mapping_get_blastn_matches.py](https://github.com/pimbongaerts/radseq/blob/master/mapping_get_blastn_matches.py) script]:

	$ MAX_E_VALUE="0.000000000000001"
	$ mapping_get_blastn_matches.py local_blastn_results/sint_2a_blastn_SymRAD16.txt $MAX_E_VALUE
	Matches: 1086 | Min.length: 47.0 bp | Min. nident: 47.0 bp | Min. pident: 79.12 %
	$ mapping_get_blastn_matches.py local_blastn_results/sint_2a_blastn_SymC.txt $MAX_E_VALUE
	Matches: 35 | Min.length: 62.0 bp | Min. nident: 60.0 bp | Min. pident: 81.82 %
	$ mapping_get_blastn_matches.py local_blastn_results/sint_2a_blastn_SymB.txt $MAX_E_VALUE
	Matches: 11 | Min.length: 59.0 bp | Min. nident: 57.0 bp | Min. pident: 85.11 %

Get list of all loci matching *Symbiodinium* references:

	$ cat local_blastn_results/sint_2a_blastn_SymRAD16_match0.000000000000001.txt \
		  local_blastn_results/sint_2a_blastn_SymC_match0.000000000000001.txt \
		  local_blastn_results/sint_2a_blastn_SymB_match0.000000000000001.txt \
		  | cut -f1 | sort | uniq > symbiodinium_loci_to_remove.txt
	$ wc -l loci_to_remove.txt # number of matching loci
	1088 loci_to_remove.txt

#####Other potential contamination
Identify potential microbial contamination by `blastn` comparison against non-redundant NCBI database.

Create output directory:

	$ mkdir nt_blastn_results

Compare (`blastn`) against the NCBI non-redundant database (this time with `taxids` and `title` included in output format and with a 10<sup>-5</sup> E-value treshold):

	$ OUTPUT_FORMAT="7 qseqid sseqid length nident pident evalue bitscore staxids stitle"
	$ blastn -query sint_1c.fa -remote -db nt -task blastn -evalue 0.0001 -max_target_seqs 10 -outfmt $OUTPUT_FORMAT -out nt_blastn_results/sint_2a_blastn_nt.txt
	
Extract matches and identify the phylum from the `taxids` (using Entrez Direct to connect to NCBI's taxonomic database) [using [mapping_identify_blast_matches.py](https://github.com/pimbongaerts/radseq/blob/master/mapping_identify_blast_matches.py) script]:

	$ mapping_identify_blast_matches.py nt_blastn_results/sint_2a_blastn_nt.txt 0.0001 [email@adress.com]
	$ cut -f 8 nt_blastn_results/sint_2a_blastn_nt_match0.0001.txt | sort | uniq -c > other_loci_hit_summary.txt

Create a list of loci that match phyla outside the Cnidaria:

	$ grep -v "Cnidaria" nt_blastn_results/sint_2a_blastn_nt_match0.0001.txt | cut -f 1 > other_loci_to_remove.txt

#####Remove potential contamination
Remove loci matching *Symbiodinium* and non-Cnidaria from `.vcf` file [using [vcf_remove_chrom.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_remove_chrom.py) script]:

	$ cat symbiodinium_loci_to_remove.txt other_loci_to_remove.txt | sort | uniq > all_loci_to_remove.txt
	$ wc -l all_loci_to_remove.txt
		1223 all_loci_to_remove.txt
	$ vcf_remove_chrom.py sint.vcf all_loci_to_remove.txt > sint_2a.vcf

###C2b - Basic SNP QC and filtering
Remove singletons, doubletons (i.e. SNPs where the minor allele only occurs in a single individual and that individual is homozygotic for that allele) and monomorphic SNPs:
	
	$ vcftools --vcf sint_2a.vcf --singletons
	$ tail -n +2 out.singletons | cut -f 1-2 > eliminate_2b.txt # remove header and cut out first two columns
	$ vcftools --vcf sint_2a.vcf --exclude-positions eliminate_2b.txt --min-alleles 2 --recode --stdout > sint_2b_temp1.vcf
	$ tail -n 4 out.log
	After filtering, kept 109 out of 109 Individuals
	Outputting VCF file...
	After filtering, kept 27149 out of a possible 48373 Sites
	Run Time = 1.00 seconds

Evaluate distribution of SNPs along position in RAD contigs [using [vcf_pos_count.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_pos_count.py) script]:

	$ vcf_pos_count.py sint_2b_temp1.vcf
	1-5		1284	*****************************
	6-10	1306	*****************************
	11-15	1340	******************************
	16-20	1275	****************************
	21-25	1343	******************************
	26-30	1291	*****************************
	31-35	1331	******************************
	36-40	1336	******************************
	41-45	1311	*****************************
	46-50	1329	******************************
	51-55	1316	*****************************
	56-60	1311	*****************************
	61-65	1373	*******************************
	66-70	1286	*****************************
	71-75	1340	******************************
	76-80	1326	******************************
	81-85	1547	***********************************
	86-90	2324	****************************************************
	91-95	1995	*********************************************
	96-100	172	***
	101-105	13

Given overrepresentation of SNPs towards end of contigs, eliminate all SNPs at POS > 80 bp [using [vcf_read_trim.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_read_trim.py) script]:

	$ vcf_read_trim.py sint_2b_temp1.vcf 80 > sint_2b.vcf
	
Create a popfile from the current `.vcf` file [using [popfile_from_vcf.py](https://github.com/pimbongaerts/radseq/blob/master/popfile_from_vcf.py) script]:

	$ popfile_from_vcf.py sint_2b.vcf 4 5 > popfile_2b.txt

###C2c - Remove low-performance individuals and check for clones
Assessed missing data for individuals [using [vcf_missing_data.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_missing_data.py) script]:

	$ vcf_missing_data.py sint_2b.vcf | tail -n +1 | sort -n --key=5 | head -n 10
	INDIVIDUAL	MISS	GENO	TOTAL	% GENOTYPED
	SIMJS6720H	20993	105	21098	0.5
	SIMJD6893H	20138	960	21098	4.55
	SIMGD6828H	17507	3591	21098	17.02
	SICB155507H	17206	3892	21098	18.45
	SIMJD6905H	16138	4960	21098	23.51
	SIMWD6962H	15237	5861	21098	27.78
	SIMPD7018H	14222	6876	21098	32.59
	SIMPD7036H	13988	7110	21098	33.7
	SIMGD6848H	13784	7314	21098	34.67

Removed two individuals with lots of missing data (`SIMJS6720H` and `SIMJD6893H`):

	$ vcftools --vcf sint_2b.vcf --remove-indv SIMJS6720H --remove-indv SIMJD6893H --recode --stdout > sint_2c_temp1.vcf

Calculate allelic similarity between all pairs of individuals from unfiltered `.vcf` file to assess potential clones [using [vcf_clone_detect.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_clone_detect.py) script]:
	
	$ vcf_clone_detect.py -v sint_2c_temp1.vcf -p popfile_2b.txt
	###1 - Pairwise comparisons of all individuals
	5671 comparisons completed
	
	###2 - Histogram (of pairwise genetic similarities)
	 75       2 **
	 76      30 ******************************
	 77      50 **************************************************
	 78      20 ********************
	 79      68 ********************************************************************
	 80     141 *********************************************************************#
	 81     721 *********************************************************************#
	 82    1849 *********************************************************************#
	 83    2059 *********************************************************************#
	 84     663 *********************************************************************#
	 85      67 *******************************************************************
	 86       0
	 87       0
	 88       1 *
	 89       0
	
	###3 - List of highest matches
	88.22	0	[B1]	SICB155507H vs SICB155378H	2183.5/2475	11239	12334
	88.0	-------------------- Potential threshold --------------------
	85.73	2.49	[GD-JD]	SIMGD6828H vs SIMJD6884H	2208.5/2576	12607	11067
	85.66	0.07	[JD-WD]	SIMWD6962H vs SIMJD6905H	2069.5/2416	18554	4960
	85.65	0.01	[PD-WD]	SIMWD6962H vs SIMPD7015H	3302.5/3856	14224	10730
	85.64	0.01	[JD-PS]	SIMJD6905H vs SIMPS6654H	3089.0/3607	12279	12426
	85.64	0.0	[GD-JD]	SIMGD6841H vs SIMJD6905H	2958.0/3454	19592	4960
	85.62	0.02	[PD-PS]	SIMPD7015H vs SIMPS6658H	6058.5/7076	15576	12598
	85.54	0.08	[GD-PD]	SIMGD6828H vs SIMPD7025H	2160.0/2525	13787	9836
	85.53	0.01	[JD-WD]	SIMWD6955H vs SIMJD6905H	2861.0/3345	19483	4960
	85.52	0.01	[GD-WD]	SIMWD6962H vs SIMGD6833H	3236.0/3784	14752	10130
	85.52	0.0	[PD-PS]	SIMPD7025H vs SIMPS6654H	5642.5/6598	15270	12426
	85.52	0.0	[GD-WD]	SIMGD6841H vs SIMWD6962H	3463.5/4050	19287	5861
	85.51	0.01	[PS]	SIMPS6652H vs SIMPS6654H	4869.0/5694	14366	12426
	85.47	0.04	[GS-JD]	SIMGS6771H vs SIMJD6905H	3076.0/3599	19737	4960
	85.46	0.01	[JD-PD]	SIMPD7036H vs SIMJD6898H	4863.5/5691	11291	15498
	85.41	0.05	[PD]	SIMPD7018H vs SIMPD7015H	3802.5/4452	14820	10730
	85.39	0.02	[GD-PD]	SIMGD6828H vs SIMPD7015H	2124.5/2488	12856	10730
	85.38	0.01	[GD-PD]	SIMPD7015H vs SIMGD6848H	3996.5/4681	18465	7314
	85.37	0.01	[PD]	SIMPD7015H vs SIMPD7022H	6822.5/7992	14534	14556
	85.34	0.03	[GD-GS]	SIMGS6768H vs SIMGD6848H	5176.0/6065	19849	7314
	85.32	0.02	[GD-PS]	SIMPS6642H vs SIMGD6848H	4564.0/5349	19133	7314
	85.32	0.0	[GS]	SIMGS6771H vs SIMGS6788H	6080.5/7127	16872	11353
	[...truncated...]
	
	###4 - Clonal groups (threshold: 88.0)
	B1: SICB155507H, SICB155378H (88.22 %)

Only the Curacao samples are pulled out as being a lot more similar to each other (likely a filtering artefact to some extent). Given that these are not included in population-based analyses, no further removal of samples:

	$ mv sint_2c_temp1.vcf sint_2c.vcf

###C2d - Minimum representation filter
Eliminate SNPs in `.vcf` that are genotyped for less than 50% of individuals in each population (rather than using an overall missing-data filter) [using [vcf_minrep_filter.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_minrep_filter.py) script]:

	$ vcf_minrep_filter.py sint_2c.vcf popfile_withCU_2d.txt 0.5 sint_2d.vcf
	13288 out of 21098 SNPs failing 0.5 threshold
	Lowest overall proportion genotyped of remaining SNPs: 0.5607476635514018
	
Before continuing with further QC (where Curacao individuals are not considered) - assess whether any alleles are fixed/private to Bermuda (note: regions were inserted into output below in square brackets) [using [vcf_contrast_samples.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_contrast_samples.py) script]:

	$ vcf_contrast_samples.py sint_2d.vcf SICB155378H SICB155507H
	7810 SNPs found in VCF file.
	52 SNPs were fixed in the non-reference [Bermuda] samples,
	of which 27 SNPs were alternatively fixed in the reference samples [Curacao]
	and 25 SNPs were variable in the reference samples [Curacao].
	
###C2e - Check for deviations from HWE
Removed Curacao individuals from dataset (and manually removed them from popfile):

	$ vcftools --vcf sint_2d.vcf --remove-indv SICB155378H --remove-indv SICB155507H --recode --stdout > sint_2e_temp1.vcf

Convert `.vcf` file to Arlequin format [using [vcf_spider.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_spider.py) script]:
	
	$ vcf_spider.py sint_2e_temp1.vcf popfile_2e.txt sint_2e.arp

Run HWE test for all 7 populations using `arlecore` (exact test using a Markov chain with chain length = 1 000 000 and dememorization steps = 100 000):

	$ arlecore sint_2e.arp arl_run.ars

Scan output for 'suspicious' SNPs, that meet one of these two criteria :

* Significantly deviating from HWE (excess or deficit) in min. 5 out of 7 populations. No Bonferroni/FDR applied as using multiple-population criteria (observed in min. 5 pops), so currently no SNPs expected to be identified just by chance (7,810 SNPs * 0.05^5). This still allows SNPs to remain in the dataset that deviate from HWE due to being under selection (associated with a particular population or habitat; never > 4 pops)
* Exhibiting a significant excess of heterozygotes in at least 1 population (to filter out potential paralogous loci). NB: Although this may be be due to overdominant selection or outbreeding.

Scrape XML output for SNPs meeting above criteria (using P-value cut-off of 0.05; `.vcf` file included to obtain original CHROM/POS identifiers) [using [E - Scripts/bermuda_arlequin_hwe_summary.py](../E\ -\ Scripts/bermuda_arlequin_hwe_summary.py) script]:
	
	$ python3 bermuda_arlequin_hwe_summary.py -i sint_2e.xml -p 0.05 -d 5 -v sint_2e_temp1.vcf
	Total of 7810 SNPs evaluated
	52 SNPs monomorphic across all populations
	2868 SNPs with a significant deficit in heterozygosity in at least one population
	
	Suspicious SNPs/loci (67):
	26 SNPs with a significant excess in heterozygosity in at least one population
	41 addditional SNPs significantly deviating from HWE in 5 or more populations

Remove loci/CHROMs (in `sint_2e_elimchrom_uniq.txt`)  and individual SNPs (in `sint_2e_elimpos.txt `) from `.vcf` file [using [vcf_remove_chrom.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_remove_chrom.py) script]:

	$ cat sint_2e_elimchrom.txt | sort | uniq > sint_2e_elimchrom_uniq.txt # reduce to only uniq loci
	$ wc -l sint_2e_elimchrom_uniq.txt sint_2e_elimpos.txt
      16 sint_2e_elimchrom_uniq.txt
      41 sint_2e_elimpos.txt
	$ vcf_remove_chrom.py sint_2e_temp1.vcf sint_2e_elimchrom_uniq.txt > sint_2e_temp2.vcf
	$ vcftools --vcf sint_2e_temp2.vcf --exclude-positions sint_2e_elimpos.txt --recode --stdout > sint_2e.vcf

###C2f - Final datasets
Dataset (no clones) multi-allelic, without Curacao samples `sint_2f.vcf` (50% max. missing data per pop; monomorphic sites [related to exclusion of Curacao] removed):
	
	$ vcftools --vcf sint_2e.vcf --mac 1  --recode --stdout > sint_2f.vcf
	$ tail -n 4 out.log
	After filtering, kept 105 out of 105 Individuals
	Outputting VCF file...
	After filtering, kept 7655 out of a possible 7707 Sites
	Run Time = 0.00 seconds

Dataset (no clones) bi-allelic, without Curacao samples `sint_bi_2f.vcf` (50% max. missing data per pop; removing 111 multi-allelic sites):
	
	$ vcftools --vcf sint_2e.vcf --min-alleles 2 --max-alleles 2 --mac 1 --recode --stdout > sint_bi_2f.vcf
	$ tail -n 4 out.log
	After filtering, kept 105 out of 105 Individuals
	Outputting VCF file...
	After filtering, kept 7547 out of a possible 7707 Sites
	Run Time = 0.00 seconds

Dataset (no clones) bi-allelic with Curacao samples `sint_withCU_2f.vcf`:

	$ tail -n +13 sint_2e.vcf | cut -f 1,2 > filtered_snps.txt
	$ vcftools --vcf sint_2d.vcf --positions filtered_snps.txt --min-alleles 2 --max-alleles 2 --keep popfile_withCU_2d.txt --recode --stdout > sint_withCU_2f.vcf
	$ tail -n 4 out.log
	After filtering, kept 107 out of 107 Individuals
	Outputting VCF file...
	After filtering, kept 7599 out of a possible 7810 Sites
	Run Time = 0.00 seconds

###C2g - Sample performance
Output number of genotyped and missing SNPs for each sample [using [vcf_missing_data.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_missing_data.py) script]:

	$ vcf_missing_data.py sint_2f.vcf > sint_2f_stats.txt

##C3 - Outlier analyses
###C3a - Lositan
Convert no-clones dataset `sint_bi_2f.vcf` to Lositan/Genepop format [using [vcf_spider.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_spider.py) script]:

	$ vcf_spider.py sint_bi_2f.vcf popfile_2f.txt sint_3a.lositan

Run the `Lositan` GUI (on Linux machine due to unresolved OSX bug) to assess outliers for the following population datasets:

|filename|description|populations|
|---|---|---|
|`fdist_shallow.txt`|within shallow between locations|PS, JS, GS|
|`fdist_deep.txt`|within deep between locations|PD, JD, GD|
|`fdist_PS-PD.txt`|within Princess between depths|PS, PD|
|`fdist_JS-JD.txt`|within JohnSmith between depths|JS, JD|
|`fdist_GS-GD.txt`|within Gurnet between depths|GS, GD|
|`fdist_overall.txt`|between all populations|PS, PD, JS, JD, GS, GD, WD|

###C3b - BayeScan
Converted no-clones dataset `sint_bi_2f.vcf` to BayeScan format [using [vcf_spider.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_spider.py) script]:

	$ vcf_spider.py sint_bi_2f.vcf popfile_2f.txt sint_3b.bayescan

Ran `bayescan` on the overall dataset:

	$ bayescan sint_3b.bayescan -threads 28
	$ cp sint_3b.baye_fst.txt bayescan_overall.txt

###C3c - Summary
Summarise outlier analyses in a single table with boolean values (`TRUE` = significant outlier) [using [E - Scripts/bermuda_outliers_summary.py](../E\ -\ Scripts/bermuda_outliers_summary.py) script]:

	$ python3 bermuda_outliers_summary.py sint_bi_2f.vcf outlier_files/fdist_overall.txt outlier_files 0.01 0.99 > outlier_summary.txt

Classify outliers types [using [E - Scripts/bermuda_outliers_classify.py](../E\ -\ Scripts/bermuda_outliers_classify.py) script]:
	
	$ python3 bermuda_outliers_classify.py outlier_summary.txt > outlier_summary_types.txt
	$ tail -n +2 outlier_summary_types.txt | cut -f 13 | sort | uniq -c | sort -n
	   4 FDIST_AND_BAYESCAN_OUTLIER
	 134 FDIST_OUTLIER
	7409 NON_OUTLIER
	
No depth-associated outliers were identified (matching all 3 criteria below):

* identified as overall outlier by both `fdist` AND `bayescan`
* identified as outliers in at least 2 individual between-depth comparisons
* never identified in any of the between-location (within-depth) comparisons

Extract one representative "overall" outlier SNP per CHROM from summary file (after assessing that SNPs on same CHROM show same pattern):

	$ grep "FDIST_AND_BAYESCAN_OUTLIER" outlier_summary_types.txt | cut -f 2,3 | sort -u -k1,1 > overall_outliers.txt
	$ wc -l overall_outliers.txt
	      4 overall_outliers.txt

Extract an additional 8 FDist outliers with the highest Fsts

	$ grep "FDIST_OUTLIER" outlier_summary_types.txt | sort -rn --key=5 | head -n 8 | cut -f 2,3 > additional_examples.txt
	$ cat overall_outliers.txt additional_examples.txt > outlier_snps_to_plot.txt

Calculate genotype frequencies of selected outliers [using [vcf_genotype_freqs.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_genotype_freqs.py) script]:
      
	$ vcf_genotype_freqs.py sint_bi_2f.vcf factors_3c.txt outlier_snps_to_plot.txt > sint_3c_freqs.txt

###C3d - Non-outlier dataset
Summarise outlier analyses but with less stringent threshold, into a single table with boolean values (`TRUE` = significant outlier) [using [E - Scripts/bermuda_outliers_summary.py](../E\ -\ Scripts/bermuda_outliers_summary.py) script and [E - Scripts/bermuda_outliers_classify.py](../E\ -\ Scripts/bermuda_outliers_classify.py) script]:

	$ python3 bermuda_outliers_summary.py sint_bi_2f.vcf outlier_files/fdist_overall.txt outlier_files 0.05 0.95 > outlier_summary.txt
	$ python3 bermuda_outliers_classify.py outlier_summary.txt > outlier_summary_types.txt
	$ tail -n +2 outlier_summary_types.txt | cut -f 13 | sort | uniq -c | sort -n
	   9 FDIST_AND_BAYESCAN_OUTLIER
	 432 FDIST_OUTLIER
	7106 NON_OUTLIER

Extract non-outlier SNPs (by eliminating all SNPs that were an overall outlier for both `fdist` and `bayescan`) and created bi-allelic dataset without-clones and outlier SNPs: `sint_neut_3d.vcf`:

	$ grep "NON_OUTLIER" outlier_summary_types.txt | cut -f 2,3 > neutral_snps.txt
	$ vcftools --vcf sint_bi_2f.vcf --positions neutral_snps.txt --recode --stdout > sint_neut_3d.vcf
	$ tail -n 4 out.log
	After filtering, kept 105 out of 105 Individuals
	Outputting VCF file...
	After filtering, kept 7106 out of a possible 7547 Sites
	Run Time = 0.00 seconds

##C4 - Genetic structure
###C4a - STRUCTURE
#####Overall dataset
Modify STRUCTURE's `mainparams` file (manually) and enter the number of unique CHROMs and total number of individuals:

	$ tail -n +13 sint_2f.vcf | cut -f 1 | sort | uniq | wc -l
	    2187 # number of unique CHROMs

Run overall dataset (`sint_2f.vcf `) through STRUCTURE (for K = 2..7) [using the multi-threading wrapper [structure_mp](https://github.com/pimbongaerts/radseq/blob/master/structure_mp.py)]:

	$ structure_mp.py sint_2f.vcf popfile_2f.txt 7 10 10
	Initialise individuals and populations...
	Subsample SNPs (one random SNP per locus)... [2187 SNPs/loci]
	Outputting 10 STRUCTURE files...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 2 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 3 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 4 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 5 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 6 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 7 ...10 reps DONE
	Running CLUMPP on replicates for K = 2 ...
	Running CLUMPP on replicates for K = 3 ...
	Running CLUMPP on replicates for K = 4 ...
	Running CLUMPP on replicates for K = 5 ...
	Running CLUMPP on replicates for K = 6 ...
	Running CLUMPP on replicates for K = 7 ...
	K = 2: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1
	K = 3: MedMeaK 0.5 MaxMeaK 1 MedMedK 0.5 MaxMedK 1
	K = 4: MedMeaK 1.0 MaxMeaK 2 MedMedK 1.0 MaxMedK 2
	K = 5: MedMeaK 1.0 MaxMeaK 2 MedMedK 1.0 MaxMedK 2
	K = 6: MedMeaK 1.0 MaxMeaK 2 MedMedK 1.0 MaxMedK 2
	K = 7: MedMeaK 0.0 MaxMeaK 2 MedMedK 0.0 MaxMedK 2

#####"Neutral" dataset
Modify STRUCTURE's `mainparams` file (manually) and enter the number of unique CHROMs and total number of individuals:

	$ tail -n +13 sint_neut_3d.vcf | cut -f 1 | sort | uniq | wc -l
	    2153	# number of unique CHROMs

Run the non-outlier ("neutral") dataset (`sint_neut_3d.vcf `) through STRUCTURE [using the multi-threading wrapper [structure_mp](https://github.com/pimbongaerts/radseq/blob/master/structure_mp.py)]:

	$ structure_mp.py sint_neut_3d.vcf popfile_2f.txt 4 10 10
	Initialise individuals and populations...
	Subsample SNPs (one random SNP per locus)... [2153 SNPs/loci]
	Outputting 10 STRUCTURE files...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 2 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 3 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 4 ...10 reps DONE
	Running CLUMPP on replicates for K = 2 ...
	Running CLUMPP on replicates for K = 3 ...
	Running CLUMPP on replicates for K = 4 ...
	K = 2: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1
	K = 3: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1
	K = 4: MedMeaK 1.0 MaxMeaK 2 MedMedK 1.0 MaxMedK 2

#####90th percentile
Extract 90th percentile using overall F<sub>ST</sub> values as calculated by Lositan (using bi-allelic datafile `sint_bi_2f.vcf`) [using [vcf_splitfst.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_splitfst.py) script]:

	$ vcf_splitfst.py sint_bi_2f.vcf fdist_overall.txt 90 100 
	Fst minimum: -0.048957
	Fst maximum: 0.281867
	Fst 90th percentile: 0.05144740000000002
	Fst 100th percentile: 0.281867
	754 SNPs between 90 and 100 percentile

Modify STRUCTURE's `mainparams` file (manually) and enter the number of unique CHROMs and total number of individuals:

	$ tail -n +13 sint_bi_2f_fst90_100.vcf | cut -f 1 | sort | uniq | wc -l
	504 	# number of unique CHROMs

Run the 90th percentile dataset (`sint_bi_2f_fst90_100.vcf`) through STRUCTURE [using the multi-threading wrapper [structure_mp](https://github.com/pimbongaerts/radseq/blob/master/structure_mp.py)]:

	$ structure_mp.py sint_bi_2f_fst90_100.vcf popfile_2f.txt 7 10 10
	Initialise individuals and populations...
	Subsample SNPs (one random SNP per locus)... [504 SNPs/loci]
	Outputting 10 STRUCTURE files...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 2 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 3 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 4 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 5 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 6 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 7 ...10 reps DONE
	Running CLUMPP on replicates for K = 2 ...
	Running CLUMPP on replicates for K = 3 ...
	Running CLUMPP on replicates for K = 4 ...
	Running CLUMPP on replicates for K = 5 ...
	Running CLUMPP on replicates for K = 6 ...
	Running CLUMPP on replicates for K = 7 ...
	K = 2: MedMeaK 2.0 MaxMeaK 2 MedMedK 2.0 MaxMedK 2
	K = 3: MedMeaK 3.0 MaxMeaK 3 MedMedK 3.0 MaxMedK 3
	K = 4: MedMeaK 4.0 MaxMeaK 4 MedMedK 4.0 MaxMedK 4
	K = 5: MedMeaK 5.0 MaxMeaK 5 MedMedK 5.0 MaxMedK 5
	K = 6: MedMeaK 6.0 MaxMeaK 6 MedMedK 6.0 MaxMedK 6
	K = 7: MedMeaK 7.0 MaxMeaK 7 MedMedK 7.0 MaxMedK 7

#####90th percentile randomized
Shuffle around population assignment of individuals [using [popfile_toggleassign.py](https://github.com/pimbongaerts/radseq/blob/master/popfile_toggleassign.py) script]:

	$ popfile_toggleassign.py popfile_2f.txt > toggle_popfile_2f.txt

Converted dataset with shuffled population file `toggle_popfile_2f.txt` to Lositan/Genepop format (using wrapper for `PGD spider`) - and calculated F<sub>ST</sub> values with Lositan GUI (output file: `shuffle_overall.txt`) [using [vcf_spider.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_spider.py) script]:

	$ vcf_spider.py sint_bi_2f.vcf toggle_popfile_2f.txt sint_rand_4a.lositan

Extract 90th percentile using overall F<sub>ST</sub> values as calculated by Lositan [using [vcf_splitfst.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_splitfst.py) script]:

	$ vcf_splitfst.py sint_bi_2f.vcf shuffle_overall.txt 90 100
	Fst minimum: -0.055631
	Fst maximum: 0.313581
	Fst 90th percentile: 0.051805
	Fst 100th percentile: 0.313581
	755 SNPs between 90 and 100 percentile
	$ mv sint_bi_2f_fst90_100.vcf sint_bi_rand_2f_fst90_100.vcf

Modify STRUCTURE's `mainparams` file (manually) and enter the number of unique CHROMs and total number of individuals:

	$ tail -n +13 sint_bi_rand_2f_fst90_100.vcf | cut -f 1 | sort | uniq | wc -l
	512	 # number of unique CHROMs

Run the shuffled 90th percentile dataset (`sint_bi_rand_2f_fst90_100.vcf`) through STRUCTURE [using the multi-threading wrapper [structure_mp](https://github.com/pimbongaerts/radseq/blob/master/structure_mp.py)]:

	$ structure_mp.py sint_bi_rand_2f_fst90_100.vcf toggle_popfile_2f.txt 7 10 10
	Initialise individuals and populations...
	Subsample SNPs (one random SNP per locus)... [512 SNPs/loci]
	Outputting 10 STRUCTURE files...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 2 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 3 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 4 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 5 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 6 ...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 7 ...10 reps DONE
	Running CLUMPP on replicates for K = 2 ...
	Running CLUMPP on replicates for K = 3 ...
	Running CLUMPP on replicates for K = 4 ...
	Running CLUMPP on replicates for K = 5 ...
	Running CLUMPP on replicates for K = 6 ...
	Running CLUMPP on replicates for K = 7 ...
	K = 2: MedMeaK 2.0 MaxMeaK 2 MedMedK 2.0 MaxMedK 2
	K = 3: MedMeaK 3.0 MaxMeaK 3 MedMedK 3.0 MaxMedK 3
	K = 4: MedMeaK 4.0 MaxMeaK 4 MedMedK 4.0 MaxMedK 4
	K = 5: MedMeaK 5.0 MaxMeaK 5 MedMedK 5.0 MaxMedK 5
	K = 6: MedMeaK 6.0 MaxMeaK 6 MedMedK 6.0 MaxMedK 6
	K = 7: MedMeaK 7.0 MaxMeaK 7 MedMedK 7.0 MaxMedK 7

###C4b - PCA
Run a PCA on the overall dataset (`sint_bi_2f.vcf`) using the `adegenet` package in R [detailed in [E - Scripts/bermuda_vcf2pca.R](../E\ -\ Scripts/bermuda_vcf2pca.R) script]:
	
	$ Rscript bermuda_vcf2pca.R sint_bi_2f.vcf popfile_pca_4b.csv

###C4c- GD matrix
Calculate GD matrix (Hamming-distance) between individual samples (of dataset including Curacao) [using [vcf_gdmatrix.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_gdmatrix.py) script]:

	$ vcf_gdmatrix.py sint_withCU_2f.vcf popfile_depth_4c.txt > sint_depth_withCU_gd_4c.txt

###C4d - DAPC
Run a DAPC on the overall dataset (`sint_bi_2f.vcf`) using the `adegenet` package in R [detailed in [E - Scripts/bermuda_vcf2dapc.R](../E\ -\ Scripts/bermuda_vcf2dapc.R) script]:
	
	$ Rscript DAPC_4d.R sint_bi_2f.vcf popfile_pca_4b.csv

##C5 - Additional analyses
###C5a - Differentiation and diversity stats
Reduced missing data by retaining only loci that are genotyped for at least 80% of individuals in each population [using [vcf_minrep_filter.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_minrep_filter.py) script]:

	$ vcf_minrep_filter.py sint_bi_2f.vcf popfile_2f.txt 0.8 sint_bi_5a.vcf
		5994 out of 7547 SNPs failing 0.8 threshold
		Lowest overall proportion genotyped of remaining SNPs: 0.8285714285714286

Calculate several statistics regarding differentiation in R using `hierfstat` and `adegenet` packages (as well as `vcfR` for importing `.vcf` file) [detailed in [E - Scripts/bermuda_diff_stats.R](../E\ -\ Scripts/bermuda_diff_stats.R) script]:

	* General F-statistics and measures of heterozygosity
	* G-statistic test for overall population structure
	* Non-parametric Kruskal-Wallis tests to assess for an overall difference between groups of populations
	* Non-parametric Wilcoxon signed rank tests to assess for differences between groups of populations


	$ Rscript bermuda_diff_stats.R sint_bi_5a.vcf popfile_pca_4b.csv "#eff3ff" "#c6dbef" "#c6dbef"
	Processed variant: 7547

	Weir and Cockerham Fst and Fis:
	$FST
	[1] 0.0008113305	
	$FIS
	[1] 0.2806109
	
	Expected heterozygosity (Hs) for each pop:
	       GD        GS        JD        JS        PD        PS        WD
	0.1618783 0.1635142 0.1622547 0.1660883 0.1622498 0.1627970 0.1620870
	
	G-statistic test for overall structure:
	Monte-Carlo test
	Observation: 107590.7	
	Based on 99 replicates
	Simulated p-value: 0.3
	Alternative hypothesis: greater
	     Std.Obs  Expectation     Variance
	5.463167e-01 1.058580e+05 1.005967e+07
	
	G-statistic test for overall structure - eastern pops only:
	Monte-Carlo test
	Observation: 91105.17
	Based on 99 replicates
	Simulated p-value: 0.17
	Alternative hypothesis: greater
	     Std.Obs  Expectation     Variance
	1.124875e+00 8.804102e+04 7.420139e+06
	
	Pairwise Fst when averaged over all loci (without NA):
	       Habt.G        Habt.J        Habt.P       Shal.GJ       Shal.JP
	-0.0010787434 -0.0014363090 -0.0009118337 -0.0002122199 -0.0011505544
	      Shal.GP       Deep.GJ       Deep.JP       Deep.GP       West.GW
	-0.0013791576 -0.0020147638  0.0003923556 -0.0006621781 -0.0025852681
	      West.JW       West.PW
	-0.0012456374 -0.0007682172
	
	Pairwise Fst:
	           1          2          3          4          5          6
	2 0.02849031
	3 0.03097026 0.02924617
	4 0.03033748 0.02866072 0.03033715
	5 0.03020273 0.02845650 0.03150790 0.02851057
	6 0.02791249 0.02463953 0.02750751 0.02655480 0.02663476
	7 0.02945874 0.02848730 0.03148605 0.03001721 0.02977629 0.02716061

###C5b - Symbiont typing
Representative samples of *Stephanocoenia intersepta* (n = 33) from the 7 populations were sequenced for the *Symbiodinium COX1* region. Two haplotypes were recovered:

	>Stephanocoenia_COX1_HAPB
	GCAACAATCTTTTCTCCATTCAAGATTATTCCAGAGAAGATAGCAATTATAGCTCCTAAAGAAAGAACAAAATGAA
	AATGTGCTACAACATAATATGTATCATGTAATCCTAGATCCACTGCACCATTTCCAAGAATTATTCCTGTTGACCC
	CCCTATCGTAAACATTAATAAAAAGAGATGTGAGAGGAAGACAGAAGTAATTCTAAGGTGTAATAATGGTGGATTT
	GAGAGATATGTAAAAAGCCAATTAAAGATTTTTGTACCAGTTGGTAAGGATATTAAGATTGTAACTCCTGTAAAAT
	AAGCTCTTGTATCACTTTCTAAACCTACAGTATACATATGATGTCCCCAAACAAGACCTCCAAGAAGAGAAATAGA
	TGACATGGCAAAGATCATTGATTGGTTAGCAAAGATTATTAACTGTAAAATACCAGAAATTATTATGGAAATGATC
	CCAAATGCAGGAATTATTAATATGTAAACTTCTGGATGTCCAAAAAACCAAAATAAATGTTGATAGAATATAGGAT
	CTCCTCCAAATATTGGATCAAAGAAAAGTGTATTAGAATGAAGATCACCCAATATTAAAAGAAGTGTACCAGATAA
	GATTGGTAATGTTAATAAAAGCATGAAAGCTGTAATCAAGAAAGCCCAAGGAAATAATGGGATAGAAGATAATATC
	AGATAATAAGATCTCAGAAAATGAATTGTTGTCCAAAAGTTAAGAGATGTAAGACATGAAGATATACCAGAGATTA
	TTAATCCAAATATAAGATTTCCTGTACTTGAAGGTGATAAAGTCA
	
	>Stephanocoenia_COX1_HAPC
	CCAACAATCTTTTCTCCATTAAAAATTATTCCAGAGAAGATAGAAATTATAGCTCCTAAAGAAAGAACAAAATGAAAA
	TGTGCTACAACATAATATGTATCATGTAATCCTAGATCTACTGCAGCATTTCCAAGAATTATTCCTGTAGAACCACCT
	ATAGTAAACATTAATAAAAAGAGATGTGAGAAGAAGACAGAACTAATTCTAAGGTGTAATAATGGTGGATTTCCAAGA
	TATGTAAAAAGCCAATTAAAGATTTTTGTACCAGTTGGTAAGGATATTAAGATTGTAACTCCTGTAAAATAAGCTCTT
	GTATCACTTTCTAAACCTACAGTATACATATGATGCCCCCAAACAAGACTTCCAAGAAGAGAAATAGATGACATGGCA
	AAGATCATTGATTGGTTACCAAAGATTATTAATTGTAAAATCCCAGAAATTATTATGGAAATGATCCCAAATGCAGGA
	ATTATTAAGATGTAAACTTCTGGATGTCCAAAAAACCAAAATAAATGTTGATAGAATATAGGATCTCCTCCAAATATA
	GGATCAAAGAAAAGTGTATTAGAATGAAGATCACCCAATATTAAAAGAAGTGTACCAGATAAGATTGGTAATGTTAAT
	AAAAGCATGAAAGCTGTAATCAAGAAAGCCCAAGGAAATAATGGGATAGTCTTTAATGTCAGATAATAAGATCTCAGA
	TTTAGGATTGTTGTCCAAAAGTTAAGAGATGTAAGACATGAAGATATACCAGAGATTAATAATCCAAATATAAGATTT
	CCTGTACTTGAAGGTGATAAAGTCA

In addition, 3 RAD loci (that were also recoved in *Agaricia fragilis*) matched *Symbiodinium* chloroplast minicircle sequences in the `blastn` search conducted in step **C2a - Remove Symbiodinium and other contamination**. 

Extract those loci that are in common from the PyRAD `.loci` files and trim to 75bp [using [pyrad_filter.py](https://github.com/pimbongaerts/radseq/blob/master/pyrad_filter.py) script and [pyrad_trim.py](https://github.com/pimbongaerts/radseq/blob/master/pyrad_trim.py) script]:

	$ echo "3967\n4580\n4729" > sint_sym_chloroplast_loci.txt
	$ pyrad_filter.py storage_server_path/sint/outfiles/sint.loci sint_sym_chloroplast_loci.txt 0 0 > sint_sym_chloroplast.loci
	$ pyrad_trim.py sint_sym_chloroplast.loci 75 > sint_sym_chloroplast_trim.loci

Visual assessment of `sint_sym_chloroplast_trim.loci` alignment identifies 3 different haplotypes (one representing an intermediate haplotype between the other two):

	>Stephanocoenia_CHL_HAPA
	TGATATATCAGCTTGGGATTCCTTTTATCTAGCAACATTCTGGATGCTTAATAGCAACACATGGATAAGCTTCTACTTCCACTACAAGCACC
	>Stephanocoenia_CHL_HAPB
	TGATATATCAGCTTGGGATTCCTTTTATCTAGCAACGTTCTGGATGCTTAATACCAACACATGGATAAGCTTCTATTTCCACTACAAGTACC
	>Stephanocoenia_CHL_HAPAB
	TGATATATCAGCTTGGGATTCCTTTTATCTAGCAACRTTCTGGATGCTTAATASCAACACATGGATAAGCTTCTAYTTCCACTACAAGYACC
	
Mitochondrial and chloroplast haplotypes for each sample are listed in: `sint_symbionts.txt`.

###B5c - Morphometrics
Corallite diameter and density were measured for representative samples of *Stephanocoenia intersepta* (n = 34) with 5 pseudoreplicate measurements per sample. Averaged values (across the 5 measurements) for each specimen are listed in `sint_morpho_data.csv`.