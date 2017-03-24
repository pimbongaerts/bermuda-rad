# B - Agaricia fragilis
*Contents:* Details the processing and analyses of the *Agaricia fragilis* dataset, from raw sequences to eventual population genomic analyses. *Notebook author:* Pim Bongaerts. 

*[Click here to go back to the overview](https://github.com/pimbongaerts/bermuda-rad/)*

* **[B1 - PyRAD clustering](#b1---pyrad-clustering)**
	* [B1a - Raw sequence data](#b1a---raw-sequence-data)
	* [B1b - PyRAD clustering](#b1b---pyrad-clustering)
	* [B1c - FASTA with reference loci](#b1c---fasta-with-reference-loci)
* **[B2 - Filtering and QC of loci and SNPs](#b2---filtering-and-qc-of-loci-and-snps)**
	* [B2a - Remove *Symbiodinium* and other contamination](#b2a---remove-symbiodinium-and-other-contamination)
	* [B2b - Basic SNP QC and filtering](#b2b---basic-snp-qc-and-filtering)
	* [B2c - Check for clones](#b2c---check-for-clones)
	* [B2d - Minimum representation filter](#b2d---minimum-representation-filter)
	* [B2e - Check for deviations from HWE](#b2e---check-for-deviations-from-hwe)
	* [B2f - Final datasets](#b2f---final-datasets)
	* [B2g - Sample performance](#b2g---sample-performance)
* **[B3 - Outlier analyses](#b3---outlier-analyses)**
	* [B3a - Lositan](#b3a---lositan)
	* [B3b - BayeScan](#b3b---bayescan)
	* [B3c - Summary](#b3c---summary)
	* [B3d - Non-outlier dataset](#b3d---non-outlier-dataset)
	* [B3e - Outlier annotation](#b3e---outlier-annotation)
* **[B4 - Genetic structure](#b4---genetic-structure)**
	* [B4a - STRUCTURE](#b4a---structure)
	* [B4b - PCA](#b4b---pca)
	* [B4c - GD matrix](#b4c---gd-matrix)
	* [B4d - DAPC](#b4d---dapc)
* **[B5 - Additional analyses](#b5---additional-analyses)**
	* [B5a - Differentiation and diversity stats](#b5a---differentiation-and-diversity-stats)
	* [B5b - Admixture](#b5b---admixture)
	* [B5c - Symbiont typing](#b5c---symbiont-typing)
	* [B5d - Morphometrics](#b5d---morphometrics)

## B1 - PyRAD clustering
Note: location of raw sequence data (hosted separately on the [NBCI SRA](https://www.ncbi.nlm.nih.gov/bioproject/361144)) is referred to as `storage_server_path/...`.

### B1a - Raw sequence data
Replace barcodes of `fastq` files by actual sample names [using [fastq_barcodes2samplenames.py](https://github.com/pimbongaerts/radseq/blob/master/fastq_barcodes2samplenames.py) script]:
	
	$ fastq_barcodes2samplenames.py storage_server_path/afra [storage_server_path/afra_barcodes.txt]

Assess number of reads for all 104 samples and calculate mean [using [fastq_seqcount.py](https://github.com/pimbongaerts/radseq/blob/master/fastq_seqcount.py) script]:

	$ fastq_seqcount.py storage_server_path/afra afra_seqcounts.txt
	mean: 1,347,958 reads
	min: 508,890 reads
	max: 4,627,981 reads

Output number of reads for 6 failed samples (excluded from analyses) [using [fastq_seqcount.py](https://github.com/pimbongaerts/radseq/blob/master/fastq_seqcount.py) script]:

	$ fastq_seqcount.py storage_server_path/afra_failed afra_failed_seqcounts.txt
	mean: 56,587 reads
	min: 725 reads
	max: 282,873 reads
	
	$ cat afra_failed_seqcounts.txt
	AFMGS6754H	282873 reads
	AFMJD6872H	2645 reads
	AFMPS6606H	2018 reads
	AFMPS6607H	947 reads
	AFMPS6609H	725 reads
	AFMWD6929H	50316 reads
	
	
### B1b - PyRAD clustering
Run PyRAD clustering pipeline (includes QC, clustering and variant calling):

	$ pyrad -p params.txt -s1234567

Params file (only modified parameters shown; MaxSH "disabled" by giving it a high value; parameter 6 can be left empty for nextRAD data):

	$ head -n 15 params.txt
	==** parameter inputs for pyRAD version 3.0.6  **======================== affected step ==
	./                        ## 1. Working directory                                 (all)
	                          ## 2. Loc. of non-demultiplexed files (if not line 16)  (s1)
	                          ## 3. Loc. of barcode file (if not line 16)             (s1)
	vsearch-1.1.3-linux-x86_64  ## 4. command (or path) to call vsearch (or usearch)    (s3,s6)
	muscle                    ## 5. command (or path) to call muscle                  (s3,s7)
	GTGTAGAGG                 ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)
	28                        ## 7. N processors (parallel)                           (all)
	6                         ## 8. Mindepth: min coverage for a cluster              (s4,s5)
	3                         ## 9. NQual: max # sites with qual < 20 (or see line 20)(s2)
	.90                       ## 10. Wclust: clustering threshold as a decimal        (s3,s6)
	rad	                      ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)
	50                        ## 12. MinCov: min samples in a final locus             (s7)
	200                       ## 13. MaxSH: max inds with shared hetero site          (s7)
	afra                      ## 14. Prefix name for final output (no spaces)         (s7)

Output number of RAD loci in `.loci` file (included as compressed file `afra.loci.gz`):

	$ tail -n 1 storage_server_path/afra/outfiles/afra.loci | cut -d '|' -f 2 # number of loci
	12145
	
Output number of SNPs (given that header of PyRAD vcf is 12 lines) in `.vcf` file (included in repository):
	
	$ wc -l afra.vcf | awk '{print $1-12}' # number of SNPs
	56244
	
### B1c - FASTA with reference loci
Extract a single reference sequence (using first sample) for each RAD locus [using [pyrad2fasta.py](https://github.com/pimbongaerts/radseq/blob/master/pyrad2fasta.py) script]:

	$ pyrad2fasta.py storage_server_path/afra/outfiles/afra.loci > afra_1c.fa

## B2 - Filtering and QC of loci and SNPs
### B2a - Remove *Symbiodinium* and other contamination
##### Symbiodinium contamination
Identify *Symbiodinium* contamination through a `blastn` comparison of RAD loci against the three references below (used rather than BWA mapping for increased sensitivity).

* *Symbiodinium* subtraction reference (isolated from *Agaricia fragilis*, *Stephanocoenia intersepta* and other species; as detailed in notebook "[A - *Symbiodinium*](https://github.com/pimbongaerts/bermuda-rad/tree/master/A%20-%20Symbiodinium)")
*  *Symbiodinium* C1 draft genome (in preparation by [CX Chan](http://researchers.uq.edu.au/researcher/1283) et al. as part of the [ReFuGe2020](http://refuge2020.com) sequencing initiative)
*  *Symbiodinium minutum* B1 genome (as published in [Shoguchi et al. (2013) *Current Biology*](http://dx.doi.org/10.1016/j.cub.2013.05.062); available through: <http://marinegenomics.oist.jp>)

Create output directory:

	$ mkdir local_blastn_results

Compare (`blastn`) against the three *Symbiodinium* references: 
	
	$ OUTPUT_FORMAT="7 qseqid sseqid length nident pident evalue bitscore"
	$ blastn -query afra_1c.fa -db SymRAD16 -task blastn -outfmt $OUTPUT_FORMAT -out local_blastn_results/afra_2a_blastn_SymRAD16.txt
	$ blastn -query afra_1c.fa -db SymC -task blastn -outfmt $OUTPUT_FORMAT -out local_blastn_results/afra_2a_blastn_SymC.txt
	$ blastn -query afra_1c.fa -db SymB -task blastn -outfmt $OUTPUT_FORMAT -out local_blastn_results/afra_2a_blastn_SymB.txt

Extract `blastn` matches with an E-value lower than 10<sup>-15</sup> [using [mapping_get_blastn_matches.py](https://github.com/pimbongaerts/radseq/blob/master/mapping_get_blastn_matches.py) script]:

	$ MAX_E_VALUE="0.000000000000001"
	$ mapping_get_blastn_matches.py local_blastn_results/afra_2a_blastn_SymRAD16.txt $MAX_E_VALUE
	Matches: 1183 | Min.length: 45.0 bp | Min. nident: 45.0 bp | Min. pident: 76.09 %
	$ mapping_get_blastn_matches.py local_blastn_results/afra_2a_blastn_SymC.txt $MAX_E_VALUE
	Matches: 105 | Min.length: 55.0 bp | Min. nident: 55.0 bp | Min. pident: 81.82 %
	$ mapping_get_blastn_matches.py local_blastn_results/afra_2a_blastn_SymB.txt $MAX_E_VALUE
	Matches: 11 | Min.length: 56.0 bp | Min. nident: 54.0 bp | Min. pident: 82.42 %

Create list of all loci matching *Symbiodinium* references:

	$ cat local_blastn_results/afra_2a_blastn_SymRAD16_match0.000000000000001.txt \
		  local_blastn_results/afra_2a_blastn_SymC_match0.000000000000001.txt \
		  local_blastn_results/afra_2a_blastn_SymB_match0.000000000000001.txt \
		  | cut -f1 | sort | uniq > symbiodinium_loci_to_remove.txt
	$ wc -l symbiodinium_loci_to_remove.txt # number of matching loci
	1188 symbiodinium_loci_to_remove.txt

##### Other potential contamination
Identify potential microbial contamination by `blastn` comparison against non-redundant NCBI database.

Create output directory:

	$ mkdir nt_blastn_results

Compare (`blastn`) against the NCBI non-redundant database (this time with `taxids` and `title` included in output format and with a 10<sup>-5</sup> E-value treshold):

	$ OUTPUT_FORMAT="7 qseqid sseqid length nident pident evalue bitscore staxids stitle"
	$ blastn -query afra_1c.fa -remote -db nt -task blastn -evalue 0.0001 -max_target_seqs 10 -outfmt "7 qseqid sseqid length nident pident evalue bitscore staxids stitle" -out nt_blastn_results/afra_2a_blastn_nt.txt
	
Extract matches and identify the phylum from the `taxids` (using Entrez Direct to connect to NCBI's taxonomic database) [using [mapping_identify_blast_matches.py](https://github.com/pimbongaerts/radseq/blob/master/mapping_identify_blast_matches.py) script]:

	$ mapping_identify_blast_matches.py nt_blastn_results/afra_2a_blastn_nt.txt 0.0001 [email@adress.com]
	$ cut -f 8 nt_blastn_results/afra_2a_blastn_nt_match0.0001.txt | sort | uniq -c > other_loci_hit_summary.txt

Create a list of loci that match phyla outside the Cnidaria:

	$ grep -v "Cnidaria" nt_blastn_results/afra_2a_blastn_nt_match0.0001.txt | cut -f 1 > other_loci_to_remove.txt

##### Remove potential contamination
Remove loci matching *Symbiodinium* and non-Cnidaria from `.vcf` file [using [vcf_remove_chrom.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_remove_chrom.py) script]:

	$ cat symbiodinium_loci_to_remove.txt other_loci_to_remove.txt | sort | uniq > all_loci_to_remove.txt
	$ wc -l all_loci_to_remove.txt
	    1485 all_loci_to_remove.txt
	$ vcf_remove_chrom.py afra.vcf all_loci_to_remove.txt > afra_2a.vcf

### B2b - Basic SNP QC and filtering
Remove singletons, doubletons (i.e. SNPs where the minor allele only occurs in a single individual and that individual is homozygotic for that allele) and monomorphic SNPs:
	
	$ vcftools --vcf afra_2a.vcf --singletons
	$ tail -n +2 out.singletons | cut -f 1-2 > eliminate_2b.txt	# remove header and cut out first two columns
	$ vcftools --vcf afra_2a.vcf --exclude-positions eliminate_2b.txt --min-alleles 2 --recode --stdout > afra_2b_temp1.vcf
	$ tail -n 4 out.log
	After filtering, kept 104 out of 104 Individuals
	Outputting VCF file...
	After filtering, kept 19899 out of a possible 49753 Sites
	Run Time = 1.00 seconds

Evaluate distribution of SNPs along position in RAD contigs [using [vcf_pos_count.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_pos_count.py) script]:

	$ vcf_pos_count.py afra_2b_temp1.vcf
	1-5		320	**
	6-10	348	**
	11-15	333	**
	16-20	322	**
	21-25	325	**
	26-30	345	**
	31-35	358	**
	36-40	354	**
	41-45	358	**
	46-50	313	**
	51-55	331	**
	56-60	370	***
	61-65	366	**
	66-70	345	**
	71-75	335	**
	76-80	541	****
	81-85	1548	************
	86-90	5414	********************************************
	91-95	6422	****************************************************
	96-100	821	******
	101-105	27
	106-107	3

Given overrepresentation of SNPs towards end of contigs, eliminate all SNPs at POS > 75 bp [using [vcf_read_trim.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_read_trim.py) script]:

	$ vcf_read_trim.py afra_2b_temp1.vcf 75 > afra_2b.vcf

Create a popfile from the current `.vcf` file [using [popfile_from_vcf.py](https://github.com/pimbongaerts/radseq/blob/master/popfile_from_vcf.py) script]:

	$ popfile_from_vcf.py afra_2b.vcf 4 5 > popfile_2b.txt

### B2c - Check for clones
Calculate allelic similarity between all pairs of individuals from unfiltered `.vcf` file to assess potential clones [using [vcf_clone_detect.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_clone_detect.py) script]:

	$ vcf_clone_detect.py -v afra_2b.vcf -p popfile_2b.txt
	###1 - Pairwise comparisons of all individuals
	5356 comparisons completed
	
	###2 - Histogram (of pairwise genetic similarities)
	 69       1 *
	 70       0
	 71      16 ****************
	 72      23 ***********************
	 73      51 ***************************************************
	 74      62 **************************************************************
	 75     145 *********************************************************************#
	 76     331 *********************************************************************#
	 77     523 *********************************************************************#
	 78     743 *********************************************************************#
	 79     963 *********************************************************************#
	 80    1026 *********************************************************************#
	 81     752 *********************************************************************#
	 82     467 *********************************************************************#
	 83     174 *********************************************************************#
	 84      45 *********************************************
	 85      12 ************
	 86       7 *******
	 87       1 *
	 88       1 *
	 89       4 ****
	 90       0
	 91       0
	 92       0
	 93       0
	 94       1 *
	 95       2 **
	 96       0
	 97       5 *****
	 98       1 *
	 99       0
	
	###3 - List of highest matches
	98.02	0	[JD]	AFMJD6868H vs AFMJD6856H	3140.5/3204	4262	4065
	97.51	0.51	[JD]	AFMJD6856H vs AFMJD6854H	3157.5/3238	4644	3717
	97.36	0.15	[GD]	AFMGD6808H vs AFMGD6823H	2622.0/2693	4580	3236
	97.34	0.02	[JD]	AFMJD6868H vs AFMJD6854H	2924.0/3004	4410	3717
	97.29	0.05	[GD]	AFMGD6808H vs AFMGD6801H	3178.5/3267	4253	4137
	97.14	0.15	[GD]	AFMGD6823H vs AFMGD6801H	2789.0/2871	3857	4137
	95.79	1.35	[PS]	AFMPS6619H vs AFMPS6618H	3253.0/3396	4486	4033
	95.39	0.4	[PX]	AFMPX6981H vs AFMPX6982H	2361.0/2475	3600	3998
	94.28	1.11	[JD]	AFMJD6857H vs AFMJD6871H	2586.0/2743	4461	3405
	94.0	-------------------- Potential threshold --------------------
	89.47	4.81	[GD-PD]	AFMPD6989H vs AFMGD6799H	3173.5/3547	4506	4164
	89.31	0.16	[JD]	AFMJD6856H vs AFMJD6878H	2801.5/3137	4590	3670
	89.24	0.07	[JD]	AFMJD6855H vs AFMJD6861H	3034.0/3400	4506	4017
	89.22	0.02	[JD]	AFMJD6854H vs AFMJD6878H	2574.0/2885	4338	3670
	88.4	0.82	[JD]	AFMJD6868H vs AFMJD6878H	2529.0/2861	4314	3670
	87.16	1.24	[JD]	AFMJD6855H vs AFMJD6877H	2641.0/3030	4590	3563
	86.75	0.41	[JD]	AFMJD6855H vs AFMJD6876H	2962.5/3415	4460	4078
	86.48	0.27	[PD]	AFMPD6998H vs AFMPD6991H	2605.5/3013	4131	4005
	86.35	0.13	[GD-PD]	AFMPD6991H vs AFMGD6803H	2508.5/2905	4563	3465
	86.32	0.03	[JD]	AFMJD6876H vs AFMJD6861H	2934.0/3399	4505	4017
	86.07	0.25	[WD]	AFMWD6935H vs AFMWD6924H	2410.0/2800	3995	3928
	86.04	0.03	[JD]	AFMJD6877H vs AFMJD6861H	2582.0/3001	4107	4017
	86.04	0.0	[JD-PD]	AFMJD6858H vs AFMPD6998H	2528.0/2938	4397	3664
	85.98	0.06	[JS]	AFMJS6682H vs AFMJS6678H	2461.5/2863	4363	3623
	85.86	0.12	[GD-PD]	AFMPD6988H vs AFMGD6803H	2596.5/3024	4682	3465
	85.79	0.07	[GD-PD]	AFMPD6990H vs AFMGD6803H	2266.5/2642	4300	3465
	85.61	0.18	[GD]	AFMGD6804H vs AFMGD6822H	2023.0/2363	4236	3250
	85.49	0.12	[JD]	AFMJD6877H vs AFMJD6854H	2415.0/2825	4231	3717
	85.46	0.03	[JS]	AFMJS6678H vs AFMJS6687H	2174.0/2544	4418	3249
	85.45	0.01	[GS]	AFMGS6736H vs AFMGS6756H	2473.0/2894	4650	3367
	85.43	0.02	[JD]	AFMJD6877H vs AFMJD6876H	2550.0/2985	4030	4078
	85.4	0.03	[JD-PD]	AFMJD6858H vs AFMPD6991H	2688.5/3148	4266	4005
	85.4	0.0	[JD]	AFMJD6856H vs AFMJD6877H	2580.0/3021	4581	3563
	85.34	0.06	[JD-PD]	AFMJD6877H vs AFMPD6991H	2526.0/2960	4078	4005
	85.04	0.3	[GD]	AFMGD6819H vs AFMGD6803H	2456.0/2888	4546	3465

There appear to be 4 sets of clones (with 94-98% match), always within a single population. Remainder of pairwise comparisons are <90% match. The `PX` clone group represents a small recruit (<1.5 cm; `AFMPX6982H`) sampled directly adjacent to an A. fragilis colony collected from a depth of 67 m depth (`AFMPX6981H`).

	###4 - Clonal groups (threshold: 94.0)
	JD: AFMJD6868H, AFMJD6856H, AFMJD6854H (97.51-98.02 %)
	GD: AFMGD6823H, AFMGD6808H, AFMGD6801H (97.29-97.36 %)
	PS: AFMPS6619H, AFMPS6618H (95.79 %)
	JD: AFMJD6857H, AFMJD6871H (94.28 %)
	PX: AFMPX6982H, AFMPX6981H (95.39 %)

Eliminate samples from each clone set so that only one individual of clonal group remains per population:

	$ echo "AFMPS6619H\nAFMJD6854H\nAFMJD6868H\nAFMPX6981H\nAFMGD6808H\nAFMGD6823H\nAFMJD6871H" > samples_to_remove.txt
	$ vcftools --vcf afra_2b.vcf --remove samples_to_remove.txt --mac 1 --min-alleles 2 --recode --stdout > afra_2c.vcf
	$ tail -n 4 out.log

Create popfile (`popfile_noclones_2c.txt`) for remaining samples (coded PX individual AFMPX6982H as belonging to population PD) [using [popfile_from_vcf.py](https://github.com/pimbongaerts/radseq/blob/master/popfile_from_vcf.py) script]:
	
	$ popfile_from_vcf.py afra_2c.vcf 4 5 | sed 's/ PX/     PD/g' > popfile_test.txt

### B2d - Minimum representation filter
Eliminate SNPs in `.vcf` that are genotyped for less than 50% of individuals in each population (rather than using an overall missing-data filter) [using [vcf_minrep_filter.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_minrep_filter.py) script]:

	$ vcf_minrep_filter.py afra_2c.vcf popfile_noclones_2c.txt 0.5 afra_2d.vcf
	1814 out of 5120 SNPs failing 0.5 threshold
	Lowest overall proportion genotyped of remaining SNPs: 0.5463917525773195
	
### B2e - Check for deviations from HWE
Convert `.vcf` file to Arlequin format [using [vcf_spider.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_spider.py) script]:
	
	$ vcf_spider.py afra_2d.vcf popfile_noclones_2c.txt afra_2e.arp

Run HWE test for all 7 populations using `arlecore` (exact test using a Markov chain with chain length = 1 000 000 and dememorization steps = 100 000):

	$ arlecore afra_2e.arp arl_run.ars

Scan output for 'suspicious' SNPs, that meet one of these two criteria :

* Significantly deviating from HWE (excess or deficit) in min. 5 out of 7 populations. No Bonferroni/FDR applied as using multiple-population criteria (observed in min. 5 pops), so currently no SNPs expected to be identified just by chance (3,306 SNPs * 0.05^5). This still allows SNPs to remain in the dataset that deviate from HWE due to being under selection (associated with a particular population or habitat; never > 4 pops)
* Exhibiting a significant excess of heterozygotes in at least 1 population (to filter out potential paralogous loci). NB: Although this may be be due to overdominant selection or outbreeding.

Scrape XML output for SNPs meeting above criteria (using P-value cut-off of 0.05; `.vcf` file included to obtain original CHROM/POS identifiers)  [using [E - Scripts/bermuda_arlequin_hwe_summary.py](../E\ -\ Scripts/bermuda_arlequin_hwe_summary.py) script]:
	
	$ python3 bermuda_arlequin_hwe_summary.py -i afra_2e.xml -p 0.05 -d 5 -v afra_2d.vcf
	Total of 3306 SNPs evaluated
	0 SNPs monomorphic across all populations
	1926 SNPs with a significant deficit in heterozygosity in at least one population
	
	Suspicious SNPs/loci (673):
	338 SNPs with a significant excess in heterozygosity in at least one population
	335 addditional SNPs significantly deviating from HWE in 5 or more populations

Remove loci/CHROMs (in `afra_2e_elimchrom_uniq.txt`) and individual SNPs (in `afra_2e_elimpos.txt`) from `.vcf` file [using [vcf_remove_chrom.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_remove_chrom.py) script]:

	$ cat afra_2e_elimchrom.txt| sort | uniq > afra_2e_elimchrom_uniq.txt # reduce to only uniq loci
	$ wc -l afra_2e_elimchrom_uniq.txt afra_2e_elimpos.txt
     188 afra_2e_elimchrom_uniq.txt
     335 afra_2e_elimpos.txt
	$ vcf_remove_chrom.py afra_2d.vcf afra_2e_elimchrom_uniq.txt > afra_2e_temp1.vcf
	$ vcftools --vcf afra_2e_temp1.vcf --exclude-positions afra_2e_elimpos.txt --recode --stdout > afra_2e.vcf

### B2f - Final datasets
Dataset no-clones, bi-allelic `afra_2f.vcf` (50% max. missing data per pop; reduced to bi-allelic as only 5 multi-allelic SNPs):
	
	$ vcftools --vcf afra_2e.vcf --min-alleles 2 --max-alleles 2 --recode --stdout > afra_2f.vcf
	$ tail -n 4 out.log
	After filtering, kept 97 out of 97 Individuals
	Outputting VCF file...
	After filtering, kept 2568 out of a possible 2573 Sites
	Run Time = 0.00 seconds

Dataset without western location (WD), no-clones, biallelic `afra_noWD_2f.vcf`:

	$ vcftools --vcf afra_2f.vcf --keep popfile_noWD_2f.txt --recode --stdout > afra_noWD_2f.vcf
	$ tail -n 4 out.log
	After filtering, kept 83 out of 97 Individuals
	Outputting VCF file...
	After filtering, kept 2568 out of a possible 2568 Sites
	Run Time = 0.00 seconds
	
Dataset with-clones, bi-allelic `afra_wclon_2f.vcf`:

	$ tail -n +13 afra_2e.vcf | cut -f 1,2 > filtered_snps.txt
	$ vcftools --vcf afra_2b.vcf --positions filtered_snps.txt --min-alleles 2 --max-alleles 2 --keep popfile_wclon_2f.txt --recode --stdout > afra_wclon_2f.vcf
	$ tail -n 4 out.log
	After filtering, kept 104 out of 104 Individuals
	Outputting VCF file...
	After filtering, kept 2568 out of a possible 5123 Sites
	Run Time = 0.00 seconds

### B2g - Sample performance
Output number of genotyped and missing SNPs for each sample [using [vcf_missing_data.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_missing_data.py) script]:

	$ vcf_missing_data.py afra_wclon_2f.vcf > afra_wclon_2f_stats.txt

## B3 - Outlier analyses
### B3a - Lositan
Convert no-clones dataset `afra_2f.vcf` to Lositan/Genepop format [using [vcf_spider.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_spider.py) script]:

	$ vcf_spider.py afra_2f.vcf popfile_noclones_2c.txt afra_3a.lositan

Run the `Lositan` GUI (on Linux machine due to unresolved OSX bug) to assess outliers for the following population datasets:

|filename|description|populations|
|---|---|---|
|`fdist_shallow.txt`|within shallow between locations|PS, JS, GS|
|`fdist_deep.txt`|within deep between locations|PD, JD, GD|
|`fdist_PS-PD.txt`|within Princess between depths|PS, PD|
|`fdist_JS-JD.txt`|within JohnSmith between depths|JS, JD|
|`fdist_GS-GD.txt`|within Gurnet between depths|GS, GD|
|`fdist_overall.txt`|between all populations|PS, PD, JS, JD, GS, GD, WD|

### B3b - BayeScan
Convert no-clones dataset `afra_2f.vcf` to BayeScan format [using [vcf_spider.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_spider.py) script]:

	$ vcf_spider.py afra_2f.vcf popfile_noclones_2c.txt afra_3b.bayescan
	
Run `bayescan` on the overall dataset:

	$ bayescan afra_3b.bayescan -threads 28
	$ cp afra_3b.baye_fst.txt bayescan_overall.txt

### B3c - Summary
Summarise outlier analyses in a single table with boolean values (`TRUE` = significant outlier) [using [E - Scripts/bermuda_outliers_summary.py](../E\ -\ Scripts/bermuda_outliers_summary.py) script]:

	$ python3 bermuda_outliers_summary.py afra_2f.vcf outlier_files/fdist_overall.txt outlier_files 0.01 0.99 > outlier_summary.txt

Classify outliers types [using [E - Scripts/bermuda_outliers_classify.py](../E\ -\ Scripts/bermuda_outliers_classify.py) script]:
	
	$ python3 bermuda_outliers_classify.py outlier_summary.txt > outlier_summary_types.txt
	$ tail -n +2 outlier_summary_types.txt | cut -f 13 | sort | uniq -c | sort -n
	  25 DEPTH_OUTLIER
	  31 FDIST_AND_BAYESCAN_OUTLIER
	 119 FDIST_OUTLIER
	2393 NON_OUTLIER

Depth-associated outliers are identified as those SNPs that meet all of these 3 criteria:

* identified as overall outlier by both `fdist` AND `bayescan`
* identified as outliers in at least 2 individual between-depth comparisons
* never identified in any of the between-location (within-depth) comparisons

Create file with all depth-outlier SNPs:

	$ grep "DEPTH_OUTLIER" outlier_summary_types.txt | cut -f 2,3 > depth_outlier_snps.txt

Extract one representative depth-outlier SNP per CHROM from summary file (after assessing that SNPs on same CHROM show the same pattern):

	$ grep "DEPTH_OUTLIER" outlier_summary_types.txt | cut -f 2,3 | sort -u -k1,1 > depth_outliers.txt
	$ wc -l depth_outliers.txt
	      12 depth_outliers.txt

Output genotype frequencies of depth-outliers [using [vcf_genotype_freqs.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_genotype_freqs.py) script]:

	$ vcf_genotype_freqs.py afra_2f.vcf factors_afra_3e.txt depth_outliers.txt > afra_b3c_freqs.txt

Extract Fdist & BayeScan outliers (that are not depth-outliers) and output genotype frequencies (for supplementary materials) [using [vcf_genotype_freqs.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_genotype_freqs.py) script]:

	$ grep "FDIST_AND_BAYESCAN_OUTLIER" outlier_summary_types.txt | cut -f 2,3 > other_overall_outliers.txt
	$ vcf_genotype_freqs.py afra_2f.vcf factors_afra_3e.txt other_overall_outliers.txt > afra_b3c_other_freqs.txt
	
### B3d - Non-outlier dataset
Summarise outlier analyses but with less stringent threshold, into a single table with boolean values (`TRUE` = significant outlier) [using [E - Scripts/bermuda_outliers_summary.py](../E\ -\ Scripts/bermuda_outliers_summary.py) script and [E - Scripts/bermuda_outliers_classify.py](../E\ -\ Scripts/bermuda_outliers_classify.py) script]:

	$ python3 bermuda_outliers_summary.py afra_2f.vcf outlier_files/fdist_overall.txt outlier_files 0.05 0.95 > outlier_summary.txt
	$ python3 bermuda_outliers_classify.py outlier_summary.txt > outlier_summary_types.txt
	$ tail -n +2 outlier_summary_types.txt | cut -f 13 | sort | uniq -c | sort -n
	  26 DEPTH_OUTLIER
	  43 FDIST_AND_BAYESCAN_OUTLIER
	 230 FDIST_OUTLIER
	2269 NON_OUTLIER

Extract non-outlier SNPs (by eliminating all SNPs that were an overall outlier for either/both `fdist` and `bayescan`) and create bi-allelic dataset without-clones and outlier SNPs: `afra_neut_3d.vcf`:

	$ grep "NON_OUTLIER" outlier_summary_types.txt | cut -f 2,3 > neutral_snps.txt
	$ vcftools --vcf afra_2f.vcf --positions neutral_snps.txt --recode --stdout > afra_neut_3d.vcf
	$ tail -n 4 out.log
	After filtering, kept 97 out of 97 Individuals
	Outputting VCF file...
	After filtering, kept 2269 out of a possible 2568 Sites
	Run Time = 0.00 seconds

### B3e - Outlier annotation
Create FASTA file with outlier sequences (using `sed` to remove gaps) [using [fasta_include.py](https://github.com/pimbongaerts/radseq/blob/master/fasta_include.py) script]:

	$ cut -f 1 depth_outliers.txt | sort > depth_outlier_chroms.txt
	$ fasta_include.py afra_1c.fa depth_outlier_chroms.txt | sed 's/-//g' > depth_outliers.fa

Blast outlier sequencess (using `blastn`) against the NCBI non-redundant database:

	$ OUTPUT_FORMAT="7 qseqid sseqid length nident pident evalue bitscore staxids stitle"
	$ blastn -query depth_outliers.fa -remote -db nt -task blastn -evalue 0.0001 -max_target_seqs 1 -outfmt $OUTPUT_FORMAT -out afra_depth_outliers_blastn_nt.txt
	$ grep -v "#" afra_depth_outliers_blastn_nt.txt | cut -f 1,6,9
	4482	1e-11	PREDICTED: Acropora digitifera phospholipase D1-like (LOC107352864), mRNA
	8374	5e-11	PREDICTED: Acropora digitifera uncharacterized LOC107332840 (LOC107332840), mRNA

## B4 - Genetic structure
### B4a - STRUCTURE
##### Overall dataset
Modify STRUCTURE's `mainparams` file (manually) and enter the number of unique CHROMs and total number of individuals:

	$ tail -n +13 afra_2f.vcf | cut -f 1 | sort | uniq | wc -l
	    1579	# number of unique CHROMs

Run overall dataset (`afra_2f.vcf`) through STRUCTURE (for K = 2..7) [using the multi-threading wrapper [structure_mp](https://github.com/pimbongaerts/radseq/blob/master/structure_mp.py)]:

	$ structure_mp.py afra_2f.vcf popfile_noclones_2c.txt 7 10 10
	Initialise individuals and populations...
	Subsample SNPs (one random SNP per locus)... [1579 SNPs/loci]
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
	K = 4: MedMeaK 3.0 MaxMeaK 3 MedMedK 3.0 MaxMedK 3
	K = 5: MedMeaK 3.0 MaxMeaK 4 MedMedK 3.0 MaxMedK 4
	K = 6: MedMeaK 3.0 MaxMeaK 4 MedMedK 3.0 MaxMedK 4
	K = 7: MedMeaK 2.5 MaxMeaK 3 MedMedK 2.5 MaxMedK 3

##### "Neutral" dataset
Modify STRUCTURE's `mainparams` file (manually) and enter the number of unique CHROMs and total number of individuals:

	$ tail -n +13 afra_neut_3d.vcf | cut -f 1 | sort | uniq | wc -l
	    1395	# number of unique CHROMs

Run the non-outlier ("neutral") dataset (`afra_neut_3d.vcf`) through STRUCTURE [using the multi-threading wrapper [structure_mp](https://github.com/pimbongaerts/radseq/blob/master/structure_mp.py)]:

	$ structure_mp.py afra_neut_3d.vcf popfile_noclones_2c.txt 7 10 10
	Initialise individuals and populations...
	Subsample SNPs (one random SNP per locus)... [1395 SNPs/loci]
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
	K = 4: MedMeaK 3.0 MaxMeaK 4 MedMedK 3.0 MaxMedK 4
	K = 5: MedMeaK 2.0 MaxMeaK 4 MedMedK 2.0 MaxMedK 4
	K = 6: MedMeaK 2.0 MaxMeaK 3 MedMedK 2.0 MaxMedK 3
	K = 7: MedMeaK 1.0 MaxMeaK 2 MedMedK 1.0 MaxMedK 2

### B4b - PCA
Run a PCA on the overall dataset (`afra_2f.vcf`) using the `adegenet` package in R [detailed in [E - Scripts/bermuda_vcf2pca.R](../E\ -\ Scripts/bermuda_vcf2pca.R) script]:

	
	$ Rscript bermuda_vcf2pca.R afra_2f.vcf popfile_pca_4b.csv

### B4c - GD matrix
Calculate GD matrix (Hamming-distance) between individual samples (using dataset that includes clones `afra_wclon_2f.vcf`) [using [vcf_gdmatrix.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_gdmatrix.py) script]:

	$ vcf_gdmatrix.py afra_wclon_2f.vcf popfile_wclon_ordered_4c.txt > afra_wclon_gd_4c.txt

### B4d - DAPC
Run a DAPC on the overall dataset (`afra_2f.vcf`) using the `adegenet` package in R [detailed in [E - Scripts/bermuda_vcf2dapc.R](../E\ -\ Scripts/bermuda_vcf2dapc.R) script]:
	
	$ Rscript bermuda_vcf2dapc.R afra_2f.vcf popfile_pca_4b.csv.txt

## B5 - Additional analyses 
### B5a - Differentiation and diversity stats
Calculate several statistics regarding differentiation in R using `hierfstat` and `adegenet` packages (as well as `vcfR` for importing `.vcf` file) [detailed in [E - Scripts/bermuda_diff_stats.R](../E\ -\ Scripts/bermuda_diff_stats.R) script]:

	* General F-statistics and measures of heterozygosity
	* G-statistic test for overall population structure
	* Non-parametric Kruskal-Wallis tests to assess for an overall difference between groups of populations
	* Non-parametric Wilcoxon signed rank tests to assess for differences between groups of populations


	$ Rscript bermuda_diff_stats.R afra_2f.vcf popfile_pca_4b.csv "#bae4b3" "#238b45" "#74c476"
	Processed variant: 2635

	Weir and Cockerham Fst and Fis:
	$FST
	[1] 0.06998256	
	$FIS
	[1] 0.488847

	Expected heterozygosity (Hs) for each pop:
	       GD        GS        JD        JS        PD        PS        WD
	0.1672369 0.1615541 0.1658039 0.1618625 0.1781947 0.1641213 0.2103866
	
	G-statistic test for overall structure:
	Monte-Carlo test
	Observation: 56314.52	
	Based on 99 replicates
	Simulated p-value: 0.01
	Alternative hypothesis: greater
	     Std.Obs  Expectation     Variance
	    20.97319  39710.50866 626753.48774
	
	G-statistic test for overall structure - eastern pops only:
	Monte-Carlo test
	Observation: 43589.48	
	Based on 99 replicates
	Simulated p-value: 0.01
	Alternative hypothesis: greater
	     Std.Obs  Expectation     Variance
	    15.24963  32330.52455 545101.69050
	
	Pairwise Fst when averaged over all loci (without NA):
	     Habt.G      Habt.J      Habt.P     Shal.GJ     Shal.JP     Shal.GP
	0.044362993 0.058421089 0.044208650 0.003175759 0.012652416 0.013741802
	    Deep.GJ     Deep.JP     Deep.GP     West.GW     West.JW     West.PW
	0.011968014 0.011920970 0.001320425 0.041809349 0.041948989 0.034801295
	
	Pairwise Fst:
	           1          2          3          4          5          6
	2 0.08013806
	3 0.04999698 0.09408470
	4 0.07687531 0.04244536 0.09540889
	5 0.03295640 0.06640810 0.04102523 0.06377280
	6 0.08750017 0.05033273 0.10496044 0.04984853 0.07136143
	7 0.06846539 0.08666910 0.06682729 0.08322243 0.05380059 0.09108684
	
	Non-parametric Kruskal-Wallis to test for overall differences
	Kruskal-Wallis rank sum test
	Kruskal-Wallis chi-squared = 499.43, df = 3, p-value < 2.2e-16
	
	Non-parametric Wilcoxon: between-depths vs within-shallow
	Wilcoxon rank sum test with continuity correction
	W = 20749000, p-value < 2.2e-16
	alternative hypothesis: true location shift is not equal to 0
	
	Non-parametric Wilcoxon: between-depths vs within-deep
	Wilcoxon rank sum test with continuity correction
	W = 21537000, p-value < 2.2e-16
	alternative hypothesis: true location shift is not equal to 0
	
	Non-parametric Wilcoxon: within-shallow vs within-deep
	Wilcoxon rank sum test with continuity correction
	W = 17487000, p-value = 0.5338
	alternative hypothesis: true location shift is not equal to 0
	
	Non-parametric Wilcoxon: between-westeast vs within-deep
	Wilcoxon rank sum test with continuity correction
	W = 17196000, p-value < 2.2e-16
	alternative hypothesis: true location shift is not equal to 0
	
	Non-parametric Wilcoxon: between-depths vs between-westeast
	Wilcoxon rank sum test with continuity correction
	W = 21257000, p-value = 0.4399
	alternative hypothesis: true location shift is not equal to 0


### B5b - Admixture
Visualising patterns of admixture between shallow and deep individuals for the most divergent SNPs. Most divergent SNPs were selected by:

* Rerunning STRUCTURE for only the eastern populations (two clusters)
* Assigning individuals with an averaged assignment of >0.98 as "parental groups", and the remaining individuals as "admixed group"
* Including only loci with a frequency difference of 0.5 between the two parental groups (corresponding to shallow and deep).

Modify STRUCTURE's `mainparams` file (manually) and enter the number of unique CHROMs and total number of individuals:

	$ tail -n +13 afra_noWD_2f.vcf | cut -f 1 | sort | uniq | wc -l
	1579 # number of unique CHROMs

Run the dataset with eastern populations (`afra_noWD_2f.vcf`) through STRUCTURE [structure_mp](https://github.com/pimbongaerts/radseq/blob/master/structure_mp.py)]:

	$ structure_mp.py afra_noWD_2f.vcf popfile_noWD_2f.txt 2 10 10
	Initialise individuals and populations...
	Subsample SNPs (one random SNP per locus)... [1579 SNPs/loci]
	Outputting 10 STRUCTURE files...10 reps DONE
	Executing 10 parallel STRUCTURE runs for K = 2 ...10 reps DONE
	Running CLUMPP on replicates for K = 2 ...

Sorted CLUMPP output by assignment:

	$ cp clumpp_K2.out.csv clumpp_K2_noWD.csv
	$ cat clumpp_K2_noWD.csv | sort -t ',' -k3n > afra_assign.txt

Manually rearranged one shallow individual that had slightly higher assignment (otherwise sorted individuals cluster by depth of origin).

Create a genotype matrix for loci that were identified as outliers (outputted first) or have an allele freq of >=0.5 between parental indivs (assignment >0.98). Reference and alternative allele are switched if necessary (so that reference allele corresponds to deep cluster) [using [vcf_ancestry_matrix.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_ancestry_matrix.py) script]:

	$ vcf_ancestry_matrix.py afra_noWD_2f.vcf afra_assign.txt 0.98 0.5 --include depth_outlier_snps.txt > afra_matrix_5b.txt
	$ tail -n +2 afra_matrix_5b.txt | wc -l  # number of SNPs in matrix
     133

### B5c - Symbiont typing
Representative samples of *Agaricia fragilis* (n = 45) from the 7 populations were sequenced for the *Symbiodinium COX1* region. Only a single haplotype was recovered:

	>Agaricia_COX1_HAPA
	GCAACAATCTTTTCTCCATTCAAGATTATTCCAGAGAAGATAGCAATTATAGCTCCTAATGAAAGAACAAAA
	TGAAAATGTGCTACAACATAATATGTATCATGTAATCCTAGATCCACTGCACCATTTCCAAGAATTATTCCT
	GTTGACCCACCTATCGTAAACATTAATAAAAAGAGATGTGAGAGGAAGACAGAAGTAATTCTAAGGTGTAAT
	AATGGTGGATTTGAGAGATATGTAAAAAGCCAATTAAAGATTTTTGTACCAGTTGGTAAGGATATTAAGATT
	GTAACTCCTGTAAAATAAGCTCTTGTATCACTTTCTAAACCTACAGTATACATATGATGTCCCCAAACAAGA
	CCTCCAAGAAGAGAAATAGATGACATGGCAAAGATCATTGATTGGTTAGCAAAGATTATTAACTGTAAAATA
	CCAGAAATTATTATGGAAATGATCCCAAATGCAGGAATTATTAATATGTAAACTTCTGGATGTCCAAAAAAC
	CAAAATAAATGTTGATAGAATATAGGATCTCCTCCAAATATTGGATCAAAGAAAAGTGTATTAGAATGAAGA
	TCACCCAATATTAAAAGAAGTGTACCAGATAAGATTGGTAATGTTAATAAAAGCATGAAAGCTGTAATCAAG
	AAAGCCCAAGGAAATAATGGGATAGAAGATAATATCAGATAATAAGATCTCAGAAAATGAATTGTTGTCCAA
	AAGTTAAGAGATGTAAGACATGAAGATATACCAGAGATTATTAATCCAAATATAAGATTTCCTGTACTTGAA
	GGTGATAAAGTCA

In addition, 3 RAD loci (that were also recoved in *Stephanocoenia intersepta*) matched *Symbiodinium* chloroplast minicircle sequences in the `blastn` search conducted in step **B2a - Remove Symbiodinium and other contamination**. 

Extract those loci that are in common from the PyRAD `.loci` files and trim to 75bp [using [pyrad_filter.py](https://github.com/pimbongaerts/radseq/blob/master/pyrad_filter.py) script and [pyrad_trim.py](https://github.com/pimbongaerts/radseq/blob/master/pyrad_trim.py) script]:

	$ echo "2199\n9878\n11373" > afra_sym_chloroplast_loci.txt
	$ pyrad_filter.py storage_server_path/afra/outfiles/afra.loci afra_sym_chloroplast_loci.txt 0 0 > afra_sym_chloroplast.loci
	$ pyrad_trim.py afra_sym_chloroplast.loci 75 > afra_sym_chloroplast_trim.loci

Visual assessment of `afra_sym_chloroplast_trim.loci` alignment confirms that no variable sites are present in these 3 loci. Mitochondrial and chloroplast haplotypes for each sample are listed in: `afra_symbionts.txt`.

### B5d - Morphometrics
Corallite diameter and density were measured for representative samples of *Agaricia fragilis* (n = 54) with 5 pseudoreplicate measurements per sample. Averaged values (across the 5 measurements) for each specimen are listed in `afra_morpho_data.csv`.