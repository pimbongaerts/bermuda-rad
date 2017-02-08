#A - *Symbiodinium*
*Contents:* Details the analysis of the sequence data for the FACS-isolated Symbiodinium that are used as subtraction datasets. *Notebook author:* Pim Bongaerts. 

*[Click here to go back to the overview](https://github.com/pimbongaerts/bermuda-rad/)*

* **[A1 - PyRAD clustering](#a1---pyrad-clustering)**
	* [A1a - Raw sequence data](#a1a---raw-sequence-data)
	* [A1b - PyRAD clustering](#a1b---pyrad-clustering)
* **[A2 - Blastn database](#a2---blastn-database)**
	* [A2a - Filter loci to FASTA](#a2a---filter-loci-to-fasta)
	* [A2b - Create blastn database](#a2b---create-blastn-database)

##A1 - PyRAD clustering
Note: location of raw sequence data is referred to as `storage_server_path/...`.

###A1a - Raw sequence data
Replace barcodes of `fastq` files by actual sample names using [fastq_barcodes2samplenames.py](https://github.com/pimbongaerts/radseq/blob/master/fastq_barcodes2samplenames.py) script] (already done in NCBI SRA):
	
	$ fastq_barcodes2samplenames.py storage_server_path/symbiont_ref storage_server_path/sym_barcodes.txt

Assess number of reads for all 15 samples [using [fastq_seqcount.py](https://github.com/pimbongaerts/radseq/blob/master/fastq_seqcount.py) script]:

	$ fast_seqcount.py storage_server_path/symbiont_ref symrad_seqcounts.txt
	mean: 975,155 reads
	min: 146,130 reads
	max: 2,189,741 reads

Details of the 15 sequenced *Symbiodinium* isolates:

|Filename|Host sample code|Host species|Location|Depth|Reads|
|---|---|---|---|---|---|
|SYAFDC6630.fastq.gz|DC6630|*Agaricia fragilis*|Bermuda|13 m|481,784|
|SYAFDC6883.fastq.gz|DC6883|*Agaricia fragilis*|Bermuda|40 m|319,520|
|SYAGDC6262.fastq.gz|DC6262|*Agaricia grahamae*|Bonaire|50 m|359,608|
|SYAGDC6262.fastq.gz|DX6262|*Agaricia grahamae* (replicate)|Bonaire|50 m|3,417,809|
|SYLGDC0780.fastq.gz|DC0780|*Leptoseris glabra*|GBR - Myrmidon|40 m|217,133|
|SYLSDC1762.fastq.gz|DC1762|*Leptoseris scabra*|Coral Sea- Osprey|10 m|450,500|
|SYLGDC2773.fastq.gz|DC2773|*Leptoseris glabra*|GBR - Tijou|10 m|147,1451|
|SYLSDC2360.fastq.gz|DC2360|*Leptoseris scabra*|Coral Sea - Osprey|10 m|1,507,896|
|SYLEDC3380.fastq.gz|DC3380|*Leptoseris cf. scabra*|GBR - Yonge|40 m|2,189,741|
|SYMMXX0007.fastq.gz|XX0007|*Madracis mirabilis*|Curacao|15 m|1,146,570|
|SYMPDC4941.fastq.gz|DC4941|*Madracis pharensis*|Curacao|15 m|1,927,412|
|SYMPDC5337.fastq.gz|DC5337|*Madracis pharensis*|Curacao|50 m|2,013,450|
|SYMPDC5343.fastq.gz|DC5343|*Madracis pharensis*|Curacao|50 m|1,049,833|
|SYSIDC6721.fastq.gz|DC6721|*Stephanocoenia intersepta*|Bermuda|13 m|146,130|
|SYSIDC6911.fastq.gz|DC6911|*Stephanocoenia intersepta*|Bermuda|40 m|371,153|

###A1b - PyRAD clustering
The first two steps of PyRAD are ran to filter reads, and discard those with >5 sites with a PHRED quality <20. The third step is ran to cluster within samples, so that in the next step singleton reads can be excluded.

	$ pyrad -p params.txt -s123

Params file (only modified parameters shown; parameter 6 can be left empty for nextRAD data):

	$ head -n 15 params.txt
	==** parameter inputs for pyRAD version 3.0.6  **======================== affected step ==
	./                        ## 1. Working directory                                 (all)
	                          ## 2. Loc. of non-demultiplexed files (if not line 16)  (s1)
	                          ## 3. Loc. of barcode file (if not line 16)             (s1)
	vsearch-1.1.3-linux-x86_64  ## 4. command (or path) to call vsearch (or usearch)    (s3,s6)
	muscle                    ## 5. command (or path) to call muscle                  (s3,s7)
	GTGTAGAGG                 ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)
	14                        ## 7. N processors (parallel)                           (all)
	5                         ## 9. NQual: max # sites with qual < 20 (or see line 20)(s2)
	.85                       ## 10. Wclust: clustering threshold as a decimal        (s3,s6)
	rad                       ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)


##A2 - Blastn database
###A2a - Filter loci to FASTA
Filtering the `.clustS.gz` output files created in **A1b - PyRAD clustering** to output a representative sequence (first one) for each cluster that had a minimum of 2 representatives [using [pyradclust2fasta.py](https://github.com/pimbongaerts/radseq/blob/master/pyradclust2fasta.py) script]. This excludes singletons: min. one unique sequence in cluster with size >= 2 or min. two distinct sequences in cluster with each size = 1.

	$ pyradclust2fasta.py storage_server_path/symbiont_ref/clust.85 2 symbiont_ref_2a.fasta
	SYAFDC6630	35,981	164,882
	SYAFDC6883	26,256	132,798
	SYLSDC2360	69,636	224,980
	SYMMXX0007	71,028	224,884
	SYMPDC4941	83,075	261,618
	SYMPDC5337	54,505	148,458
	SYMPDC5343	36,762	108,664
	SYSIDC6721	16,046	61,788
	SYSIDC6911	28,045	103,594
	SYAGDC6262	13,564	44,434
	SYAGDX6262	73,776	226,108
	SYLEDC3380	71,945	208,672
	SYLGDC0780	19,141	114,426
	SYLGDC2773	64,624	210,244
	SYLSDC1762	32,762	128,210
	Total seqs in symbiont_ref.fasta: 697,146
	
Created a separate fasta just for *A. fragilis* and *S. intersepta* [using [pyradclust2fasta.py](https://github.com/pimbongaerts/radseq/blob/master/pyradclust2fasta.py) script]:

	$ mkdir storage_server_path/symbiont_ref/clust.85/afra_sint_only
	$ cp storage_server_path/symbiont_ref/clust.85/SYAFDC6630.clustS.gz /
		 storage_server_path/symbiont_ref/clust.85/SYAFDC6883.clustS.gz /
		 storage_server_path/symbiont_ref/clust.85/SYSIDC6721.clustS.gz /
		 storage_server_path/symbiont_ref/clust.85/SYSIDC6911.clustS.gz /
		 storage_server_path/symbiont_ref/clust.85/afra_sint_only
	$ python3 pyradclust2fasta.py storage_server_path/symbiont_ref/clust.85/afra_sint_only] 2 symbiont_ref_afra_sint.fasta
	SYAFDC6630	35,981	164,882
	SYAFDC6883	26,256	132,798
	SYSIDC6721	16,046	61,788
	SYSIDC6911	28,045	103,594
	Total seqs in symbiont_ref_afra_sint.fasta: 106,328

###A2b - Create blastn database
Remove gaps from FASTA before indexing:

	$ sed 's/-//g' symbiont_ref_2a.fasta > symbiont_ref_nogaps_2b.fasta

Create blastn database using FASTA file:

	./makeblastdb -in symbiont_ref_nogaps_2b.fasta -title SymRAD16 -dbtype 'nucl' -out SymRAD16