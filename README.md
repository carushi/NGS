# NGS
* Scripts and snippets for NGS 

## script
* read\_count\_samtools.py
	* <b>format: python read\_count\_samtools.py [bam or sam file] [outputfile] refGene.txt </b>
	* to compare the result between samtools and featureCounts
		* corresponding commands:
		* featureCounts -s 1 -O -M -a tss.gtf -o [outputfile] [bam or sam file]
	* It calculates read counts mapped around tss (+/- 500) in the same way of featureCounts.

	* If printMode is active, it outputs samtools command lines which are used to count up.
	
