# NGS
* scripts for NGS 

## script
* read\_count\_samtools.py
	* to compare the result of samtools and featureCounts
	* calculates read counts mapped around tss (+/- 500) in the same way of featureCounts.
	* python read\_count\_samtools.py bamfile bedfile reffile
	* If printMode is on, print samtools command lists used to count up.
	
