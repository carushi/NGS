# NGS
* Scripts and snippets for NGS 

## script
### read\_count\_samtools.py
* to calculate read counts mapped around tss (+/- 500) using samtools
* in the same way of featureCounts

```
python read_count_samtools.py [bam or sam file] [outputfile] refGene.txt
```

* corresponding command for featureCounts

```
featureCounts -s 1 -O -a tss.gtf -o [outputfile] [bam or sam file]
```

* <b>samtools\_commands\_for\_featureCounts</b>
* function
* arguments
	* bamfile
	* outputfile
	* reffile (refGene.txt)
	* flagopt ("-F 256 ")
		* "-F 256" means "excluding multi-mapped reads".
		* This flag is applied for both strands.
		* If you would like to include multi-mapped reads, please add "-M" option to featureCounts.
	* printMode (False)
		* If printMode is active, it outputs samtools command lines which are used to count up.
	* single-end mode specific


