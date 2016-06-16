cat file_list.txt | xargs -I{} echo "python read_count_samtools.py {} {}_sam.bed ../sample/refGene.txt &"  | bash
python read_count_samtools.py ../sample/refGene.txt > ../sample/tss.gtf
cat file_list.txt | xargs -I{} echo "featureCounts -T 5 -s 1 -O -a ../sample/tss.gtf -o {}.bed {} "  | bash
