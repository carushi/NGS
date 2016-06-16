import os
import sys
import subprocess

# The format of end position is different between refseq and gtf
# gtf -> include end position
# refGene -> exclude end position
# samtools -> include end position
# sam, refseq.gtf -> 1-based?
# bam, refGene -> 0-based

def start_position(contents, dist = 500):
    if contents[3] == "+":
        return int(contents[4])-dist+1
    else:
        return int(contents[5])-dist+1

def end_position(contents, dist = 500):
    if contents[3] == "+":
        return int(contents[4])+dist+1
    else:
        return int(contents[5])+dist+1


def get_ref_list(reffile):
    dict = {}
    with open(reffile) as f:
        for line in f.readlines():
            contents = line.rstrip('\n').split('\t')
            if contents[1] not in dict.keys():
                dict[contents[1]] = []
            dict[contents[1]].append(line)
    return dict

def is_overlap(rstart, cigar, start, end):
    pre = 0
    for i in range(1, len(cigar)):
        c = cigar[i]
        if c.isalpha() or c == "=":
            if i == pre:
                pos = 1
            else:
                pos = int(cigar[pre:i])
            if c == "M" or c == "=":
                if start <= rstart+pos and rstart <= end: #include start and end
                    return True
                rstart += pos
            elif c == "D" or c == "N" or c == "X":
                rstart += pos
            else: #P or H or S or I
                pass
            pre = i+1
    return False

def uniq_read_count(output, start, end):
    wc = 0
    reads = str(output.decode("utf-8")).split('\n')
    for i in range(0, len(reads))[::-1]:
        if reads[i] in reads[0:i] or len(reads[i]) == 0:
            reads.pop(i)
    for line in reads:
        if line[0:2] == "[E":
            raise Exception(line)
        contents = line.split('\t')
        cigar = contents[5]
        if cigar[0:-1].isdigit():
            wc = wc+1
        else:
            if is_overlap(int(contents[3]), cigar, start, end):
                wc = wc+1
    return wc

def samtools_commands_for_featureCounts(bamfile, bedfile, reffile, flagopt = "-F 256 ", printMode = False):
    with open(bedfile, 'w')  as f:
        dict = get_ref_list(reffile)
        count = 0
        for gene in sorted(dict.keys()):
            output = b""
            f.write(gene+"\t")
            for line in dict[gene]:
                contents = line.rstrip('\n').split('\t')
                start = start_position(contents)
                end = end_position(contents)
                chr = contents[2]
                flag = ""
                if contents[3] == "+":
                    flag = "-F 16 "
                else:
                    flag = "-f 16 "
                cmd = "samtools view -m 1 "+flag+flagopt+" -t SP "+bamfile+" \""+chr+":"+str(start)+"-"+str(end)+"\" | grep NH:i:1 | grep -v NH:i:1[0-9] "
                if printMode:
                    print(cmd+"# "+gene)
                else:
                    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    output += process.communicate()[0]
            count = count+1
            if count%(len(dict.keys())/10) == 0:
                print("#", count%(len(dict.keys())/10), "completed...")
            if printMode:
                continue
            try:
                wc = uniq_read_count(output, start, end)
                f.write(str(wc)+"\n")
                f.flush()
            except:
                print(sys.exc_info()[0])
                raise

def convert_refseq_to_tss_gtf(reffile):
    with open(reffile) as f:
        for line in f.readlines():
            contents = line.rstrip('\n').split('\t')
            start = start_position(contents)
            end = end_position(contents)
            print("\t".join([contents[2], contents[1], "gene", str(start+1), str(end), ".", contents[3], ".", "gene_id "+contents[1]]))
            print("\t".join([contents[2], contents[1], "exon", str(start+1), str(end), ".", contents[3], ".", "gene_id "+contents[1]]))

if __name__ == '__main__':
    if len(sys.argv) <= 3:
        if len(sys.argv) == 2:
            convert_refseq_to_tss_gtf(sys.argv[1])
        elif len(sys.argv) == 1:
            convert_refseq_to_tss_gtf(sys.argv[1], int(sys.argv[2]))
    elif len(sys.argv) > 3:
        argvs = sys.argv
        samtools_commands_for_featureCounts(argvs[1], argvs[2], argvs[3])
    else:
        bamout = "temp.bam"
        bedfile = bamout.replace(".bam", ".bed")
        reffile = "../example/refGene.txt"
        samtools_commands_for_featureCounts(bamout, bedfile, reffile)
