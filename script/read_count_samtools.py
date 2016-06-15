import os
import sys
import subprocess


def start_position(contents, dist = 500):
    if contents[3] == "+":
        return int(contents[4])-dist
    else:
        return int(contents[5])-dist

def end_position(contents, dist = 500):
    if contents[3] == "+":
        return int(contents[4])+dist
    else:
        return int(contents[5])+dist

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
            pos = int(cigar[pre:i])
            if c == "M" or c == "=":
                if start < rstart+pos and rstart < end:
                    return True
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
        for gene in sorted(dict.keys()):
            output = b""
            for line in dict[gene]:
                contents = line.rstrip('\n').split('\t')
                start = start_position(contents)
                end = end_position(contents)
                chr = contents[2]
                f.write(gene+"\t")
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
            if printMode:   continue
            try:
                wc = uniq_read_count(output, start, end)
                f.write(str(wc)+"\n")
                f.flush()
            except:
                print(sys.exc_info()[0])
                raise

if __name__ == '__main__':
    if len(sys.argv) > 3:
        argvs = sys.argv
        samtools_commands_for_featureCounts(argvs[1], argvs[2], argvs[3])
    else:
        bamout = "temp.bam"
        bedfile = bamout.replace(".bam", ".bed")
        reffile = "refGene.txt"
        samtools_commands_for_featureCounts(bamout, bedfile, reffile)
