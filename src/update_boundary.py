#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print """
    update_boundary -i <bedfile>
    
    Parameters:
        -i/--input-bed-file    [string  :    path to BED with reads                                          ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./                              ]
        -p/--percent           [float   :    percentage of reads at each side to consider                    ]
        -c/--cutoff            [float   :    increase cutoff to change boundary                              ]
        -k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep                        ]
    
    Version:                    1.0.0
          """
#parameters
input_bed_file = ''
output_dir = ''
is_rm_tmp=None
percent = 0.15
cutoff = 0.33

def setDefault():
    global output_dir
    output_dir = './'
    global percent
    percent = 0.1
    global cutoff
    cutoff = 0.33
    global is_rm_tmp
    is_rm_tmp=True

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hki:o:p:c:",["help",
                                                      "keep-tmp",
                                                      "input-bed-file=",
                                                      "output-dir=",
                                                      "percent=",
                                                      "cutoff="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        elif opt in ("-k","--keep-tmp"):
            global is_rm_tmp
            is_rm_tmp = False
        elif opt in ("-i", "--input-bed-file"):
            global input_bed_file
            input_bed_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-p", "--percent"):
            global percent
            percent=float(arg)
        elif opt in ("-c", "--cutoff"):
            global cutoff
            cutoff = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0755 ) 

def getFilename(filename):
    bn=basename(filename)
    tmp=bn.split(".")
    fn=tmp[0]
    for x in range(1, len(tmp)-1):
        fn=fn+"."+tmp[x]
    return fn

def getOutputName2(filename,append_format):
    fn=getFilename(filename)
    return fn+"."+append_format

def getNextPosition(p,sigar):
    if 'S' in sigar or 'H' in sigar: ##### not supporting clipped reads
        return -1
    else:
        s=''
        nums=[]
        num_str=''
        for x in range(len(sigar)):
            if sigar[x] in {'M','I','D','S','H'}:
                s=s+sigar[x]
                nums.append(int(num_str))
                num_str=''
            else:
                num_str=num_str+sigar[x]
        totalLen=0
        for x in range(len(s)):
            if s[x]=='M' or s[x]=='D':
                totalLen=totalLen+nums[x]
        return p+totalLen-1

def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

long_reads_sorted = []

def sort_reads(reads):
    long_reads = []
    for x in range(len(reads)):
        tmp=reads[x].split(",")
        start=tmp[3][1:]
        cigar=tmp[4]
        one_long = []
        one_long.append(reads[x])
        one_long.append(int(start))
        one_long.append(getNextPosition(int(start),cigar)-int(start)+1)
        long_reads.append(one_long)
    global long_reads_sorted
    long_reads_sorted=sorted(long_reads, key = lambda x: (x[1], x[2]))

def cal_bpl(num_reads,num_rm):
    long_reads=long_reads_sorted[num_rm:num_reads-num_rm]
    s=long_reads[0][1]
    e=long_reads[len(long_reads)-1][1]+long_reads[len(long_reads)-1][2]-1
    t=0
    for x in range(len(long_reads)):
        t=t+long_reads[x][2]
    return t*1.0/(e-s+1),s-1,e

def get_reads_str(num_reads,num_rm):
    long_reads=long_reads_sorted[num_rm:num_reads-num_rm]
    reads_str = ""
    for x in range(len(long_reads)):
        reads_str=reads_str+long_reads[x][0]
        if x!=len(long_reads)-1:
            reads_str=reads_str+";"
    return reads_str

def get_update_line(line):
    tmp=line.split("\t")
    chr0 = tmp[0]
    start0 = tmp[1]
    end0 = tmp[2]
    reads = tmp[3].rstrip('\n').split(";")
    num_reads = len(reads)
    num_rm = int(num_reads*percent)
    if num_rm > 0 and num_reads-2*num_rm > 0:
        sort_reads(reads)
        basePerLength0,start0,end0=cal_bpl(len(long_reads_sorted),0)
        basePerLength1,start1,end1=cal_bpl(num_reads,num_rm)
        if basePerLength1 > basePerLength0*(1.0+cutoff):
            reads_str = get_reads_str(num_reads,num_rm)
            return chr0+"\t"+str(start1)+"\t"+str(end1)+"\t"+reads_str+"\n"
        else:
            return line
    return line
    

def run_update():
    infile = "%s" % input_bed_file
    f=open(infile,"r")
    outfile2 = "%s" % output_dir+'/'+getOutputName2(input_bed_file,"update.bed")
    f2=open(outfile2, "w")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            f2.write(get_update_line(line))
    f.close()
    f2.close()

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_file=='':
        usage()
        exit(1);
    else:
        use_real_path()
        make_dir(output_dir)
        make_dir(output_dir+'/tmp')
        run_update()
        remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
