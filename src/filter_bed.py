#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print ("""
    filter_bed -i <input-bed-file> -n <num_reads>
    
    Parameters:
        -i/--input-bed-file    [string  :    path to bed file with column 4 as ; seperated alignments                       ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./                                             ]
        -n--number-of-reads    [int     :    min number of reads to keep a bed record   Default: 2                          ]
    Version:                   1.0.0
         """)
#parameters
input_bed_file = ''
output_dir = ''
#is_rm_tmp=None
num_of_reads = 2

def setDefault():
    global output_dir
    output_dir = './'
    #global is_rm_tmp
    #is_rm_tmp=True
    global num_of_reads
    num_of_reads = 2

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:",["help",
                                                    "input-bed-file=",
                                                    "output-dir="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        #elif opt in ("-k","--keep-tmp"):
        #    global is_rm_tmp
        #    is_rm_tmp = False
        elif opt in ("-i", "--input-bed-file"):
            global input_bed_file
            input_bed_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-n", "--number-of-reads"):
            global num_of_reads
            num_of_reads = int(arg)

#def make_dir(path):
#    if not os.path.exists(path):
#        os.mkdir( path, 0o755 ) 

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

def run_convert():
    infile = "%s" % input_bed_file
    f=open(infile,"r")
    outfile2 = "%s" % output_dir+'/tmp/'+getOutputName2(input_bed_file,"f.bed")
    f2=open(outfile2, "w")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            tmp=line.split("\t")
            alignments=tmp[3].split(";")
            nms = {}
            for l in range(len(alignments)):
                nm=alignments[l].split(",")[0]
                if nm not in nms:
                    nms[nm]=1
            if len(nms) >= num_of_reads:
                f2.write(line)
    f.close()
    f2.close()

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_file=='':
        usage()
        exit(1);
    else:
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        run_convert()
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
