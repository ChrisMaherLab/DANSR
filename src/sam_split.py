#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print ("""
    sam_split -i <input-sam-file>
    
    Parameters:
        -i/--input-sam-file    [string  :    path to SAM file from  BWA                      ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./              ]
        -c/--chromosomes       [string  :    chr1-22,X                                       ]
        -n/--number            [int     :    max of number of hit on chr1-22,X Default:  5   ]
        -k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep        ]
    
    Version:                    1.0.0
          """)
#parameters
input_sam_file = ''
output_dir = ''
#is_rm_tmp=True
chromosomes= ''
number=0

def setDefault():
    global output_dir
    output_dir = './'
    chromosomes='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X'+ \
                'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16'+ \
                'chr17,chr18,chr19,chr20,chr21,chr22,chrX'
    number=5

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:c:n:",["help",
                                                     "input-sam-file=",
                                                     "output-dir=",
                                                     "chromosomes=",
                                                     "number="])
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
        elif opt in ("-i", "--input-sam-file"):
            global input_sam_file
            input_sam_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-c", "--chromosomes"):
            global chromosomes
            chromosomes = arg
        elif opt in ("-n", "--number"):
            global number
            chromosomes = arg

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

chrDict = {}

def createDict():
    global chromosomes
    tmp=chromosomes.split(",")
    for x in range(len(tmp)):
        chrDict[tmp[x]]=1

def isMulti(chr1,others):
    global number
    tmp=others.split(";")
    chrlist = []
    chrlist.append(chr1)
    for x in range(len(tmp)-1):
        chrlist.append((tmp[x].split(","))[0])
    num=0
    for x in range(len(chrlist)):
        if chrlist[x] in chrDict:
            num=num+1
    if num <=number:
        return True
    else:
        return False


def run_split():
    infile = "%s" % input_sam_file
    f=open(infile,"r")
    outfile1 = "%s" % output_dir+'/tmp/'+getFilename(input_sam_file) + ".uniq.sam"
    f1=open(outfile1, "w")
    outfile2 = "%s" % output_dir+'/tmp/'+getFilename(input_sam_file) + ".unmapped.sam"
    f2=open(outfile2, "w")
    outfile3 = "%s" % output_dir+'/tmp/'+getFilename(input_sam_file) + ".repeats.sam"
    f3=open(outfile3, "w")
    outfile4 = "%s" % output_dir+'/tmp/'+getFilename(input_sam_file) + ".multihits.sam"
    f4=open(outfile4, "w")    
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            iskeep=False
            whichFile=0
            if len(line)>0:
                if line[0]=='#' or line[0]=='@':
                    iskeep=True
                else:
                    nf=len(line.split("\t"))
                    if nf==19:
                        UorR=(((line.split("\t"))[11]).split(":"))[2]
                        if UorR=='U':           #### mapped uniq and one location
                            whichFile=1
                            iskeep=True
                        elif UorR=='R':         #### if the -n/N for BWA setted very high, then here means repeats 
                            whichFile=3         #### if -n/N is small, no info to decide the number of hit on chr1-22,X
                            iskeep=True
                    elif nf==12:                #### unmapped
                        whichFile=2
                    elif nf==20:                #### (1) best hit (X0=1) is uniq, but has other hits (2) R: have many hits 
                        chr1=(line.split("\t"))[0]
                        others=(((line.split("\t"))[19]).split(":"))[2]
                        if isMulti(chr1,others):
                            whichFile=4
                            iskeep=True
                        else:
                            whichFile=3
                            iskeep=True 
            if iskeep==True:
                if whichFile==0:
                    f1.write(line)
                    f2.write(line)
                    f3.write(line)
                    f4.write(line)
                elif whichFile==1:
                    f1.write(line)
                elif whichFile==2:
                    f2.write(line)
                elif whichFile==3:
                    f3.write(line)
                elif whichFile==4:
                    f4.write(line)
    f.close()
    f1.close()
    f2.close()
    f3.close()
    f4.close()

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_sam_file=='':
        usage()
        exit(1);
    else:
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        createDict()
        run_split()   
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
