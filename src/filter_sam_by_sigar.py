#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print ("""
    filter_sam_by_sigar -i <input-sam-file>
    
    Parameters:
        -i/--input-sam-file    [string  :    path to SAM file from  BWA                       ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./               ]
        -s/--sigars-allow      [string  :    file of allowed sigars Default: sigars_allow.txt ]
    
    Version:                    1.0.0
          """)
#parameters
input_sam_file = ''
output_dir = ''
#is_rm_tmp=True
sigars_allow_file= ''

def setDefault():
    global output_dir
    output_dir = './'
    cur=os.path.dirname(os.path.abspath(__file__))
    global sigars_allow_file
    sigars_allow_file=cur+"/"+"sigars_allow.txt"

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:",["help",
                                                    "input-sam-file=",
                                                    "output-dir=",
                                                    "sigars-allow="])
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
        elif opt in ("-s", "--sigars-allow"):
            global sigars_allow_file
            sigars_allow_file = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0o755 ) 

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

def isStringInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

sigar_range = []

def getSigarRange():
    infile = "%s" % sigars_allow_file
    f=open(infile,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            tmp=line.split("\t")
            sigar_str=tmp[0]
            aa = []
            aa.append(sigar_str)
            ranges=tmp[1].split(";")
            for x in range(len(ranges)):
                tt=ranges[x].split(",") 
                bb =[]
                bb.append(int(tt[0]))
                if not isStringInt(tt[1]):
                    bb.append(sys.maxsize)
                else:
                    bb.append(int(tt[1]))
                aa.append(bb)
            sigar_range.append(aa)

def isGoodSigar(sigar):
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
    is_in_file=False
    for x in range(len(sigar_range)):
        if s==sigar_range[x][0]:
            is_in_file=True
            isGood = True
            for y in range(len(sigar_range[x][0])):
                l=sigar_range[x][y+1][0]
                r=sigar_range[x][y+1][1]
                if nums[y]<l or nums[y]>r:
                    isGood=False
            if isGood ==False:
                return False
    if is_in_file==False:
        return False
    return True

def run_filter():
    infile = "%s" % input_sam_file
    f=open(infile,"r")
    outfile = "%s" % output_dir+"/tmp/"+getFilename(input_sam_file)+".no.bad.sigar.sam"
    f2=open(outfile, "w")    
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            iskeep=False
            if len(line)>0:
                if line[0]=='#' or line[0]=='@':
                    iskeep=True
                else:
                    iskeep=isGoodSigar((line.split("\t"))[5])
            if iskeep==True:
                f2.write(line)
    f.close()    
    f2.close()

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_sam_file=='' or sigars_allow_file=='':
        usage()
        exit(1);
    else:
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        getSigarRange()
        run_filter()
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
