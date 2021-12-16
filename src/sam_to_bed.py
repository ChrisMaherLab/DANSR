#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print ("""
    sam_to_bed -i <input-sam-file> -t <type>
    
    Parameters:
        -i/--input-sam-file    [string  :    path to SAM file from  BWA                                                     ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./                                             ]
        -t/--type              [string  :    "single" or "pair"                                                             ]
        -p/--pair-type         [string  :    "fr-unstranded"/"fr-firststranded"/"fr-secondstranded"  Default: fr-unstranded ]
        -s/--single-type       [string  :    "forward"/"reverse"/"both"                              Default: forward       ]
    
    Version:                    1.0.0
          """)
#parameters
input_sam_file = ''
output_dir = ''
#is_rm_tmp=None
isSingle=None
pair_type=""
singe_type=""

def setDefault():
    global output_dir
    output_dir = './'
    global pair_type
    pair_type="fr-unstranded"
    global singe_type
    singe_type="forward"
    #global is_rm_tmp
    #is_rm_tmp=True

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:t:p:s:",["help",
                                                       "input-sam-file=",
                                                       "output-dir=",
                                                       "type=",
                                                       "pair-type=",
                                                       "single-type="])
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
        elif opt in ("-t", "--type"):
            global isSingle
            if arg=="single":
                isSingle=True
            elif arg=="pair":
                isSingle=False
            else:
                print ("-t, choose single or pair")
                exit(1)
        elif opt in ("-p", "--pair-type"):
            global pair_type
            pair_type=arg
            if pair_type !="fr-unstranded" and pair_type !="fr-firststranded" and pair_type !="fr-secondstranded":
                print ("-p, choose from fr-unstranded/fr-firststranded/fr-secondstranded")
                exit(1)
        elif opt in ("-s", "--single-type"):
            global single_type
            single_type = arg
            if single_type !="forward" and single_type !="reverse" and single_type !="both":
                print ("-s, choose from forward/reverse/both")
                exit(1)

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

def denary2Binary(n):
    bStr = ''
    if n < 0:  raise ValueError("must be a positive integer")
    if n == 0: return '0'
    while n > 0:
        bStr = str(n % 2) + bStr
        n = n >> 1
    return bStr

def get_strand_from_flag(f):
    bStr=denary2Binary(int(f))
    if len(bStr)<5 or bStr[len(bStr)-5]=='0':
        return "+"
    else:
        return "-"

def get_isFirst_from_flag(f):
    bStr=denary2Binary(int(f))
    if len(bStr)<7:
        return None
    else:
        if bStr[len(bStr)-7]=='1':
            return True
        else:
            return False

def sam_to_compact(s):
    tmp=s.split("\t")
    chr1=tmp[2]
    strand=get_strand_from_flag(tmp[1])
    position=tmp[3]
    cigar=tmp[5]
    edit=tmp[12].split(":")[2]
    return chr1+","+strand+position+","+cigar+","+edit

def getNextPosition(p,cigar):
    if 'S' in cigar or 'H' in cigar: ##### not supporting clipped reads
        return -1
    else:
        s=''
        nums=[]
        num_str=''
        for x in range(len(cigar)):
            if cigar[x] in {'M','I','D','S','H'}:
                s=s+cigar[x]
                nums.append(int(num_str))
                num_str=''
            else:
                num_str=num_str+cigar[x]
        totalLen=0
        for x in range(len(s)):
            if s[x]=='M' or s[x]=='D':
                totalLen=totalLen+nums[x]
        return p+totalLen-1

def getSDB(isFirst):
    if isSingle==True:
        if single_type == "forward":
            return "S"
        elif single_type == "reverse":
            return "D"
        elif single_type == "both":
            return "B"
    else:
        if pair_type=="fr-unstranded":
            if isFirst==True:
                return "B"
        elif pair_type=="fr-firststranded":
            if isFirst==True:
                return "D"
            else:
                return "S"
        elif pair_type=="fr-secondstranded":
            if isFirst==True:
                return "S"
            else:
                return "D"
    print ("Warning: Something is wrong at from getSDB")
    return "S"

def one_compact_to_bed(name,c,f):
    tmp=c.split(",")
    chr1=tmp[0]
    position=tmp[1][1:]
    cigar=tmp[2]
    nextPos=getNextPosition(int(position),cigar)
    if nextPos==-1:
        return ""
    else:
        return chr1+"\t"+str(int(position)-1)+"\t"+str(nextPos)+"\t"+name+","+getSDB(get_isFirst_from_flag(f))+","+c

def compact_to_bed(name,f,compacts):
    final=''
    tmp=compacts.split(";")
    for x in range(len(tmp)-1):
        cc=one_compact_to_bed(name,tmp[x],f) 
        if cc!='':
            final=final+cc+"\n"
    return final

def run_convert():
    infile = "%s" % input_sam_file
    f=open(infile,"r")
    outfile2 = "%s" % output_dir+'/tmp/'+getOutputName2(input_sam_file,"bed")
    f2=open(outfile2, "w")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            if len(line)>0:
                if line[0]!='#' and line[0]!='@':
                    UorR=(((line.split("\t"))[11]).split(":"))[2]
                    if UorR=='U' and len(line.split("\t"))==19:      #### mapped uniq and one location
                        f2.write(one_compact_to_bed(line.split("\t")[0],sam_to_compact(line),line.split("\t")[1])+"\n")
                    if len(line.split("\t"))==20: #### (1) best hit (X0=1) is uniq, but has other hits (2) R: have many hits
                        others=(((line.split("\t"))[19]).split(":"))[2]
                        all_compact=sam_to_compact(line)+";"+others
                        f2.write(compact_to_bed(line.split("\t")[0],line.split("\t")[1],all_compact))
    f.close()
    f2.close()

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
        run_convert()
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
