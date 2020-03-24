#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print """
    filter_result_bed -i <bedfile> -0 <cluster-strand0-file> -1 <cluster-strand1-file>
    
    Parameters:
        -i/--input-bed-file    [string  :    path to BED with reads                                          ]
        -0/--cluster-strand0   [string  :    path to clusters file on strand 0                               ]
        -1/--cluster-strand1   [string  :    path to clusters file on strand 1                               ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./                              ]
        -p/--percent           [float   :    min percentage of unique reads                                  ]
        -m/--min-unique        [int:    :    min number of unique reads                                      ]
        -c/--cutoff            [int:    :    min number of reads for a cluster                               ]
        -k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep                        ]
    
    Version:                    1.0.0
          """
#parameters
input_bed_file = ''
output_dir = ''
is_rm_tmp=None
percent = 0.75
min_uniq = 1
cutoff = 2
file0 = ""
file1 = ""

def setDefault():
    global output_dir
    output_dir = './'
    global percent
    percent = 0.3
    global min_uniq
    min_uniq = 1
    global cutoff
    cutoff = 2
    global is_rm_tmp
    is_rm_tmp=True

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hki:o:p:c:m:0:1:",["help",
                                                            "keep-tmp",
                                                            "input-bed-file=",
                                                            "output-dir=",
                                                            "percent=",
                                                            "cutoff=",
                                                            "min-unique="
                                                            "cluster-strand0=",               
                                                            "cluster-strand1="])
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
            cutoff = int(arg)
        elif opt in ("-m", "--min-unique"):
            global min_uniq
            min_uniq = int(arg)
        elif opt in ("-0", "--cluster-strand0"):
            global file0
            file0 = arg
        elif opt in ("-1", "--cluster-strand1"):
            global file1
            file1 = arg

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

records = []
m_index = {} #"M:a;b" is 4,5 etc
header = ""
m_names = {} #a,b etc

    #if for clusters with M:, then share/M is the adjusted share
    #if the best past threshold, then keep all, otherwise filter all
    #for others, just filter

def run_load():   
    f=open(input_bed_file, "r")
    index=0
    global header
    header = f.readline()
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            tmp=line.rstrip("\n").split("\t")
            tmp.append(0)
            records.append(tmp)
            if tmp[13][0]=="M" and tmp[13][1]==":":
                m_name = tmp[13]
                strand = "0"
                if tmp[5]=="-":
                    strand = "1"
                cname = tmp[0]+"_"+str(int(tmp[1])+1)+"_"+str(int(tmp[2])-int(tmp[1]))+"_"+strand
                m_names[cname]=1
                if m_name not in m_index:
                    m_index[m_name]=str(index)
                else:
                    m_index[m_name]=m_index[m_name]+";"+str(index)
            index=index+1
    f.close()

m_reads = {}

def run_get_share_counts():
    global m_reads
    files = [file0,file1]
    for x in range(len(files)):
        filename = files[x]
        f=open(filename,"r")
        while True:
            line = f.readline()
            if line=="":
                break
            else:
                tmp = line.rstrip("\n").split("\t")
                cname = tmp[0]+"_"+str(int(tmp[1])+1)+"_"+str(int(tmp[2])-int(tmp[1]))+"_"+str(x)
                #print cname
                if cname in m_names:
                    read_names = []
                    aa=tmp[3].split(";")
                    for t in range(len(aa)):
                        read_names.append(aa[t].split(",")[0])
                    m_reads[cname]=read_names
        f.close()

def get_share_counts(m_line):
    global m_reads
    m_line=m_line.split(";")[0][2:]
    #m_line=m_line[2:]
    names=m_line.split(",")
    shared = {}
    zero=m_reads[names[0]]
    for x in range(len(zero)):
        shared[zero[x]]=1
    for y in range(len(names)):
        if y!=0:
            one = m_reads[names[y]]
            for t in range(len(one)):
                nm=one[t]
                if nm in shared:
                    shared[nm] = shared[nm]+1     
    num=0
    target=len(names)
    for x in shared:
        if shared[x]==target:
            num=num+1
    return num

def run_filter():
    not_m_pass = 0
    not_m_all = 0
    m_pass = 0
    m_all = 0
    for x in range(len(records)):
        tmp=records[x]
        tmp2=tmp[14].split("|")
        cluster_records=tmp2[0].split(";")
        alignments = int(cluster_records[2].split(":")[1])
        #print "alignments=",alignments
        unique = int(cluster_records[3].split(":")[1])
        #print "unique=",unique
        is_m = False
        if tmp[13][0]=="M" and tmp[13][1]==":":
            is_m = True
        if is_m == True:
            unique = get_share_counts(tmp[13])
        pass1 = False
        if is_m == True:
            m_all=m_all+1
            #print "unique=",unique,"alignments=",alignments,"unique*1.0/alignments=",unique*1.0/alignments
        else:
            not_m_all=not_m_all+1
            #print "unique=",unique,"alignments=",alignments,"unique*1.0/alignments=",unique*1.0/alignments
        if unique >= min_uniq and alignments >= cutoff and unique*1.0/alignments >= percent:
            pass1 = True
        if pass1 == True and is_m == False:
            records[x][15]=1
            not_m_pass=not_m_pass+1
        if pass1 == True and is_m == True:
            m_pass=m_pass+1
            indeces = m_index[tmp[13]].split(";")
            for i in range(len(indeces)):
                ii = int(indeces[i])
                records[ii][15]=1
    outfile1 = "%s" % output_dir+'/'+getOutputName2(input_bed_file,"passed.bed")
    outfile2 = "%s" % output_dir+'/'+getOutputName2(input_bed_file,"filtered.bed")
    #print "not_m pass, all = ",not_m_pass,not_m_all, "m pass, all = ",m_pass,m_all
    f1 = open(outfile1,"w")
    f2 = open(outfile2,"w")
    for i in range(len(records)):
        tmp = records[i]
        f=f1
        if tmp[15]==0:
            f=f2 
        for x in range(len(tmp)-2):
            f.write("%s\t" % (tmp[x]))
        f.write("%s\n" % (tmp[len(tmp)-2]))
    f1.close()
    f2.close()

def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_file=='' or file0=='' or file1=='':
        usage()
        exit(0);
    else:
        use_real_path()
        make_dir(output_dir)
        make_dir(output_dir+'/tmp')
        run_load()
        run_get_share_counts()
        run_filter()
        remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
