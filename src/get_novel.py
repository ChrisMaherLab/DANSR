#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print """
    get_novel -i <bedfile> -g <gtf_list>
    
    Parameters:
        -i/--input-bed-file    [string  :    path to BED with reads                                          ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./                              ]
        -g/--gtf               [string  :    path to gtf file of small RNA. (Or , seperated paths)           ]
        -d/--degree            [int     :    max degree allowed                              Default:   5    ]
        -a/--alignments        [int     :    min num of alignments required                  Default:   5    ]
        -u/--unique            [int     :    min num of uniq alignments required             Default:   2    ]
        -l/--overlap           [int     :    min num of reads overlapping with other reads   Default:   3    ]
        -p/--percent-unique    [float   :    min percent of uniq alignments                  Default:   0.75 ]
        -q/--precent-overlap   [float   :    min percent of overlapping reads                Default:   0.75 ]
        -k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep                        ]
    
    Version:                    1.0.0
          """

#parameters
input_bed_file = ''
output_dir = ''
is_rm_tmp=None
gtf_file = ''
n_degree = 100
n_alignments = 0
n_unique = 0
u_overlap = 0
p_unique = 0.0
p_overlap = 0.0

def setDefault():
    global output_dir
    output_dir = './'
    global is_rm_tmp
    is_rm_tmp=True
    global n_degree
    n_degree = 5
    global n_alignments
    n_alignments = 5
    global n_unique
    n_unique = 2
    global n_overlap
    n_overlap = 3
    global p_unique
    p_unique =0.75
    global p_overlap
    p_overlap = 0.75

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hki:o:g:d:a:u:l:p:q:",["help",
                                                                "keep-tmp",
                                                                "input-bed-file=",
                                                                "output-dir=",
                                                                "gtf=", 
                                                                "degree=",
                                                                "alignments=",
                                                                "unique=",
                                                                "overlap=",
                                                                "percent-unique=",
                                                                "percent-overlap="])
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
        elif opt in ("-g", "--gtf"):
            global gtf_file
            gtf_file=arg
        elif opt in ("-d", "--degree"):
            global n_degree
            n_degree=int(arg)
        elif opt in ("-a", "--alignments"):
            global n_alignments
            n_alignments=int(arg)
        elif opt in ("-u", "--unique"):
            global n_unique
            n_unique=int(arg)
        elif opt in ("-l", "--overlap"):
            global n_overlap
            n_overlap=int(arg)
        elif opt in ("-p", "--percent-unique"):
            global p_unique
            p_unique=float(arg)
        elif opt in ("-q", "--percent-overlap"):
            global p_overlap
            p_overlap=float(arg)
        elif opt in ("-p", "--percent"):
            global percent
            percent=arg

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

def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

def rm_quoat(aa):
    if aa[0]=='\"' and aa[len(aa)-1]==";":
        return aa[1:len(aa)-2]
    else:
        return aa

all_annot_segs = []

def get_all_small_RNA():
    global all_annot_segs
    all_annot_segs=[]
    gtfs=gtf_file.split(",")
    for x in range(len(gtfs)):
        f=open(gtfs[x],"r")
        while True:
            line=f.readline()
            if line=="":
                break
            else:
                tmp=line.split("\t")
                if tmp[2]!="exon":
                    continue      
                tmp2=tmp[8].split(" ")
                isGeneType=False
                isGeneName=False
                tt=''
                nn=''
                for x in range(len(tmp2)):
                    if tmp2[x]=='gene_name':
                        isGeneName=True
                        continue
                    if tmp2[x]=='gene_biotype':
                        isGeneType=True
                        continue
                    if isGeneName:
                        nn=rm_quoat(tmp2[x])
                        isGeneName=False
                    if isGeneType:
                        tt=rm_quoat(tmp2[x])
                        isGeneType=False
                if tt!='' and nn!='':
                    record=[]
                    record.append(tmp[0])
                    record.append(tmp[3])
                    record.append(tmp[4])
                    record.append(nn)
                    record.append("0")
                    record.append(tmp[6])
                    record.append(tt)
                    all_annot_segs.append(record)
        print "Here the size of all_annot_segs is",len(all_annot_segs)
        f.close()

tmp_bed_annot_file = ""

def print_records():
    name=""
    gtfs=gtf_file.split(",")
    for x in range(len(gtfs)):
        path, filename = os.path.split(gtfs[x])
        filename = getFilename(filename)
        if name=="":
            name=filename
        else:
            name=name+"_"+filename

    global tmp_bed_annot_file
    tmp_bed_annot_file=output_dir +'/tmp/'+name+'.bed'

    print "Name of the bed from gtf: ",tmp_bed_annot_file

    if os.path.exists(tmp_bed_annot_file)==False:
        get_all_small_RNA()
        f=open(tmp_bed_annot_file,"w")
        print "I am writing the bed file",tmp_bed_annot_file,"with",str(len(all_annot_segs)),"records"
        for x in range(len(all_annot_segs)):
            tmp=all_annot_segs[x]
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6]))
        f.close()
        cmd = 'bedtools sort -i '+tmp_bed_annot_file+' > '+ output_dir+'/tmp/'+getOutputName2(tmp_bed_annot_file,"sort.bed")
        print cmd
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        print output
        p.communicate()
        if p.returncode:
            print('{} returned with exit code {}, aborting.'.format(cmd, p.returncode), file=sys.stderr)
            sys.exit(1)


def intersect():
    cmd = 'bedtools intersect -a '+input_bed_file+' -b '+output_dir+'/tmp/'+getOutputName2(tmp_bed_annot_file,"sort.bed")
    cmd = cmd + ' -wa -wb'
    cmd = cmd + ' -f '+str(percent)+ ' -r'
    if no_strand != False:
        cmd = cmd + " -s"
    cmd = cmd + ' > '+output_dir+"/tmp/"+getOutputName2(input_bed_file,"append.bed")
    print cmd
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output
    p.communicate()
    if p.returncode:
        print('{} returned with exit code {}, aborting.'.format(cmd, p.returncode), file=sys.stderr)
        sys.exit(1)


cluster_to_bed_record = {}


# chage this to a different logic: if in overlapping, not print
# if not pass threshold, not print
def get_best_annot():
    f=open(output_dir+"/tmp/"+getOutputName2(input_bed_file,"append.bed"),"r")
    while True:
        line = f.readline()
        if line == "":
            break
        else:
            tmp = line.split("\t")
            s1 = int(tmp[1])
            e1 = int(tmp[2])
            s2 = int(tmp[16])
            e2 = int(tmp[17])
            l1 = e2-s1
            l2 = e1-s2
            lm = min(l1,l2)
            l1 = e1-s1
            l2 = e2-s2
            ol = 1.0*lm/(l1+l2-lm)
            if tmp[3] not in cluster_to_bed_record:
                cluster_to_bed_record[tmp[3]]=(line).rstrip('\n')+"\t"+str(ol)
            else:
                if float(cluster_to_bed_record[tmp[3]].rstrip('\n').split("\t")[22]) < ol:
                    cluster_to_bed_record[tmp[3]]=(line).rstrip('\n')+"\t"+str(ol)
    f.close()
    outfile = "%s" % output_dir+'/'+getOutputName2(input_bed_file,"annotated.bed")
    f = open(outfile,"w")
    for c in cluster_to_bed_record:
        tmp=cluster_to_bed_record[c].split("\t")
        for x in range(len(tmp)-1):
            f.write("%s\t" % (tmp[x]))
        f.write("%s\n" % tmp[len(tmp)-1])
    f.close()

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_file=='' or gtf_file=='':
        usage()
        exit(0);
    else:
        use_real_path()
        make_dir(output_dir)
        make_dir(output_dir+'/tmp')
        print_records()
        intersect() 
        get_best_annot()
        remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))






