#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print """
    call_clusters -i <input-bed-file-1,input-bed-file-1,...> -n <name> 
    
    Parameters:
        -i/--input-bed-files   [string  :    path to bed files                               ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./              ]
        -n/--name              [string  :    name of the merged bed Default: all             ]
        -b/--bedtools          [string  :    path to bedtools       Default: bedtools        ]
        -k/--keep-tmp          [             keep the tmp folder.   Default: not keep        ]
    
    Version:                    1.0.0
          """
#parameters
input_bed_files = ''
output_dir = ''
is_rm_tmp=None
merge_name=''
path_to_bedtools=''

def setDefault():
    global output_dir
    output_dir = './'
    global path_to_bedtools
    path_to_bedtools = 'bedtools'
    global merge_name
    merge_name='all'
    global is_rm_tmp
    is_rm_tmp=True

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hki:o:b:n:",["help",
                                                      "keep-tmp",
                                                      "input-bed-files=",
                                                      "output-dir=",
                                                      "bedtools=",
                                                      "name="])
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
        elif opt in ("-i", "--input-bed-files"):
            global input_bed_files
            input_bed_files = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-b", "--bedtools"):
            global path_to_bedtools   
            path_to_bedtools=arg
        elif opt in ("-n", "--name"):
            global merge_name
            merge_name = arg

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


bed_files = []

def get_bed_files():
    tmp=input_bed_files.split(",")
    for x in range(len(tmp)):
        bed_files.append(tmp[x])

def merge_sort_all():
    cmd = 'cp ' + bed_files[0] + ' ' + output_dir + '/'+merge_name+'.bed'   
    print cmd
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output
    for x in range(len(bed_files)):
        if x != 0:
            cmd = 'cat ' + bed_files[x] + '>>'+output_dir+'/'+merge_name+'.bed'
            print cmd
            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            output = p.stdout.read()
            print output  
    cmd = path_to_bedtools + ' sort -i '+output_dir+'/'+merge_name+'.bed > '+output_dir+'/'+merge_name+'.sort.bed'
    print cmd
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

def get_strand_from_bed(s):
    ss=(s.split("\t")[3]).split(",")[3][0]
    return ss

def get_SDB_from_bed(s):
    ss=(s.split("\t")[3]).split(",")[1][0]
    return ss

def get_where_to_go(s):
    ss=get_strand_from_bed(s)
    sdb=get_SDB_from_bed(s)
    if sdb=="B":
        return 2
    elif sdb=="S":
        if ss=='+':
            return 0
        else:
            return 1
    elif sdb=="D":
        if ss=='+':
            return 1
        else:
            return 0

def split_with_strand():
    outfile = "%s" % output_dir+'/'+merge_name+'.sort.bed'
    f=open(outfile, "r")
    outfile2 = "%s" % output_dir+'/'+merge_name+'.sort.0.bed'
    f2=open(outfile2, "w")
    outfile3 = "%s" % output_dir+'/'+merge_name+'.sort.1.bed'
    f3=open(outfile3, "w")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            if get_where_to_go(line)==0:
                f2.write(line)
            elif get_where_to_go(line)==1:
                f3.write(line)
            elif get_where_to_go(line)==2:
                f2.write(line)
                f3.write(line)
    f.close()
    f2.close()
    f3.close()
    
def run_call_cluster(file_name):
    cmd = path_to_bedtools + ' merge -i '+output_dir + '/'+file_name+' -nms > '+ getOutputName2(file_name,'merge.bed')  
    print cmd
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print output

def call():
    run_call_cluster(merge_name+'.sort.0.bed')
    run_call_cluster(merge_name+'.sort.1.bed')
    
def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_files=='':
        usage()
        exit(1);
    else:
        use_real_path()
        make_dir(output_dir)
        make_dir(output_dir+'/tmp')
        get_bed_files()
        merge_sort_all()
        split_with_strand()
        call()
        remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
