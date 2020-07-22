#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print ("""
    bed_split_and_merge.py -i <comma separated bed_files>
    
    Parameters:
        -i/--input-bed-file-list    [string  :    path to SAM file from  BWA                   ]
        -o/--output-dir             [string  :    path to output dir.    Default: ./           ]
        -s/--sample-name            [string  :    sample name            Default: sample       ]
        -b/--bedtools               [string  :    path to bedtoos.       Default: bedtools     ]
    
    Version:                    1.0.0
          """)
#parameters
input_bed_files_list = ''
output_dir = ''
#is_rm_tmp=None
path_to_bedtools = ''
sample_name = ''

def setDefault():
    global output_dir
    output_dir = './'
    #global is_rm_tmp
    #is_rm_tmp=True
    global path_to_bedtools
    path_to_bedtools = 'bedtools'
    global sample_name
    sample_name = "sample"

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:b:o:s:",["help",
                                                      "input-bed-file-list=",
                                                      "bedtools=",
                                                      "output-dir=",
                                                      "sample-name=",])
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
        elif opt in ("-i", "--input-bed-files-list"):
            global input_bed_files_list
            input_bed_files_list = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-b", "--bedtools"):
            global path_to_bedtools
            path_to_bedtools = arg
        elif opt in ("-s", "--sample-name"):
            global sample_name
            sample_name = arg

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
    return fn

def run_convert():
    infiles_names = "%s" % input_bed_files_list
    infiles = infiles_names.split(",")
    outfile1 = "%s" % (output_dir+"/tmp/"+sample_name+".0.bed")
    outfile2 = "%s" % (output_dir+"/tmp/"+sample_name+".1.bed")
    f1=open(outfile1, "w")
    f2=open(outfile2, "w")
    for index in range(len(infiles)):
        f = open(infiles[index],"r")
        while True:
            line=f.readline()
            if line=="":
                break
            else:
                tmp=line.split("\t")
                tp = tmp[3].split(",")[1]
                strand = tmp[3].split(",")[3][0]
                index_out_0 = 0
                index_out_1 = 0
                if strand=="+" and (tp=="S" or tp=="B"):
                    index_out_0=1
                if strand=="+" and (tp=="D" or tp=="B"):
                    index_out_1=1
                if strand=="-" and (tp=="S" or tp=="B"):
                    index_out_1=1
                if strand=="-" and (tp=="D" or tp=="B"):
                    index_out_0=1
                if index_out_0 == 1:
                    f1.write(line)
                if index_out_1 == 1:
                    f2.write(line)                        
        f.close()
    f1.close()
    f2.close()

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'

def run_bedtools():
    bedfile1 = "%s" % (output_dir+"/tmp/"+sample_name+".0.bed")
    bedfile2 = "%s" % (output_dir+"/tmp/"+sample_name+".1.bed")
    bedoutfile1 = "%s" % (output_dir+"/tmp/"+sample_name+".sort.merge.0.bed")
    bedoutfile2 = "%s" % (output_dir+"/tmp/"+sample_name+".sort.merge.1.bed")    

    #cmd = path_to_bedtools + " sort" + " -i "+ bedfile1 + " | " + path_to_bedtools + " merge -nms > " + bedoutfile1 
    cmd = path_to_bedtools + " sort" + " -i "+ bedfile1 + " | " + path_to_bedtools + " merge -i - -c 4 -o collapse -delim \";\" > " + bedoutfile1   ### modified by Eteleeb
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)
    p.communicate()
    if p.returncode:
        print('{} returned with exit code {}, aborting.'.format(cmd, p.returncode), file=sys.stderr)
        sys.exit(1)


    #cmd = path_to_bedtools + " sort" + " -i "+ bedfile2 + " | " + path_to_bedtools + " merge -nms > " + bedoutfile2
    cmd = path_to_bedtools + " sort" + " -i "+ bedfile2 + " | " + path_to_bedtools + " merge -i - -c 4 -o collapse -delim \";\" > " + bedoutfile2   ### modified by Eteleeb
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)
    p.communicate()
    if p.returncode:
        print('{} returned with exit code {}, aborting.'.format(cmd, p.returncode), file=sys.stderr)
        sys.exit(1)


def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_files_list=='':
        usage()
        exit(1);
    else:
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        run_convert()
        run_bedtools() 
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
