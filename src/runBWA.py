#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print ("""
    runBWA -r <reference> -i <input-file> -o <output-dir>
    runBWA -r <reference> -i <input-file> -j <input-file-2> -o <output-dir>

    Parameters:
        -i/--input-file        [string  :    path to trimmed input               ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./  ]
        -s/--sample-name       [string  :    sample name, use if set             ]
        -r/--reference         [string  :    path to the reference               ]
        -x/--bwa-options-aln   ["string":    BWA options aln.       Default:     ]
                               [                                                 ]
        -y/--bwa-options-sam   ["string":    BWA options sam.       Default:     ]
                               [                                                 ]
        -j/--input-file-2      [string  :    path to trimmed input2 Default: ''  ]
        -b/--bwa               [string  :    path to BWA                         ] 
    
    Version:                    1.0.0
          """)

#parameters
input_file = ''
input_file_2 = ''
output_dir = ''
path_to_bwa = ''
#is_rm_tmp=True
bwa_options_aln='-t 4 -q 5 -l 17 -k 1'
bwa_options_sam='-n 100 -N 100'
reference=''
sample_name = ''

def setDefault():
    global output_dir
    output_dir = './'
    global path_to_bwa
    path_to_bwa = 'bwa'
    global bwa_options_aln
    bwa_options_aln = '-q 5 -l 17 -k 1'
    global bwa_options_sam
    bwa_options_sam = ''

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)
    global reference
    reference=os.path.realpath(reference)   

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:j:o:s:b:x:y:r:",["help",
                                                              "input-file=",
                                                              "input-file-2=",
                                                              "output-dir=",
                                                              "sample-name=",
                                                              "path_to_bwa=",
                                                              "bwa_options_aln=",
                                                              "bwa_options_sam=",
							      "reference="])
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
        elif opt in ("-i", "--input-file"):
            global input_file
            input_file = arg
        elif opt in ("-j", "--input-file-2"):
            global input_file_2
            input_file_2 = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-s", "--sample-name"):
            global sample_name
            sample_name = arg
        elif opt in ("-b", "--path_to_bwa"):
            global path_to_bwa
            path_to_bwa = arg
        elif opt in ("-x", "--bwa-options-aln"):
            global bwa_options_aln
            bwa_options_aln = arg
        elif opt in ("-y", "--bwa-options-sam"):
            global bwa_options_sam
            bwa_options_sam = arg
        elif opt in ("-r", "--reference"):
            global reference
            reference = arg

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

def run_bwa_aln():
    cmd = path_to_bwa + ' aln '
    if bwa_options_aln!='':
        cmd = cmd + bwa_options_aln + ' '
    if sample_name =='':
        cmd = cmd + '-f ' + output_dir + '/alignment/'+getOutputName2(input_file,'sai')+' '
    else:
        cmd = cmd + '-f ' + output_dir + '/alignment/'+sample_name+'.sai '
    cmd = cmd + reference +' '+ input_file;
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)

def run_bwa_aln_2():
    run_bwa_aln()
    cmd = path_to_bwa + ' aln '
    if bwa_options_aln!='':
        cmd = cmd + bwa_options_aln + ' '
    if sample_name =='':
        cmd = cmd + '-f ' + output_dir + '/alignment/'+getOutputName2(input_file_2,'sai')+' '
    else:
        cmd = cmd + '-f ' + output_dir + '/alignment/'+sample_name+'.2.sai '
    cmd = cmd + reference +' '+ input_file_2;
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)

def run_bwa_sam():
    if sample_name =='':
        cmd = path_to_bwa+ ' samse '+ bwa_options_sam +' -f ' + output_dir +'/alignment/'+getOutputName2(input_file,'sam') +' '
        cmd = cmd + reference +' '+ output_dir +'/alignment/'+getOutputName2(input_file,'sai')+' '+input_file
    else:
        cmd = path_to_bwa+ ' samse '+ bwa_options_sam +' -f ' + output_dir +'/alignment/'+sample_name+'.sam ' 
        cmd = cmd + reference +' '+ output_dir +'/alignment/'+sample_name+'.sai '+input_file
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)

def run_bwa_sam_2():
    if sample_name =='':
       cmd = path_to_bwa+ ' sampe '+ bwa_options_sam +' -f ' + output_dir +'/alignment/'+getOutputName2(input_file,'sam') +' '
       cmd = cmd + reference +' '+ output_dir +'/alignment/'+getOutputName2(input_file,'sai') + ' '
       cmd = cmd + output_dir +'/alignment/'+getOutputName2(input_file_2,'sai') + ' '
    else:
       cmd = path_to_bwa+ ' sampe '+ bwa_options_sam +' -f ' + output_dir +'/alignment/'+sample_name+'.sam '
       cmd = cmd + reference +' '+ output_dir +'/alignment/'+sample_name+'.sai '
       cmd = cmd + output_dir +'/alignment/'+sample_name+'.2.sai '
    cmd = cmd + input_file + ' ' + input_file_2
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)

def run_samtools():
    if sample_name == '':
        cmd = 'samtools view -Sb -o '+output_dir+'/alignment/'+getOutputName2(input_file,'bam')+' '
        cmd = cmd + output_dir +'/alignment/'+getOutputName2(input_file,'sam')
    else:
        cmd = 'samtools view -Sb -o '+output_dir+'/alignment/'+sample_name+'.bam '
        cmd = cmd + output_dir +'/alignment/'+sample_name+'.sam'
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)
    if sample_name == '':
        cmd = 'samtools sort -o '+output_dir +'/alignment/'+getFilename(input_file)+'.sort.bam '+output_dir+'/alignment/'+getOutputName2(input_file,'bam')+' '
    else:
        cmd = 'samtools sort -o 'output_dir +'/alignment/'+sample_name+'.sort.bam '+output_dir+'/alignment/'+sample_name+'.bam '
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)
    if sample_name == '':
        cmd = 'samtools index '+output_dir+'/alignment/'+getFilename(input_file)+'.sort.bam'
    else:
        cmd = 'samtools index '+output_dir+'/alignment/'+sample_name+'.sort.bam'
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_file_2=='':
        if input_file=='' or reference=='':
            usage()
            exit(1);
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        run_bwa_aln()
        run_bwa_sam()
        run_samtools()
        #remove_tmp()
    else:
        if input_file=='' or input_file_2=='' or reference=='':    
            usage()
            exit(1);
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        run_bwa_aln_2()
        run_bwa_sam_2()
        run_samtools()
        #remove_tmp()    

if __name__ == '__main__':
    sys.exit(main(sys.argv))
