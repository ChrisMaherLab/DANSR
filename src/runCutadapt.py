#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print ("""
    runCutadapt -a ADAPTER -i input-file -o output-dir -m min-length -x "<cutadapt-options>"
    runCutadapt -a ADAPT1 -A ADAPT2 -i input-file-1 -j input-file-2 -o output-dir -x "<cutdapt-options>" 

    Parameters:
        -a/--adapter           [string:    3' ADAPTER                                        ]
        -i/--input-file        [string:    input file                                        ]
        -o/--output-dir        [string:    path to output dir                                ]
        -m/--min-read-length   [int:       minimum read length left.   Default:17            ]
        -x/--cutadapter-opts   ["string:"  cutadaptor options.         Default: ''           ]
        -A/--adapter2          [string:    3' adapter for reads 2      Default: ''           ]
        -j/--input-file-2      [string:    input file 2                Default: ''           ]
        -c/--cutadapt          [string:    path to cutadapt            Default: cutadapt     ]

    Version:                    1.0.0
          """)

#parameters
adapter = ''
input_file= ''
output_dir = ''
min_read_length=0
path_to_cutadapt = ''
adapter2 = ''
input_file_2 = ''
#is_rm_tmp=True
cutadapt_opts = ''

def setDefault():
    global output_dir
    output_dir = './'
    global path_to_cutadapt
    path_to_cutadapt = 'cutadapt'
    global min_read_length
    min_read_length=17

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir);


def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"ha:i:o:m:x:A:j:c:",["help",
                                                              "adapter=",
                                                              "input-file=",
                                                              "output-dir=",
                                                              "min-read-length=",
                                                              "cutadaper-opts=",
                                                              "adapter2=",
                                                              "input_file_2=",
                                                              "path-to-cutadpt="])
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
        elif opt in ("-a", "--adapter"):
            global adapter
            adapter = arg
        elif opt in ("-i", "--input-file"):
            global input_file
            input_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-m", "--min-read-length"):
            global min_read_length
            min_read_length = arg
        elif opt in ("-x", "--cutadapt-opts"):
            global cutadapt_opts
            cutadapt_opts = arg
        elif opt in ("-j", "--input-file-2"):
            global input_file_2
            input_file_2 = arg
        elif opt in ("-A", "--path-to-bwa"):
            global adapter2
            adapter2 = arg
        elif opt in ("-c", "--path-to-cutadapt"):
            global path_to_cutadapt
            path_to_cutadapt = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0o755 )

def getFormat(filename):
    tmp=filename.split(".")     
    fm=tmp[len(tmp)-1]
    return fm

def getFilename(filename):
    bn=basename(filename) 
    tmp=bn.split(".")
    fn=tmp[0]
    for x in range(1, len(tmp)-1):
        fn=fn+"."+tmp[x]
    return fn

def getOutputName(filename,insert):
    fn=getFilename(filename)
    fm=getFormat(filename)
    return fn+"."+insert+"."+fm
    

def run_cutadapt_1():
    cmd = path_to_cutadapt +' -a '+ adapter +' -o '+ output_dir +'/tmp/'+getOutputName(input_file,'trim')
    if cutadapt_opts!='':
        cmd=cmd +' '+ cutadapt_opts + ' '
    cmd = cmd + ' -m '+str(min_read_length)+' '+ input_file
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output) 

def run_cutadapt_2():
    cmd = path_to_cutadapt +' -a '+ adapter +' -A '+adapter2+' -o '+ output_dir +'/tmp/'+getOutputName(input_file,'trim')
    cmd = cmd + ' -p '+ output_dir +'/tmp/'+getOutputName(input_file_2,'trim')
    if cutadapt_opts!='':
        cmd=cmd +' '+ cutadapt_opts + ' '
    cmd = cmd + ' -m '+ str(min_read_length) +' '+ input_file + ' ' + input_file_2
    p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    print (output)

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'
#        p = Popen(cmd, cwd=output_dir, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
#        output = p.stdout.read()
#        print (output)

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_file_2=='':
        if input_file!='' and adapter=='':
            usage()
            exit(1)
        if adapter=='' or input_file=='' or output_dir=='':
            usage()
            exit(1)
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        run_cutadapt_1()
        #remove_tmp()
    else:    
        if (input_file!='' and adapter=='') or (input_file_2!='' and adapter2==''):
            usage()
            exit(1)
        if adapter=='' or input_file=='' or output_dir=='' or adapter2=='' or input_file=='':
            usage()
            exit(1)
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        run_cutadapt_2()
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
