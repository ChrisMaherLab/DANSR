#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename
from difflib import get_close_matches

def usage():
    print ("""
    dansr -i <reads.fastq>

    Parameters:
        
        -h/--help              [        :    Print help information                 Optional                                  ]
        -o/--output-dir        [string  :    Path to output dir                     Optional        Default:  ./              ]
        -k/--keep-tmp          [        :    Keep the tmp folder                    Optional        Default: Do not keep      ]
        -S/--sample-name       [string  :    Sample name                            Optional        Default: sample           ]   
        -z/--setup-file        [string  :    Path to setup file                     Optional        Default: setup.small.ini  ] 
        -v/--verbose           [        :    Print details                          Optional                                  ]
        -w/--skip              [        :    Skip a step if target output exists    Optional        Default: Not skip         ] 

        -b/--begin-no-trimming [        :    Input already trimmed                  Optional        Default: Trim             ]
        -a/--adapter           [string  :    Adapter sequence to cut                Required                                  ]
        -i/--input-file        [string  :    Path to single reads (FASTQ)           Required                                  ]
                               [             Or first reads in pair                                                           ]
        -j/--input-file-2      [string  :    Path to second reads in pair (FASTQ)   Optional                                  ]
        -x/--cutadapter-opts   [string  :    cutadaptor options in ""               Optional        Default: ""               ]
        -A/--adapter2          [string:      3\' adapter for reads 2                Optional        Default: ""               ]   

        -r/--reference         [string  :    Path human reference genome            Required                                  ]
        -X/--bwa-options-aln   [string  :    bwa-aln options in ""                  Optional        Default: "-q 5 -l 17 -k 1"]
        -Y/--bwa-options-sam   [string  :    bwa-samse/sampe option in ""           Optional        Default: ""               ]

        -c/--chromosomes       [string  :    Comma separated chr names              Optional        Default: chr1-22,X        ]
        -n/--number-hits       [int     :    max number of hits on chr1-22,X        Optional        Default:  5               ]

        -p/--pair-type         [string  :    "fr-unstranded"/"fr-firststranded"/"fr-secondstranded" Default: fr-unstranded    ]
        -s/--single-type       [string  :    "forward"/"reverse"/"both"                             Default: forward          ]

        -N/--number-reads      [int     :    Min number of reads for a cluster                      Default: 2                ]    

        -P/--percent-cur       [float   :    Max percentage of reads from precurser Optional        Default: 0.3              ]  
        -f/--cutoff            [float   :    RPM increase cutoff for boudary update Optional        Default: 0.33             ] 

        -l/--list-of-gtf       [string  :    Comma separated GTF files for all      Required                                  ]
    
        -U/--percent-uniq      [float   :    Min percent uniq reads in a cluster    Optional        Default: 0.5              ] 
        -R/--uniq-reads        [int:    :    Min number of uniq reads in a cluster  Optional        Default: 2                ]  
        -V/--ov-with-largest   [float   :    Min % of overlap with largest cluster  Optional        Default 0.75              ]

        -g/--gtf-small         [string  :    Comma separated GTFs for small RNA     Required                                  ]
        -e/--small-types       [string  :    Comma separated small RNA types        Required                                  ]
        -u/--unstranded        [        :    No strand awareness                    Optional        Default: stranded         ]
        -J/--jaccard-index     [float   :    Min Jaccard similarity score.          Optional        Default: 0.3              ]

    Version:                   1.0.1
          """)

#parameters
output_dir = ''
is_rm_tmp = True
sample_name =  ''
setup_file = ''
silence = True
skip_step = False
begin_no_cut = False
adapter = ''
input_file = ''
input_file_2 = ''
cutadapter_opts = ''
adapter2 = ''
reference = ''
bwa_options_aln = ''
bwa_options_sam = ''
chromosomes = ''
number_hits = ''          
type_str = ''
pair_type = ''
single_type = ''
number_reads = ''

percent_cur = ''
cutoff = ''

percent_uniq = ''
uniq_reads = ''
ov_with_largest = ''

gtf_small = ''
#percent_over = ''
small_types = ''
#min_novo_reads = ''
no_strand = False
jaccard_index = ''

path_to_cutadapt = ''
path_to_bwa = ''
path_to_bedtools = ''

python_path="python"
cur_dir="."

def setDefault():
    global output_dir
    output_dir = './'
    global cur_dir
    cur_dir=os.path.dirname(os.path.abspath(__file__))
    global is_rm_tmp
    is_rm_tmp = True
    global sample_name
    sample_name = 'sample'
    global setup_file
    setup_file=cur_dir+"/"+"setup.small.ini"
    global begin_no_cut
    begin_no_cut = False
    global silence
    silence = True
    global skip_step
    skip_step = False
    global python_path
    python_path=sys.executable
    global pair_type
    pair_type = "fr-unstranded"
    global single_type
    single_type = "forward"
    global type_str
    type_str = "single" # can update at run BWA
    global number_reads
    number_reads = '2'

def get_tool_path():
    tmp_dict = {}
    f=open(setup_file,"r")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            tmp=line.rstrip('\n').split("\t")
            tmp_dict[tmp[0]]=tmp[1]
    global path_to_cutadapt
    global path_to_bwa
    global path_to_bedtools
    path_to_cutadapt = tmp_dict['cutadapt']
    path_to_bwa = tmp_dict['bwa']
    path_to_bedtools = tmp_dict['bedtools']
    f.close()

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)
    global reference
    reference=os.path.realpath(reference)   

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"ho:vwkS:z:ba:i:j:x:A:r:X:Y:c:n:p:s:N:P:f:l:g:e:uR:U:V:J:",
       ["help",
        "output-dir",
        "verbose",
        "skip",
        "keep-tmp",
        "sample-name=",
        "setup-file=",
        "begin-no-trimming",
        "adapter=",
        "input-file=",
        "input-file-2=",
        "cutadapter-opts=",
        "adapter2=",
        "reference=",
        "bwa-options-aln=",
        "bwa-options-sam=",
        "chromosomes=",
        "number-hits=",
        "pair-type=",
        "single-type=",
        "number-reads=",
        "percent-cur=",
        "cutoff=",
        "list-of-gtf=",
        "gtf-small=",
        "small-types=",
        "unstranded",
        "uniq-reads=",
        "percent-uniq=",
        "ov-with-largest=",
        "jaccard-index="])

    except getopt.GetoptError:
        print ("getopt.GetoptError")
        usage()
        sys.exit(1)
    for opt, arg in opts:
        print (opt,arg)
        if opt in ("-h","--help"):
            usage()
            sys.exit(0)
        elif opt in ("-o","--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-v","--verbose"):
            global silence
            silence = False
        elif opt in ("-w","--skip"):
            global skip_step
            skip_step = True
        elif opt in ("-k","--keep-tmp"):
            global is_rm_tmp
            is_rm_tmp  = False
        elif opt in ("-S","--sample-name"):
            global sample_name
            sample_name = arg
        elif opt in ("-z","--setup-file"):
            global setup_file
            setup_file = arg
        elif opt in ("-b","--begin-no-trimming"):
            global begin_no_cut
            begin_no_cut = True
        elif opt in ("-a","--adapter"):
            global adapter
            adapter = arg
        elif opt in ("-i","--input-file"):
            global input_file
            input_file = arg
        elif opt in ("-j","--input-file-2"):
            global input_file_2
            input_file_2 = arg
        elif opt in ("-x","--cutadapter-opts"):
            global cutadapter_opts
            cutadapter_opts = arg
        elif opt in ("-A","--adapter2"):
            global adapter2
            adapter2 = arg
        elif opt in ("-r","--reference"):
            global reference
            reference = arg
        elif opt in ("-X","--bwa-options-aln"):
            global bwa_options_aln
            bwa_options_aln = arg
        elif opt in ("-Y","--bwa-options-sam"):
            global bwa_options_sam
            bwa_options_sam = arg
        elif opt in ("-c","--chromosomes"):
            global chromosomes
            chromosomes = arg
        elif opt in ("-n","--number-hits"):
            global number_hits
            number_hits = arg
        elif opt in ("-p","--pair-type"):
            global pair_type
            pair_type = arg
        elif opt in ("-s","--single-type"):
            global single_type
            single_type = arg
        elif opt in ("-N","--number-reads"):
            global number_reads
            number_reads = arg
        elif opt in ("-P","--percent-cur"):
            global percent_cur
            percent_cur = arg
        elif opt in ("-f","--cutoff"):
            global cutoff
            cutoff = arg
        elif opt in ("-l","--list-of-gtf"):
            global list_of_gtf
            list_of_gtf = arg
        elif opt in ("-g","--gtf-small"):
            global gtf_small
            gtf_small = arg
        #elif opt in ("-q","--percent-over"):
        #    global percent_over
        #    percent_over = arg
        elif opt in ("-e","--small-types"):
            global small_types
            small_types = arg
        #elif opt in ("-m","--min-novo-reads"):
        #    global min_novo_reads
        #    min_novo_reads = arg
        elif opt in ("-u","--unstranded"):
            global no_strand
            no_strand = 'True'
        elif opt in ("-U","--percent-uniq"):
            global percent_uniq
            percent_uniq = arg
        elif opt in ("-R","--uniq-reads"):
            global uniq_reads
            uniq_reads = arg
        elif opt in ("-V","--ov-with-largest"):
            global ov_with_largest
            ov_with_largest = arg
        elif opt in ("-J","--jaccard-index"):
            global jaccard_index
            jaccard_index = arg
    if not input_file_2:
        type_str="single"
        if single_type != "forward":
            print('Error: this version of DANSR supports only forward stranded single read libraries', file=sys.stderr)
            sys.exit(1)
    else:
        type_str="pair"
        print('Warning: paired end mode is still under beta development, take care to check results for accuracy.', file=sys.stderr)

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

def getFormat(filename):
    tmp=filename.split(".")
    fm=tmp[len(tmp)-1]
    return fm

def getOutputName(filename,insert):
    fn=getFilename(filename)
    fm=getFormat(filename)
    return fn+"."+insert+"."+fm

def getOutputName2(filename,append_format):
    fn=getFilename(filename)
    return fn+"."+append_format

def run_cutadapt():
    global input_file_2
    global input_file
    if silence == False:
        print
        print ("########### in run_cutadapt ###############")
        print ("""
runCutadapt.py parameters:

-a/--adapter           [string:    3' ADAPTER                                        ]
-i/--input-file        [string:    input file                                        ]
-o/--output-dir        [string:    path to output dir                                ]
-m/--min-read-length   [int:       minimum read length left.   Default:17            ]
-x/--cutadapter-opts   [string:    cutadaptor options.         Default: ''           ]
-A/--adapter2          [string:    3' adapter for reads 2      Default: ''           ]
-j/--input-file-2      [string:    input file 2                Default: ''           ]
-c/--cutadapt          [string:    path to cutadapt            Default: cutadapt     ]
-k/--keep-tmp          [      :    keep the tmp folder.        Default: not using -k ]

        """)
        print ("Set as the following:")
        print ("adapter:",adapter)
        print ("input-file:",input_file)
        print ("output-dir:",output_dir)
        print ("min-read-length:","17")
        print ("cutadapter-opts:",cutadapter_opts)
        print ("adapter2:",adapter2)
        print ("input-file-2",input_file_2)
        print ("cutadapt:",path_to_cutadapt)
        #ot="Yes"
        #if is_rm_tmp == True:
        #    ot="No"
        #print ("keep-tmp:",ot)
        #print

    if begin_no_cut:
        print ("Note    ###############     ' ","User have already done trimming, return"," ' ############### exsits, skip ...")
        return

    is_return = False

    # refer runCutadapt.py: line 132
    target_file = output_dir+"/tmp/"+getOutputName(input_file,'trim')
    if skip_step and os.path.exists(target_file):
        print ("Warning ############### >>>> ",target_file," <<<< ############### exsits, skip ...")
        is_return = True
    if input_file_2 !='':
        target_file = output_dir+"/tmp/"+getOutputName(input_file_2,'trim')
        if skip_step and os.path.exists(target_file):
            print ("Warning ############### >>>> ",target_file_2," <<<< ############### exsits, skip ...")
            is_return = True
    if is_return:
        return
    
    cmd = python_path + ' ' + cur_dir+"/runCutadapt.py" + ' -a ' + adapter + ' -i ' + input_file + ' -o ' + output_dir + ' -m 17 ' + '-c ' + path_to_cutadapt
    #if is_rm_tmp == False:
    #    cmd = cmd + ' -k'
    if adapter2!='':
        cmd = cmd + ' -A ' + adapter2 + ' -j ' + input_file_2
    if cutadapter_opts!='':
        cmd = cmd + ' -x ' + cutadapter_opts
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if silence == False:
        print (output)

def run_bwa():
    global input_file
    global input_file_2
    if silence == False:
        print
        print ("########### in run_bwa ###############")
        print ("""
runBWA.py parameters:

-i/--input-file        [string  :    path to trimmed input               ]
-o/--output-dir        [string  :    path to output dir.    Default: ./  ]
-s/--sample-name       [string  :    sample name, use if set             ]
-r/--reference         [string  :    path to the reference               ]
-x/--bwa-options-aln   ["string":    BWA options aln.       Default:     ]
-y/--bwa-options-sam   ["string":    BWA options sam.       Default:     ]
-k/--keep-tmp          [        :    keep the tmp folder.                ]
-j/--input-file-2      [string  :    path to trimmed input2 Default: ''  ]
-b/--bwa               [string  :    path to BWA                         ]

        """)
        print ("Set as the following:")
        print ("input-file:", "IS FROM HERE", input_file)
        print ("input-file-2:", "IS FROM HERE", input_file_2)
        print ("output-dir:", output_dir)
        print ("sample-name:", sample_name)
        print ("reference:", reference)
        print ("bwa-options-aln:", bwa_options_aln)
        print ("bwa-options-sam:", bwa_options_sam)
        print ("bwa:", path_to_bwa)
        #ot="Yes"
        #if is_rm_tmp == True:
        #    ot="No"
        #print ("keep-tmp:",ot)
        print
    #refere to runCutadpp.py line 132
    if not begin_no_cut:
        input_file = output_dir+"/tmp/"+getOutputName(input_file,'trim')
        if input_file_2!='':
            input_file_2 = output_dir+"/tmp/"+getOutputName(input_file_2,'trim')
            #for sam to bed
            global type_str
            type_str = "pair"

    print ("+++++++++++++ updated input_file to, ",input_file)

    #### refer to runBWA.py line 144, but the target is using sample_name in runBWA
    target_file = output_dir+"/alignment/"+getOutputName2(sample_name,'sam')
    print ("$$$$$$$$$$$$$$$",target_file)
    if skip_step and os.path.exists(target_file):
        print ("Warning ############### >>>> ",target_file," <<<< ############### exsits, skip ...")
        return 

    cmd = python_path + ' ' + cur_dir+"/runBWA.py" + ' -r ' + reference + ' -i ' + input_file + ' -o ' + output_dir + ' -b ' + path_to_bwa
    #if is_rm_tmp == False:
    #    cmd = cmd + ' -k'
    if input_file_2!='':
        cmd = cmd + ' -j ' + input_file_2
    if bwa_options_aln!='':
        cmd = cmd + ' -x ' + '\"' + bwa_options_aln + '\"'
    if bwa_options_sam!='':
        cmd = cmd + ' -y ' + '\"' + bwa_options_sam + '\"'
    if sample_name!="":
        cmd = cmd + " -s "  + sample_name
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if silence == False:
        print (output)

def run_filter_sam_by_sigar():
    global input_file
    if silence == False:
        print
        print ("########### in run_filter_sam_by_sigar ###############")
        print ("""

run_filter_sam_by_sigar parameters:

-i/--input-sam-file    [string  :    path to SAM file from  BWA                       ]
-o/--output-dir        [string  :    path to output dir.    Default: ./               ]
-s/--sigars-allow      [string  :    file of allowed sigars Default: sigars_allow.txt ]
-k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep

        """)

        print ("Set as the following:")
        print ("input-sam-file:", "IS FROM HERE")
        print ("output-dir:", output_dir)
        print ("sigars-allow:", "Use sigars-allow.txt")
        #ot="Yes"
        #if is_rm_tmp == True:
        #    ot="No"
        #print ("keep-tmp:",ot)
        print


    #### refer to runBWA.py lines 137 and 145 # but in run BWA it is changed to sample_name
    input_file = getOutputName2(sample_name,'sam')

    print ("+++++++++++++ updated input_file to, ",input_file)

    #### refer to run filter_sam_by_sigar.py line 144
    target_file = output_dir+"/tmp/"+getFilename(input_file)+".no.bad.sigar.sam"
    print ("$$$$$$$$$$$$",target_file)
    if skip_step and os.path.exists(target_file):
        print ("Warning ############### >>>> ",target_file," <<<< ############### exsits, skip ...")
        return 

    cmd = python_path + ' ' + cur_dir+"/filter_sam_by_sigar.py" +  ' -i ' + output_dir+'/alignment/'+input_file + ' -o ' + output_dir
    #if is_rm_tmp == False:
    #    cmd = cmd + ' -k'
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if silence == False:
        print (output)

def run_sam_split():
    global input_file
    if silence == False:
        print ("########### in run_sam_split ###############")
        print ("""

run_sam_split parameters:

-i/--input-sam-file    [string  :    path to SAM file from  BWA                      ]
-o/--output-dir        [string  :    path to output dir.    Default: ./              ]
-c/--chromosomes       [string  :    chr1-22,X                                       ]
-n/--number            [int     :    max of number of hit on chr1-22,X Default:  5   ]
-k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep        ]

        """)

        print ("Set as the following:")
        print ("input-sam-file:", "IS FROM HERE")
        print ("output-dir:", output_dir)
        print ("chromosomes:", chromosomes)
        print ("number:", number_hits)
        #ot="Yes"
        #if is_rm_tmp == True:
        #    ot="No"
        #print ("keep-tmp:",ot)
        print

        print ("Command:")

    #### refer to run filter_sam_by_sigar.py line 144
    input_file = getFilename(input_file)+".no.bad.sigar.sam"

    #### refer to sam_split.py line 122
    target_file = output_dir+"/tmp/"+getFilename(input_file) + ".uniq.sam"
    print ("$$$$$$$$$$$$",target_file)
    if skip_step and os.path.exists(target_file):
        print ("Warning ############### >>>> ",target_file," <<<< ############### exsits, skip ...")
        return
	
    cmd = python_path + ' ' + cur_dir+"/sam_split.py" +  ' -i ' + output_dir+'/tmp/'+input_file + ' -o ' + output_dir

    if chromosomes!='':
        cmd = cmd + ' -c ' + chromosomes
    if number_hits!="": 
        cmd = cmd + ' -n ' + number_hits
    #if is_rm_tmp == False:
    #    cmd = cmd + ' -k'
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if silence == False:
        print (output)


def run_sam_to_bed():
    global input_file
    global input_file_2
    if silence == False:
        print ("########### in run_sam_to_bed ###############")
        print ("""

run_sam_to_bed parameters:

-i/--input-sam-file    [string  :    path to SAM file from  BWA                                                     ]
-o/--output-dir        [string  :    path to output dir.    Default: ./                                             ]
-t/--type              [string  :    "single" or "pair"                                                             ]
-p/--pair-type         [string  :    "fr-unstranded"/"fr-firststranded"/"fr-secondstranded"  Default: fr-unstranded ]
-s/--single-type       [string  :    "forward"/"reverse"/"both"                              Default: forward       ]
-k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep

        """)

        print ("Set as the following:")
        print ("input-sam-file:", "IS FROM HERE")
        print ("output-dir:", output_dir)
        print ("type:", type_str)
        print ("pair-type:", pair_type)
        print ("single-type:", single_type)
        #ot="Yes"
        #if is_rm_tmp == True:
        #    ot="No"
        #print ("keep-tmp:",ot)
        print

    #### refer to sam_split.py line 122
    tmp_name = input_file
    input_file = getFilename(input_file) + ".uniq.sam"
    input_file_2 = getFilename(tmp_name) + ".multihits.sam"  

    #### refer to sam_to_bed.py line 211
    target_file_1 = output_dir+"/tmp/"+getOutputName2(input_file,"bed")
    print ("$$$$$$$$$$$$",target_file_1)
    if skip_step and os.path.exists(target_file_1):
        print ("Warning ############### >>>> ",target_file_1," <<<< ############### exsits, skip ...")
    else:
        cmd = python_path + ' ' + cur_dir+"/sam_to_bed.py" +  ' -i ' + output_dir+'/tmp/'+input_file + ' -o ' + output_dir + ' -t ' + type_str

        if pair_type!='':
            cmd = cmd + ' -p ' + pair_type
        if single_type!='':
            cmd = cmd + ' -s ' + single_type
        #if is_rm_tmp == False:
        #    cmd = cmd + ' -k'
        print (cmd)
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        if silence == False:
            print (output)

    target_file_2 = output_dir+"/tmp/"+getOutputName2(input_file_2, "bed")
    print ("$$$$$$$$$$$$",target_file_2)
    if skip_step and os.path.exists(target_file_2):
        print ("Warning ############### >>>> ",target_file_2," <<<< ############### exsits, skip ...")
    else:
        cmd = python_path + ' ' + cur_dir+"/sam_to_bed.py" +  ' -i ' + output_dir+'/tmp/'+input_file_2 + ' -o ' + output_dir + ' -t ' + type_str

        if pair_type!='':
            cmd = cmd + ' -p ' + pair_type
        if single_type!='':
            cmd = cmd + ' -s ' + single_type
        #if is_rm_tmp == False:
        #    cmd = cmd + ' -k'
        print (cmd)
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        if silence == False:
            print (output)


def run_bed_split_and_merge():
    global input_file
    global input_file_2
    if silence == False:
        print ("########### in run_bed_split_and_merge ###############")
        print ("""

run_sam_to_bed parameters:

-i/--input-bed-file-list    [string  :    path to SAM file from  BWA                   ]
-o/--output-dir             [string  :    path to output dir.    Default: ./           ]
-s/--sample-name            [string  :    sample name            Default: sample       ]
-b/--bedtools               [string  :    path to bedtoos.       Default: bedtools     ]
-k/--keep-tmp               [        :    keep the tmp folder.   Default: not keep     ]

        """)

        print ("Set as the following:")
        print ("input-bed-file-list:", "IS FROM HERE")
        print ("output-dir:", output_dir)
        print ("sample-name:", sample_name)
        print ("bedtools", path_to_bedtools)
        #ot="Yes"
        #if is_rm_tmp == True:
        #    ot="No"
        #print ("keep-tmp:",ot)
        print

    #### refer to sam_to_bed.py line 211
    input_file = getOutputName2(input_file,"bed")
    input_file_2 = getOutputName2(input_file_2,"bed")

    bedoutfile1 = "%s" % (output_dir+"/"+sample_name+".sort.merge.0.bed")

    #### refer to bed_split_and_merge.py line 132
    target_file_1 = output_dir+"/tmp/"+getOutputName2(sample_name,"sort.merge.0.bed")
    target_file_2 = output_dir+"/tmp/"+getOutputName2(sample_name,"sort.merge.1.bed")
    print ("$$$$$$$$$$$$",target_file_1,target_file_2)
    if skip_step and os.path.exists(target_file_1) and os.path.exists(target_file_2):
        print ("Warning ############### >>>> ",target_file_1," and ",target_file_1," <<<< ############### exsit, skip ...")
        return    
    cmd = python_path + ' ' + cur_dir+"/bed_split_and_merge.py" +  ' -i ' + output_dir+'/tmp/'+input_file+','+output_dir+'/tmp/'+input_file_2 + ' -o ' + output_dir + ' -s ' + sample_name + ' -b ' + path_to_bedtools
    #if is_rm_tmp == False:
    #    cmd = cmd + ' -k'
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if silence == False:
        print (output)


def run_filter_bed():
    global input_file
    global input_file_2
    if silence == False:
        print ("########### in run_filter_bed ###############")
        print ("""

run_filter_bed parameters:

-i/--input-bed-file    [string  :    path to bed file with column 4 as ; seperated alignments                       ]
-o/--output-dir        [string  :    path to output dir.    Default: ./                                             ]
-k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep                                       ]
-n--number-of-reads    [int     :    min number of reads to keep a bed record   Default: 2                          ]

        """)

        print ("Set as the following:")
        print ("input-bed-file:", "IS FROM HERE")
        print ("output-dir:", output_dir)
        print ("number-of-reads", number_reads)
        #ot="Yes"
        #if is_rm_tmp == True:
        #    ot="No"
        #print ("keep-tmp:",ot)
        print

    #### refer to bed_split_and_merge.py line 132
    input_file = getOutputName2(sample_name,"sort.merge.0.bed")
    input_file_2 = getOutputName2(sample_name,"sort.merge.1.bed")

    #### refer to filter_bed.py line 83
    target_file_1 = output_dir+"/tmp/"+getOutputName2(input_file,"f.bed")
    target_file_2 = output_dir+"/tmp/"+getOutputName2(input_file_2,"f.bed")
   
    print ("$$$$$$$$$$$$",target_file_1)
    if skip_step and os.path.exists(target_file_1):
        print ("Warning ############### >>>> ",target_file_1," <<<< ############### exsits, skip ...")
    else:
        cmd = python_path + ' ' + cur_dir+"/filter_bed.py" +  ' -i ' + output_dir+'/tmp/'+input_file + ' -o ' + output_dir + ' -n ' + number_reads
        #if is_rm_tmp == False:
        #    cmd = cmd + ' -k'
        print (cmd)
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        if silence == False:
            print (output)

        print ("$$$$$$$$$$$$",target_file_2)
    if skip_step and os.path.exists(target_file_2):
        print ("Warning ############### >>>> ",target_file_2," <<<< ############### exsits, skip ...")
    else:
        cmd = python_path + ' ' + cur_dir+"/filter_bed.py" +  ' -i ' + output_dir+'/tmp/'+input_file_2 + ' -o ' + output_dir + ' -n ' + number_reads 
        #if is_rm_tmp == False:
        #    cmd = cmd + ' -k'
        print (cmd)
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        if silence == False:
            print (output)

def run_update_boundary_2():
    global input_file
    global input_file_2
    if silence == False:
        print ("########### in run_filter_bed ###############")
        print ("""

        -i/--input-bed-file    [string  :    path to BED with reads                                          ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./                              ]
        -p/--percent           [float   :    percentage of reads at each side to consider                    ]
        -c/--cutoff            [float   :    increase cutoff to change boundary                              ]
        -k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep                        ]

        """)
        print ("Set as the following:")
        print ("input-bed-file:", "IS FROM HERE")
        print ("output-dir:", output_dir)
        print ("percent:", percent_cur)
        print ("cutoff:", cutoff)
        #ot="Yes"
        #if is_rm_tmp == True:
        #    ot="No"
        #print ("keep-tmp:",ot)
        print

    #### refer to filter_bed.py line 83
    input_file = getOutputName2(input_file,"f.bed")
    input_file_2 = getOutputName2(input_file_2,"f.bed")

    #### refer to update_boundary_2.py line 347
    target_file_1 = output_dir+"/tmp/"+getOutputName2(input_file,"update.bed")
    target_file_2 = output_dir+"/tmp/"+getOutputName2(input_file_2,"update.bed")

    print ("$$$$$$$$$$$$",target_file_1)
    if skip_step and os.path.exists(target_file_1):
        print ("Warning ############### >>>> ",target_file_1," <<<< ############### exsits, skip ...")
    else:
        cmd = python_path + ' ' + cur_dir+"/update_boundary_2.py" +  ' -i ' + output_dir+'/tmp/'+input_file + ' -o ' + output_dir
        #if is_rm_tmp == False:
        #    cmd = cmd + ' -k'
        if percent_cur !='':
            cmd = cmd + ' -p ' + percent_cur
        if cutoff !='':
            cmd = cmd + ' -c ' + cutoff
        print (cmd)
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        if silence == False:
            print (output)

        print ("$$$$$$$$$$$$",target_file_2)
    if skip_step and os.path.exists(target_file_2):
        print ("Warning ############### >>>> ",target_file_2," <<<< ############### exsits, skip ...")
    else:
        cmd = python_path + ' ' + cur_dir+"/update_boundary_2.py" +  ' -i ' + output_dir+'/tmp/'+input_file_2 + ' -o ' + output_dir
        #if is_rm_tmp == False:
        #    cmd = cmd + ' -k'
        if percent_cur !='':
            cmd = cmd + ' -p ' + percent_cur
        if cutoff !='':
            cmd = cmd + ' -c ' + cutoff        
        print (cmd)
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        if silence == False:
            print (output)

def run_bed_to_graph():
    global input_file
    global input_file_2
    if silence == False:
        print ("########### in run_bed_to_graph ###############")
        print ("""
        inputs order: XX.0.f.updagte.bed XX.1.f.updagte.bed list_of_gtf result.txt result.bed

        """)
        print ("Other parameters:","IS FROM HERE")
        print ("list-of-gtf:", list_of_gtf)
        print ("Command:")

    #### refer to update_boundary_2.py line 347
    input_file = output_dir+'/tmp/'+getOutputName2(input_file,"update.bed")
    input_file_2 = output_dir+'/tmp/'+getOutputName2(input_file_2,"update.bed")

    #refer to bedToGraph inputs
    target_file_1 = output_dir+"/tmp/"+"result."+sample_name+".txt"
    target_file_2 = output_dir+"/tmp/"+"result."+sample_name+".bed"

    if skip_step and os.path.exists(target_file_1) and os.path.exists(target_file_2):
        print ("Warning ############### >>>> ",target_file_1,target_file_2," <<<< ############### exsit, skip ...")
        return

    cmd = cur_dir+"/"+"bedToGraph "+input_file+" "+input_file_2+" "+list_of_gtf+" "+target_file_1+" "+target_file_2
    print (cmd)       
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if silence == False:
        print (output)

def run_filter_result():

    input_file = output_dir+"/tmp/"+"result."+sample_name+".bed"

    #if skip_step and os.path.exists(target_file_1) and os.path.exists(target_file_2):
    #    print "Warning ############### >>>> ",target_file_1,target_file_2," <<<< ############### exsit, skip ..."
    #    return

    #cmd = python_path + ' ' + cur_dir+'/filter_result_bed.py -i '+input_file+' -0 '+input_file_bed1+' -1 '+input_file_bed2 + ' -o ' + output_dir
    cmd = python_path + ' ' + cur_dir+'/filter_result.py -i '+input_file+' -o '+output_dir

    #if is_rm_tmp == False:
    #    cmd = cmd + ' -k'
    if percent_uniq!='':
        cmd = cmd + ' -p ' + percent_uniq
    if uniq_reads!='':
        cmd = cmd + ' -u ' + uniq_reads
    if number_reads!='':
        cmd = cmd + ' -c ' + number_reads
    if ov_with_largest!='':
        cmd = cmd + ' -v ' + ov_with_largest

    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if silence == False:
        print (output)

def run_assign_best_annotation():
        
    input_file = output_dir+"/tmp/"+"result."+sample_name+".annotated.tsv"
    
    cmd = python_path + ' ' + cur_dir +'/assign_best_annotation.py -i ' + input_file +' -o ' + output_dir + ' -b ' + path_to_bedtools

    #if is_rm_tmp == False:
    #    cmd = cmd + ' -k'
    if list_of_gtf!='':
        cmd = cmd + ' -g ' + list_of_gtf
    if jaccard_index!='':
        cmd = cmd + ' -j ' + jaccard_index
    if no_strand:
        cmd = cmd + ' -u '
    if small_types:
        cmd += '-s '+ small_types
        
    print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    if silence == False:
        print (output)

#def run_filter_result_bed():
#    global input_file
#    global input_file_2
#    if silence == False:
#        print "########### in run_filter_bed ###############"
#        print """
#
#-i/--input-bed-file    [string  :    path to BED with reads                                          ]
#-0/--cluster-strand0   [string  :    path to clusters file on strand 0                               ]
#-1/--cluster-strand1   [string  :    path to clusters file on strand 1                               ]
#-o/--output-dir        [string  :    path to output dir.    Default: ./                              ]
#-p/--percent           [float   :    min percentage of unique reads                                  ]
#-m/--min-unique        [int:    :    min number of unique reads                                      ]
#-c/--cutoff            [int:    :    min number of reads for a cluster                               ]
#-k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep                        ]
#
#        """
#
#        print "Set as the following:"
#        print "input-bed-file:", "IS FROM HERE"
#        print "cluster-strand0", "IS FROM HERE"
#        print "cluster-strand1", "IS FROM HERE"
#        print "output-dir:", output_dir
#        print "percent:", percent_uniq, "from -U/--percent-uniq here"  # -p
#        print "uniq-reads:",uniq_reads, "from -R/--uniq-reads here" # -m 
#        print "cutoff:", number_reads, "from -N/--number-reads, shared with run_filter_bed" # -c
#        ot="Yes"
#        if is_rm_tmp == True:
#            ot="No"
#        print "keep-tmp:",ot
#        print
#
#    ####
#    input_file_bed1 = input_file
#    input_file_bed2 = input_file_2
#    input_file = output_dir+"/"+"result."+sample_name+".bed"
#
#    #### target
#    target_file_1 = output_dir+"/"+"result."+sample_name+".filtered.bed"
#    target_file_2 = output_dir+"/"+"result."+sample_name+".passed.bed"
#
#    if skip_step and os.path.exists(target_file_1) and os.path.exists(target_file_2):
#        print "Warning ############### >>>> ",target_file_1,target_file_2," <<<< ############### exsit, skip ..."
#        return
#
#    cmd = python_path + ' ' + cur_dir+'/filter_result_bed.py -i '+input_file+' -0 '+input_file_bed1+' -1 '+input_file_bed2 + ' -o ' + output_dir
#    if is_rm_tmp == False:
#        cmd = cmd + ' -k'
#    if percent_uniq!='':
#        cmd = cmd + ' -p ' + percent_uniq
#    if uniq_reads!='':
#        cmd = cmd + ' -m ' + uniq_reads
#    if number_reads!='':
#        cmd = cmd + ' -c ' + number_reads
#    print cmd
#    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
#    output = p.stdout.read()
#    if silence == False:
#        print output

#def run_assign_best_annotation():
#    if silence == False:
#        print
#        print "########### in run_assign_best_annotation ###############"
#        print """
#
#assign_best_annotation.py:
#
#-i/--input-bed-file    [string  :    path to BED with reads                                                       ]
#-o/--output-dir        [string  :    path to output dir.    Default: ./                                           ]
#-g/--gtf               [string  :    path to gtf file of small RNA. (Or , seperated paths)                        ]
#-p/--percent           [float   :    min overlapping ratio                                                        ]
#-s/--small-types       [string  :    comma separted small RNA types                                               ]
#-m/--min-novo-reads    [int     :    minimum number of reads to call unannotated novel small RNA. Default: 5.     ]
#-n/--no-strand         [        :    no strand awareness                                                          ]
#-k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep                                     ]
#        """
#
#        print "Set as the following:"
#        print "input-bed-file:", "IS FROM HERE"
#        print "output-dir:", output_dir
#        print "gtf:",gtf_small, "from -g/--gtf-small here"
#        print "percent:", percent_over, "from -q/--percent-over here"
#        print "small-types:", small_types
#        print "min-novo-reads:",min_novo_reads, "from -m/--min-novo-reads here"
#        ns="No"
#        if no_strand == True:
#            ns = "Yes"
#        print "no-strand:",ns, "from -u/--no-strand here"
#        ot="Yes"
#        if is_rm_tmp == True:
#            ot="No"
#        print "keep-tmp:",ot
#        print
#
#    ####
#    input_file = output_dir+"/"+"result."+sample_name+".passed.bed"
#
#    #### run with out checking existence
#
#    cmd = python_path + ' ' + cur_dir+'/assign_best_annotation.py -i '+input_file+' -g '+gtf_small+' -s '+small_types +' -o ' + output_dir
#    if is_rm_tmp == False:
#        cmd = cmd + ' -k'
#    if percent_over!='':
#        cmd = cmd + ' -p ' + percent_over
#    if min_novo_reads!='':
#        cmd = cmd + ' -m ' + min_novo_reads
#    print cmd
#    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
#    output = p.stdout.read()
#    if silence == False:
#        print output
#
def remove_tmp():
    if is_rm_tmp:
        cmd = 'rm -rf ' + output_dir +'/tmp'

def warnings_and_usage():
    if input_file == '':
        print ("-i/--input-file is required")
        exit(1)
        
    # cutadapt
    if begin_no_cut == False:
        if adapter == '':
            print ("-a/--adapter is required")
            exit(1)
        if not((input_file_2 =='' and adapter2 =='') or (input_file_2!='' and adapter2!='' )):
            print ("-A/--adapter2 and -j/--input-file-2 have to be set or not set at the same time")
            exit(1)
    # bwa
    if reference == '':
        print ("-r/--reference is required")
        exit(1)
    # if begin_no_cut == True:
    #     print ("here",input_file)
    #     if input_file == '':
    #         print ("-i/--input-file is required")
    #         exit(1)
    # other
    if type_str == '':
        print ("-t/--type is required")
        exit(1)
    if list_of_gtf == '':
        print ("-l/--list-of-gtf is required")
        exit(1)
    if gtf_small == '':
        print ("-g/--gtf-small is required")
        exit(1)
    # if small_types == '':
    #     print ("-e/--small-types is required")
    #     exit(1)

    # classification 

    # assign best annotation 

def main(argv):
    setDefault()
    get_tool_path()
    getParameters(argv[1:])
    warnings_and_usage()
    use_real_path()
    make_dir(output_dir)
    make_dir(output_dir+'/tmp')
    make_dir(output_dir+'/alignment')
    make_dir(output_dir+'/results')
    run_cutadapt()
    run_bwa()
    run_filter_sam_by_sigar()
    run_sam_split()
    run_sam_to_bed()
    run_bed_split_and_merge()
    run_filter_bed()
    run_update_boundary_2()
    run_bed_to_graph()
    run_filter_result()
    run_assign_best_annotation()
    remove_tmp()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
