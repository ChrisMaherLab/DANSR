#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns; sns.set()

def usage():
    print ("""
    assign_best_annotation -i <bedfile> -g <gtf> -s <str:small1,small2>
    
    Parameters:
        -i/--input-bed-file    [string  :    path to BED with reads                                                       ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./                                           ]
        -g/--gtf               [string  :    path to gtf file of small RNA. (Or , seperated paths)                        ]
        -j/--j-score           [float   :    min Jaccard similarity. Default: 0.3                                         ]
        -s/--small-types       [string  :    comma separated small RNA types                                               ]
        -m/--min-novo-reads    [int     :    minimum number of reads to call unannotated novel small RNA. Default: 5.     ]
        -u/--unstranded        [        :    no strand awareness. Default: stranded                                       ]
        -b/--bedtools          [string  :    path to bedtoos.       Default: bedtools                                     ]    
    Version:                    1.0.0
          """)
#parameters
input_bed_file = ''
output_dir = ''
#is_rm_tmp=None
gtf_file = ''
j_score = 0.3
unstranded = False
small_types_str = ''
path_to_bedtools = ''

def setDefault():
    global output_dir
    output_dir = './'
    #global is_rm_tmp
    #is_rm_tmp=True
    global j_score
    j_score = 0.3
    global unstranded 
    unstranded = False
    global path_to_bedtools
    path_to_bedtools = 'bedtools'

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:g:j:b:s:u",["help",
                                                       "input-bed-file=",
                                                       "output-dir=",
                                                       "gtf=", 
                                                       "j-score=",
                                                       "bedtools=",
                                                       "small-types=",
                                                       "unstranded"])
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
        elif opt in ("-u","--unstranded"):
            global unstranded
            unstranded = True
        elif opt in ("-i", "--input-bed-file"):
            global input_bed_file
            input_bed_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-g", "--gtf"):
            global gtf_file
            gtf_file=arg
        elif opt in ("-j", "--j-score"):
            global j_score
            j_score=float(arg)
        elif opt in ("-b", "--bedtools"):
            global path_to_bedtools
            path_to_bedtools = arg
        elif opt in ("-s", "--small-types"):
           global small_types_str
           small_types_str=arg
        #elif opt in ("-m", "--min-novo-reads"):
        #    global min_novo
        #    min_novo=int(arg)

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

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'

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
            if line.startswith("#!"):
                continue  
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
        #print ("Here the size of all_annot_segs is",len(all_annot_segs))
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

    #print ("Name of the bed from gtf: ",tmp_bed_annot_file)

    if os.path.exists(tmp_bed_annot_file)==False:
        get_all_small_RNA()
        f=open(tmp_bed_annot_file,"w")

        #print ("Writing the bed file",tmp_bed_annot_file,"with",str(len(all_annot_segs)),"records")
        for x in range(len(all_annot_segs)):
            tmp=all_annot_segs[x]
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('chr'+tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6]))
        f.close()
        cmd = 'bedtools sort -i '+tmp_bed_annot_file+' > '+ output_dir+'/tmp/'+getOutputName2(tmp_bed_annot_file,"sort.bed")
        #print cmd
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        p.communicate()
        if p.returncode:
            print('{} returned with exit code {}, aborting.'.format(cmd, p.returncode), file=sys.stderr)
            sys.exit(1)

        #print (output)        

def intersect():
    cmd = 'grep -w -v chrom '+input_bed_file+' | awk -v OFS="\t" \'{print $1,$2,$3,$1"_"$2"_"$3,0,$4}\' | ' + path_to_bedtools + ' intersect -a - -b '+output_dir+'/tmp/'+getOutputName2(tmp_bed_annot_file,"sort.bed")

    cmd = cmd + ' -wa -wb'

    if unstranded == False:
        cmd = cmd + " -s"

    cmd = cmd + ' > '+output_dir+"/tmp/"+getOutputName2(input_bed_file,"append.bed")
    
    #print (cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    output = p.stdout.read()
    p.communicate()
    if p.returncode:
        print('{} returned with exit code {}, aborting.'.format(cmd, p.returncode), file=sys.stderr)
        sys.exit(1)

    #print (output)

cluster_to_bed_record = {}

annotated_clusters = {}

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
            s2 = int(tmp[7])
            e2 = int(tmp[8])
            l1 = e2-s1
            l2 = e1-s2
            lm = min(l1,l2)
            l1 = e1-s1
            l2 = e2-s2
            ol = 1.0*lm/(l1+l2-lm)
            if tmp[3] not in cluster_to_bed_record:
                cluster_to_bed_record[tmp[3]]=(line).rstrip('\n')+"\t"+str(ol)
            else:
                if float(cluster_to_bed_record[tmp[3]].rstrip('\n').split("\t")[13]) < ol:
                    cluster_to_bed_record[tmp[3]]=(line).rstrip('\n')+"\t"+str(ol)
    f.close()
    
    outfile = "%s" % output_dir+'/tmp/'+getOutputName2(input_bed_file,"with_jaccard_scores.tsv")
    f = open(outfile,"w")
    for c in cluster_to_bed_record:
        tmp=cluster_to_bed_record[c].split("\t")
        annotated_clusters[tmp[3]]=1
        for x in range(len(tmp)-1):
            f.write("%s\t" % (tmp[x]))
        f.write("%s\n" % tmp[len(tmp)-1])
    f.close()
    #print (len( annotated_clusters))

small_types = {}

#def set_small_types():
#    if small_types_str =="":
#        print ("Please define small RNA types")
#        exit(1)
#    sms=small_types_str.split(",")
#    for x in range(len(sms)):
#        small_types[sms[x]]=1

def is_small(tmp):
    length = int(tmp[2])-int(tmp[1])
    #num=int(tmp[12])
    if length <= 200 and length>=17:
        return True
    else:
        return False    

def print_results():
    #set_small_types()
    small_types = ['misc', 'piRNA', 'miRNA', 'rRNA', 'snRNA', 'snoRNA', 'tRNAs', 'Mt_rRNA', 'Mt_tRNA', 'misc_RNA', 'snRNA',
                   'siRNA', 'vaultRNA','hg19_F_misc','hg19_F_piRNA','hg19_miRNA','hg19_rRna','hg19_snRna','hg19_snoRna','hg19_tRNAs']
    # Append user-supplied biotypes, if available
    if small_types_str:
        small_types += small_types_str.strip().split(',')

    outfile1 = "%s" % output_dir+'/results/'+getOutputName2(input_bed_file,"smallRNAs.tsv")
    outfile2 = "%s" % output_dir+'/results/'+getOutputName2(input_bed_file,"smallRNAs.close.proximity.tsv")    
    outfile3 = "%s" % output_dir+'/results/'+getOutputName2(input_bed_file,"rejected.tsv")

    f1 = open(outfile1,"w")
    f2 = open(outfile2,"w")
    f3 = open(outfile3,"w")
    #f4 = open(outfile4,"w")
    #f5 = open(outfile5,"w")
   
    #f42 = open(outfile42,"w")
    f1.write('chrom\tstart\tend\tstrand\ttotal_mapped_reads\tunique_mapped_reads\tshared_reads\tfeature\tannotation_group\tbest_feature\tbest_feature_JS\n')
    f2.write('chrom\tstart\tend\tstrand\ttotal_mapped_reads\tunique_mapped_reads\tshared_reads\tfeature\tbest_feature\tbest_feature_JS\n') 
    f3.write('chrom\tstart\tend\tstrand\ttotal_mapped_reads\tunique_mapped_reads\tshared_reads\tfeature\tannotation_group\tbest_feature\tbest_feature_JS\n')
 
    ### loop through best annotated ones 
    bestclusterIDs = {}
    bestFeature = {}
    JC_score = {}
    #f=open(input_bed_file,"r")
    f=open(output_dir+"/tmp/"+getOutputName2(input_bed_file,"with_jaccard_scores.tsv"),"r")
    while True:
        line = f.readline()
        if line == "":
            break
        else:
            tmp = line.split("\t")
            bestclusterIDs[tmp[3]]= tmp[12]
            bestFeature[tmp[3]] = tmp[9]
            JC_score[tmp[3]] = tmp[13]
            
    f.close()

    ff=open(input_bed_file,"r")
    while True:
        line = ff.readline()
        line = line.strip() 
        if line == "":
            break
        else:
            tmp = line.split("\t")
            del tmp[8] 
            clusterID = tmp[0]+"_"+str(tmp[1])+"_"+str(tmp[2])

            if clusterID in bestclusterIDs:
                clusterBioType = bestclusterIDs[clusterID]
                best_feature = bestFeature[clusterID]
                f_JC_score = float('%.2f' % float(JC_score[clusterID].rstrip("\n")))
                
                if clusterBioType in small_types and float(f_JC_score) >= j_score and is_small(tmp):
                    f1.write("\t".join(tmp)+"\t"+best_feature+"\t"+str(f_JC_score)+"\n")
                elif tmp[8] == "SM-ONLY" and is_small(tmp):
                    #f2.write("\t".join(tmp[:-1])+"\n")
                    f2.write("\t".join(tmp[:-1])+"\t"+best_feature+"\t"+str(f_JC_score)+"\n")
                else: 
                    f3.write("\t".join(tmp)+"\t"+best_feature+"\t"+str(f_JC_score)+"\n")
                    #f3.write("\t".join(tmp)+"\n")
            else:
                if tmp[8] == "SM-ONLY" and is_small(tmp):
                    #f2.write("\t".join(tmp[:-1])+"\n")
                    f2.write("\t".join(tmp[:-1])+"\t"+best_feature+"\t"+str(f_JC_score)+"\n") 
   
    ff.close()
    f1.close() 
    f2.close()
    f3.close()

def generate_leng_ranges(output):
    ### compute clusters length
    cl_width = []
    with open(output+'/'+getOutputName2(input_bed_file, 'best.tsv'), 'r') as f:
        first_line = f.readline()
        for line in f:
            line = line.rstrip('\n')
            cl_width.append(abs(int(line.split('\t')[1]) - int(line.split('\t')[2])) + 1)

    ### generate length ranges
    max_length = max(cl_width)
    l_bound = []
    u_bound = []
    for i in range(1, max_length, 5):
        l_bound.append(i)
        u_bound.append(i+5-1)

    ### check the last upper-bound values 
    if (abs(u_bound[-1] - max_length) == 1):
        u_bound[-1] = max_length

    cl_ranges = {}
    for x in range(0, len(cl_width)):
       cl_len = cl_width[x]
       for z in range(0, len(l_bound)):
           if cl_len >= l_bound[z] and cl_len <= u_bound[z]:
               cl_ranges[cl_len] = str(l_bound[z])+'-'+str(u_bound[z])

    ### extract features 
    #df = pd.DataFrame(columns=['Feature', 'Biotype', 'Cluster_Length'])
    outfile = open(output+'/tmp/features.tsv', 'w+')
    outfile.write('Length_range\tCluster_Length\tBiotype\n')
    #pos = 0
    with open(output+'/'+getOutputName2(input_bed_file, 'best.tsv'), 'r') as inFile:
        first_line = inFile.readline()
        for line in inFile:
            line = line.rstrip('\n')
            cluster_width = abs(int(line.split('\t')[1]) - int(line.split('\t')[2])) + 1
            features = line.split('\t')[7].split(',')
            cl_len_range = cl_ranges[cluster_width]

            f_biotype = []
            for x in range(0, len(features)):
                f_bio = features[x].split("(")[-1].replace(")", "")
                f_bio = f_bio.replace("hg19_F_", "")
                f_bio = f_bio.replace("hg19_", "")
                f_bio = f_bio.replace("misc_RNA", "misc")
                f_bio = f_bio.replace("Mt_", "")
                f_bio = f_bio.replace("tRNAs", "tRNA")
                f_bio = f_bio.replace("Rna", "RNA")

                if len(f_biotype) ==0:
                    outfile.write(str(cl_len_range)+'\t'+str(cluster_width)+'\t'+f_bio+'\n')
                    f_biotype.append(f_bio.upper())
                elif f_bio.upper() not in f_biotype:
                    outfile.write(str(cl_len_range)+'\t'+str(cluster_width)+'\t'+f_bio+'\n')
                    f_biotype.append(f_bio.upper())

    outfile.close()

# def plot_annot_length_dist(output):
#     ### generate length ranges 
#     generate_leng_ranges(output)

#     ff = pd.read_csv(output+'/tmp/features.tsv', sep='\t')
#     ff = ff.sort_values(by=['Cluster_Length'])
#     ff_cts = ff.groupby(['Length_range', 'Biotype'], sort=False).agg('size').reset_index()
#     ff_cts.columns = ['Length_range', 'Biotype', 'Counts']
#     ff_cts2 = ff_cts.pivot('Biotype','Length_range','Counts')

#     ### reorder columns 
#     column_order = ff_cts['Length_range'].tolist()
#     cols = []
#     for x in column_order:
#         if x not in cols:
#             cols.append(x)
#     ff_cts2 = ff_cts2.reindex(cols, axis=1)

#     sns.set_style("ticks")
#     #fig = plt.figure()
#     fig, ax = plt.subplots(1,1, figsize=(20, 8))
#     sns.heatmap(data=ff_cts2,
#                 cmap=sns.color_palette("YlGnBu"),
#                 annot=True, annot_kws={"size": 7},
#                 cbar_kws={"shrink": .45, 'label': 'Number of small RNAs'},
#                 square=True, vmin = 0, vmax= ff_cts['Counts'].max(),
#                 fmt=".0f", ax=ax)

#     ax.set_xticklabels(ax.get_xticklabels(), rotation =90, fontsize=10)
#     ax.set_yticklabels(ax.get_yticklabels(), rotation =0, fontsize=10)
#     #ax.set_title("Length distribution of annotated small RNAs (n="+str(len(annot_res))+")", fontsize=15)
#     ax.set_xlabel('Length ranges')
#     ax.set_ylabel('')
#     plt.xticks(fontsize=10, ha="center")
#     plt.yticks(fontsize=10, verticalalignment="center")

#     for tick in ax.get_xticklabels():
#         tick.set_rotation(90)
#     for _, spine in ax.spines.items():
#         spine.set_visible(True)

#     plt.savefig(output+'/annotated_smallRNAs_length_dist.png', dpi=400)
 
def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_file=='' or gtf_file=='':
        usage()
        exit(0);
    else:
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        print_records()
        intersect() 
        get_best_annot()
        print_results()
        #plot_annot_length_dist(output_dir)
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))






