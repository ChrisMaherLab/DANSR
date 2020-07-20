#!/usr/bin/python
from __future__ import division
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename
from difflib import get_close_matches
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns; sns.set()

def usage():
    print("""
    filter_result -i <bedfile>
    
    Parameters:
        -i/--input-bed-file    [string  :  path to BED with reads                                                            ]
        -o/--output-dir        [string  :  path to output dir.    Default: ./                                                ]
        -c/--cutoff            [int:    :  min number of reads for a clusters. Default: 2                                    ]
        -p/--uniq-percent      [float   :  min percentage of unique reads for good-quality clusters. Default 50%             ]
        -u/--min-unique        [int:    :  min number of unique reads for good-quality clusters. Default: 2                  ]
        -v/--ov-with-largest   [float   :  min percentage of overlap with largest node for low-quality clusters. Default 75% ]
    Version:                    1.0.0
          """)


# parameters
input_bed_file = ''
output_dir = ''
cutoff = 2
uniq_perc = 0.50
min_uniq = 2
ov_with_largest = 0.75

def setDefault():
    global output_dir
    output_dir = './'
    global uniq_perc
    uniq_perc = 0.50
    global min_uniq
    min_uniq = 2
    global cutoff
    cutoff = 2
    global ov_with_largest
    ov_with_largest = 0.75
    global is_rm_tmp
    is_rm_tmp = True


def use_real_path():
    global output_dir
    output_dir = os.path.realpath(output_dir)


def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:o:c:p:u:v:", ["help",
                                                           "input-bed-file=",
                                                           "output-dir=",
                                                           "cutoff=",
                                                           "uniq-percent=",
                                                           "min-unique="
                                                           "ov-with-largest="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(1)
        #elif opt in ("-k", "--keep-tmp"):
        #    global is_rm_tmp
        #    is_rm_tmp = False
        elif opt in ("-i", "--input-bed-file"):
            global input_bed_file
            input_bed_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-c", "--cutoff"):
            global cutoff
            cutoff = int(arg)
        elif opt in ("-p", "-uniq-percent"):
            global uniq_perc
            uniq_perc = float(arg)
        elif opt in ("-u", "--min-unique"):
            global min_uniq
            min_uniq = int(arg)
        elif opt in ("-v", "--ov-with-largest"):
            global ov_with_largest
            ov_with_largest = float(arg)


def make_dir(path):
    if not os.path.exists(path):
        os.mkdir(path, 0o755)


def getFilename(filename):
    bn = basename(filename)
    tmp = bn.split(".")
    fn = tmp[0]
    for x in range(1, len(tmp)-1):
        fn = fn+"."+tmp[x]
    return fn


def getOutputName2(filename, append_format):
    fn = getFilename(filename)
    return fn+"."+append_format

#def set_small_types():
#    if small_types_str =="":
#        print ("Please define small RNA types")
#        exit(1)
#    sms=small_types_str.split(",")
#    for x in range(len(sms)):
#        small_types[sms[x]]=1

def determine_feature_class(f):

    #set_small_types() 
   
    smallRNA = ['misc', 'piRNA', 'miRNA', 'rRNA', 'snRNA', 'snoRNA', 'tRNAs', 'Mt_rRNA', 'Mt_tRNA', 'misc_RNA', 'snRNA',
                'siRNA', 'vaultRNA','hg19_F_misc','hg19_F_piRNA','hg19_miRNA','hg19_rRna','hg19_snRna','hg19_snoRna','hg19_tRNAs']

    protein_coding = ["protein_coding", "polymorphic_pseudogene", "TR_V_gene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene",
                      "TR_C_gene", "IG_D_gene", "TR_D_gene", "IG_LV_gene", "IG_M_gene", "IG_Z_gene", "nonsense_mediated_decay",
                      "non_stop_decay", "nontranslating_CDS", "TR_gene"]

    pseudogene = ["IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "polymorphic_pseudogene", "processed_pseudogene",
                  "pseudogene", "disrupted_domain", "IG_pseudogene", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "unitary_pseudogene",
                  "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "unprocessed_pseudogene", "transcribed_unitary_pseudogene"]

    lncrna = ["3prime_overlapping_ncrna", "antisense", "lincRNA", "processed_transcript","sense_intronic", "sense_overlapping", "non_coding", "retained_intron", 
              "hg19_lincRNAsTranscripts"]

    sm_ctr = pc_ctr = ps_ctr = ln_ctr = all_ctr = 0

    for x in range(len(f)):
        #feature = f[x].replace("hg19_F_", "")
        #feature = feature.replace("hg19_", "")
        feature = f[x].split("(")[-1].replace(")", "")

        sm = feature in smallRNA
        pc = feature in protein_coding
        ps = feature in pseudogene
        ln = feature in lncrna

        if sm is True:
            sm_ctr += 1
        if pc is True:
            pc_ctr += 1
        if ps is True:
            ps_ctr += 1
        if ln is True:
            ln_ctr += 1
    
    if sm_ctr > 0:
        if pc_ctr == 0 and ps_ctr == 0 and ln_ctr == 0:
            cl_class = 'SM-ONLY'
        elif pc_ctr > 0 and ps_ctr == 0 and ln_ctr == 0:
            cl_class = 'SM-PC'
        elif pc_ctr == 0 and ps_ctr > 0 and ln_ctr == 0:
            cl_class = 'SM-PS'
        elif pc_ctr == 0 and ps_ctr == 0 and ln_ctr > 0:
            cl_class = 'SM-LN'
        elif pc_ctr > 0 and ps_ctr > 0 and ln_ctr == 0:
            cl_class = "SM-PC-PS"
        elif pc_ctr > 0 and ps_ctr == 0 and ln_ctr > 0:
            cl_class = "SM-PC-LN"
        elif pc_ctr == 0 and ps_ctr > 0 and ln_ctr > 0:
            cl_class = "SM-PS-LN"
        elif pc_ctr > 0 and ps_ctr > 0 and ln_ctr > 0:
            cl_class = "SM-PC-PS-LN"
        else:
            cl_class = "NC"
    else:
        cl_class = 'NO-SM'

    return (cl_class)


def decision_tree(num_nodes, cluster):
    
    nodes, t_reads, u_reads, ovl, ovany, has_annot = [], [], [], [], [], []
    annot, cl_group = [], []
   
    for x in range(len(cluster)):
        if x == 0:
            nodes.append(cluster[x].split(";")[0].split(":")[1])
        else:
            nodes.append(cluster[x].split(";")[0])

        # n_degree.append(int(cluster[x].split(";")[1].replace('degree:','')))
        t_reads.append(int(cluster[x].split(";")[2].replace('alignments:', '')))
        u_reads.append(int(cluster[x].split(";")[3].replace('uniq:', '')))
        ovl.append(int(cluster[x].split(";")[4].replace('overlap:', '')))
        ovany.append(int(cluster[x].split(";")[2].replace('alignments:', '')) - int(cluster[x].split(";")[3].replace('uniq:', '')))
        if (len(cluster[x].split(";")) == 6):
            has_annot.append("Y")
            annot.append(cluster[x].split(";")[5])
            cl_group.append(determine_feature_class(cluster[x].split(";")[5].split(",")))
        else:
            has_annot.append("N")
            annot.append("NA")
            cl_group.append("Novel")

    #determine largest node(s)
    maxValue = max(t_reads)
    idx = [i for i, j in enumerate(t_reads) if j == maxValue]
      
    #define response variable
    resp = [0] * int(num_nodes)

    #(1) is largest has annotation
    for i in idx:
        if has_annot[i] == 'Y':
            resp[i] = 'A'
        else:
            if t_reads[i] > cutoff:
                resp[i] = 'U'
            else:
            	resp[i] = 'L'

    #(2) loop through other nodes and determine the status
    for i in range(0, int(num_nodes)):
        if i in idx:
            continue

        if has_annot[i] == "Y":
            resp[i] = "A"
        else:
            if t_reads[i] <= cutoff:
                resp[i] = "L"
            else:
                if ovany[i] == 0:
                    if u_reads[i] > min_uniq:
                        resp[i] = "U"
                    else:
                        resp[i] = "L"
                else:
                    if ovl[i]/t_reads[i] > ov_with_largest and u_reads[i] ==0:
                        resp[i] = "L"
                    else:
                        if u_reads[i]/t_reads[i] > uniq_perc:
                            resp[i] = "U"
                        else:
                            if ovany[i]/t_reads[i] > 0.50:
                                if u_reads[i] > min_uniq:
                                    resp[i] = "U"
                                else:
                                    resp[i] = "L"
                            else:
                                resp[i] = "U"

    # write results
    for z in range(len(nodes)):
        tmp = nodes[z].rstrip("\n").split("_")
        c_chr = "chr" + tmp[0]
        c_start = int(tmp[1]) - 1
        c_end = int(tmp[1]) - 1 + int(tmp[2])
        if tmp[3] == "0":
            c_strand = "+"
        else:
            c_strand = "-"

        f2.write(c_chr+"\t"+str(c_start)+"\t"+str(c_end)+"\t"+c_strand+"\t" +
                 str(t_reads[z])+"\t"+str(u_reads[z])+"\t" + str(ovany[z])+"\t" +
                 annot[z]+"\t"+resp[z]+"\t"+cl_group[z]+"\n")


def output_results(output_dir):
    exists = os.path.isfile(output_dir+'/results/'+getOutputName2(input_bed_file, 'all.clusters.tsv'))
    if exists:
        cl_file = pd.read_csv(output_dir+'/results/'+getOutputName2(input_bed_file, 'all.clusters.tsv'), sep='\t')
        global annot_res, novel_res
        annot_res= cl_file.loc[~cl_file['annotation_group'].isin(['NO-SM','Novel'])]
        #annot_res = annot_res.drop('class', axis=1)
        #annot_res = annot_res.drop('annotation_group', axis=1) 
        
        #annot_ov_pc = cl_file.loc[cl_file['annotation_group'].isin(['SM-PC','SM-LN','SM-PS','SM-PC-LN','SM-PC-PS','SM-PS-LN','SM-PC-PS-LN']) ]
        #annot_ov_pc = annot_ov_pc.drop('class', axis=1)
        
        novel = cl_file.loc[cl_file['class'] == 'U']       
        novel = novel.drop(['class','feature','annotation_group'], axis=1)
        #novel = novel.drop('feature', axis=1)
        novel_small = novel.loc[(novel['end'] - novel['start'] <= 200) & (novel['end'] - novel['start'] >=17)  ] 
        novel_other = novel.loc[(novel['end'] - novel['start'] > 200) | (novel['end'] - novel['start'] < 17)  ]

        annot_res.to_csv(output_dir+'/tmp/'+getOutputName2(input_bed_file, 'annotated.tsv'), sep='\t', index=False)
        #annot_ov_pc.to_csv(output_dir+'/'+getOutputName2(input_bed_file, 'smallRNAs.overlap.with.other.features.tsv'), sep='\t', index=False)
        novel_other.to_csv(output_dir+'/results/'+getOutputName2(input_bed_file, 'unannotated.rejected.tsv'), sep='\t', index=False, na_rep="NA")
        novel_small.to_csv(output_dir+'/results/'+getOutputName2(input_bed_file, 'unannotated.smallRNAs.tsv'), sep='\t', index=False, na_rep="NA")

    else:
        sys.exit('File '+output_dir+'/'+getOutputName2(input_bed_file,'all.clusters.tsv')+' was not found!.')

# def plot_clusters_width_dist(output):
#     sns.set_style("whitegrid")
#     for file in ['annot', 'novel']:
#         if file == 'annot':
#             tmp = pd.read_csv(output+'/'+getOutputName2(input_bed_file, 'annotated.smallRNAs.tsv'), sep='\t') 
#             mytitle = "Distribution of the length of annotated small RNA clusters (n="+str(len(annot_res))+")"
#             fileName = "annotated_smallRNA_clusters_width_dist.png"
#         else: 
#             tmp = pd.read_csv(output+'/'+getOutputName2(input_bed_file, 'unannotated.smallRNAs.tsv'), sep='\t')
#             mytitle = "Distribution of the length of unannotated small RNA clusters(n="+str(len(novel_res))+")"
#             fileName = "unannotated_smallRNA_clusters_width_dist.png"
        
#         tmp['length'] = (tmp['end'] - tmp['start']) + 1
#         tmp['group'] = '> 500'
#         tmp.loc[(tmp['length'] <= 500) & (tmp['length'] >200) , 'group'] = '<=500 & >200'
#         tmp.loc[(tmp['length'] <= 200) & (tmp['length'] >100) , 'group'] = '<=200 & >100'
#         tmp.loc[(tmp['length'] <= 100) & (tmp['length'] >50) , 'group'] = '<=100 & >50'
#         tmp.loc[(tmp['length'] <= 50) & (tmp['length'] >20), 'group'] = '<=50 & >20'
#         tmp.loc[(tmp['length'] <= 20), 'group'] = '<=20'
        
#         cl_counts = tmp.groupby('group', sort=False).agg('size').reset_index()
#         cl_counts.columns = ['group', 'num_clusters']
        
#         ### create breaks
#         brks = ["<=20", "<=50 & >20","<=100 & >50","<=200 & >100","<=500 & >200","> 500"]
#         for x in brks:
#             if x not in cl_counts['group'].tolist():
#                 brks.remove(x)
        
#         grps = []
#         for x in brks:
#             grps.append(cl_counts[cl_counts.group==x].num_clusters.item())
        
#         fig, ax = plt.subplots(1,1, figsize=(1.6* len(cl_counts['group']),5))
#         ax = sns.barplot(x="group", y = "num_clusters", data=cl_counts, order=brks)            
#         ax.set_title(mytitle, fontsize=15)
#         ax.set_xlabel('')
#         ax.set_ylabel('Number of clusters', fontsize=12)
#         ax.set(yticks= np.arange(0, cl_counts['num_clusters'].max(), 100))
#         for index, row in cl_counts.iterrows():
#             ax.text(row.name,grps[index], grps[index], color='black', ha="center")
    
#         plt.savefig(output+'/'+fileName)

def generate_leng_ranges(output):
    ### compute clusters length
    cl_width = [] 
    with open(output+'/'+getOutputName2(input_bed_file, 'annotated.smallRNAs.tsv'), 'r') as f:
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
    with open(output+'/'+getOutputName2(input_bed_file, 'annotated.smallRNAs.tsv'), 'r') as inFile:
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
#     ax.set_title("Length distribution of annotated small RNAs (n="+str(len(annot_res))+")", fontsize=15)
#     ax.set_xlabel('Length ranges')
#     ax.set_ylabel('')
#     plt.xticks(fontsize=10, ha="center")
#     plt.yticks(fontsize=10, verticalalignment="center")
    
#     for tick in ax.get_xticklabels():
#         tick.set_rotation(90)
#     for _, spine in ax.spines.items():
#         spine.set_visible(True)
    
#     plt.savefig(output+'/annotated_smallRNAs_length_dist.png', dpi=400)
        
def run_classification():
    f = open(input_bed_file, "r")
    #lines = f.readline()
    outfile2 = "%s" % output_dir+'/results/' + getOutputName2(input_bed_file, "all.clusters.tsv")
    global f2
    f2 = open(outfile2, "w")
    f2.write("chrom\tstart\tend\tstrand\ttotal_mapped_reads\tunique_mapped_reads\tshared_reads\tfeature\tclass\tannotation_group\n")
    ct = 1
    for line in f:
        if line == "":
            break
        else:
            #print ("Running for cluster: " + str(ct))
            group_nodes = line.split("\t")
            num_nodes = group_nodes[len(group_nodes)-1].split(":")[0]

            if num_nodes.isdigit() is True:
                new_line = group_nodes[len(group_nodes)-1].rstrip('\n').split("|")
                decision_tree(num_nodes, new_line)
            else:
                continue

        # if ct == 3:
        #	break
        ct = ct + 1

    f.close()
    f2.close()

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_file == '':
        usage()
        exit(1)
    else:
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        run_classification()
        # filter Low-Quality clsuters out
        output_results(output_dir)
        #plot_clusters_width_dist(output_dir)
        #plot_annot_length_dist(output_dir)
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
