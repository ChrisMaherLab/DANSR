#!/usr/bin/python
from __future__ import division
import sys
import os
import math
import getopt
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE, STDOUT, call
from os.path import basename
from difflib import get_close_matches

def usage():
    print ("""
    merge_clusters.py -i <comma separated bed_files>
    
    Parameters:
        -i/--sample-info-file    [string  :    path to samples information file             				]
        -o/--output-dir          [string  :    path to output dir.    Default: ./           				]
        -s/--sample-name         [string  :    sample name            Default: sample       				]
        -b/--bedtools            [string  :    path to bedtoos.       Default: bedtools     				]
        -d/max-dist              [integer :    Maximum distance between features to be merged. Default 0	]
        -k/--keep-tmp            [        :    keep the tmp folder.   Default: not keep     				]
    
    Version:                    1.0.0
          """)
#parameters
sample_info_file = ''
output_dir = ''
is_rm_tmp=None
path_to_bedtools = ''
sample_name = ''
dist = 0

def setDefault():
    global output_dir
    output_dir = './'
    global is_rm_tmp
    is_rm_tmp=True
    global path_to_bedtools
    path_to_bedtools = 'bedtools'
    global sample_name
    sample_name = "sample"
    global dist
    dist = 0
	
def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hki:b:o:s:d:",["help",
                                                      "keep-tmp",
                                                      "sample-info-file=",
                                                      "bedtools=",
                                                      "output-dir=",
                                                      "sample-name=",
                                                      "max-dist",])
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
            global sample_info_file
            sample_info_file = arg
        elif opt in ("-o", "--output-dir"):
            global output_dir
            output_dir = arg
        elif opt in ("-b", "--bedtools"):
            global path_to_bedtools
            path_to_bedtools = arg
        elif opt in ("-s", "--sample-name"):
            global sample_name
            sample_name = arg
        elif opt in ("-d", "--max-dist"):
            global dist
            dist = arg

def make_dir(path):
    if not os.path.exists(path):
        os.mkdir( path, 0o755 ) 

def write_files(files, cmbFile):
	# For the header 
	with open(files[0], 'r') as first:
		cmbFile.write(first.read())
	
	for i in range(1, len(files)):
		with open(files[i], 'r') as f:
			# ignore the rest of headers 
			next(f, None)
			for line in f:
				cmbFile.write(line) 
 
def combine_files(samplesFile, outputPath):
	# count the n8umber of files 
	global num_files 
	num_files = len(open(samplesFile).readlines())
	
	# read the file contains all files 
	filenames = list()            
	with open(samplesFile,"r") as f:
		for line in f:
			filenames.append(line)
	
	with open(output_dir+"/tmp/combined_clusters.tsv", 'w+') as outfile:
		# First file 
		with open(filenames[0].rstrip("\n"), 'r') as first:
			for line1 in first:
				outfile.write(line1)
					
		# Rrest of files
		for i in range(1, num_files):
			with open(filenames[i].rstrip("\n"), 'r') as infile:
				next(infile, None)
				for line in infile:
					outfile.write(line)
	
def run_merge_clusters(output_dir):
	combined_file = "%s" % (output_dir+"/tmp/combined_clusters.tsv")
	
	mergedFileFeature = "%s" % (output_dir+"/tmp/mergedClustersFeature.tsv")
	mergedFileClass = "%s" % (output_dir+"/tmp/mergedClustersClass.tsv")
	mergedFileAnnotClass = "%s" % (output_dir+"/tmp/mergedClustersAnnotClass.tsv")
	
	cmd1 = "grep -w -v chrom "+ combined_file + " | "+path_to_bedtools+" sort -i - | " + path_to_bedtools + " merge -i - -d "+str(dist)+" -c 8 -o collapse -delim \",\" > " + mergedFileFeature
	cmd2 = "grep -w -v chrom "+ combined_file + " | "+path_to_bedtools+" sort -i - | " + path_to_bedtools + " merge -i - -d "+str(dist)+" -c 9 -o collapse -delim \",\" > " + mergedFileClass
	cmd3 = "grep -w -v chrom "+ combined_file + " | "+path_to_bedtools+" sort -i - | " + path_to_bedtools + " merge -i - -d "+str(dist)+" -c 10 -o collapse -delim \",\" > " + mergedFileAnnotClass
	#print(cmd)
	p1 = call(cmd1, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
	p2 = call(cmd2, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
	p3 = call(cmd3, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
	#output = p1.stdout.read()

def determine_clusterClassMV(output_dir):
	#Function to determine cluster class based on majority vote
	merged_file = "%s" % output_dir+"/tmp/mergedClustersClass.tsv"
	fout = "%s" % output_dir+"/merged_clusters/mergedClusters.tsv"
	outFile=open(fout, "w")
	outFile.write("chrom\tstart\tend\tmerged_classes\tcluster_class\n")
	
	f = open(merged_file,"r")

	for line in f:
		line = line.rstrip("\n")
		cl_class = line.split("\t")[3].split(",")
		
		if cl_class.count("A") > 0:
			outFile.write(line+"\t"+"A\n")
		else:
			other_class_cts = {"U": cl_class.count("U"), "L":cl_class.count("L")}
			overall_class = max(other_class_cts, key=other_class_cts.get)
			outFile.write(line+"\t"+overall_class+"\n")
	
	f.close()
	outFile.close()

def determine_clusterClass(output_dir):
	#Function to determine cluster class based on percentage
	merged_file = "%s" % output_dir+"/tmp/mergedClustersClass.tsv"
	fout = "%s" % output_dir+"/tmp/mergedClustersClassProcessed.tsv"
	outFile=open(fout, "w")
	outFile.write("chrom\tstart\tend\tmerged_cluster_classes\tcluster_class\n")
	
	### computer the percentage of un-annotated samples 
	pct_samples=1 if num_files*0.05 < 0.5 else round(num_files*0.05) 
	
	f = open(merged_file,"r")
	for line in f:
		line = line.rstrip("\n")
		cl_class = line.split("\t")[3].split(",")
		
		if cl_class.count("A") > 0:
			outFile.write(line+"\t"+"A\n")
		else:
			other_class_cts = {"U": cl_class.count("U"), "L":cl_class.count("L")}			
			if other_class_cts["U"] > pct_samples: 
				overall_class = 'U'
			else:
				overall_class = 'L'
			
			outFile.write(line+"\t"+overall_class+"\n")
	
	f.close()
	outFile.close()

def determine_clusterAnnotClass(output_dir):
    mergedFile = "%s" % output_dir+"/tmp/mergedClustersAnnotClass.tsv"
    f = open(mergedFile,"r")
    fout = "%s" % output_dir+"/tmp/mergedClustersAnnotClassProcessed.tsv"
    outFile=open(fout, "w")
    outFile.write("chrom\tstart\tend\tmerged_annotation_classes\tannotation_class\n")
	
    for line in f:
        line = line.rstrip()
        chrom = line.split("\t")[0]
        start = line.split("\t")[1]
        end = line.split("\t")[2]
        meregedClasses = line.split("\t")[3]
        annotClass = line.split("\t")[3].split(",")
        annotClass = [x for x in annotClass if x != 'Novel']
        annotClass = ",".join(np.unique(annotClass))
    
        outFile.write(chrom+'\t'+start+'\t'+end+'\t'+meregedClasses+'\t'+annotClass+'\n')
        #break
    
    f.close()
    outFile.close()
    
def determine_feature_class(f):
    smallRNA = ['misc', 'piRNA', 'miRNA', 'rRNA', 'snRNA', 'snoRNA', 'tRNAs', 'Mt_rRNA', 'Mt_tRNA', 'misc_RNA', 'snRNA',
                'siRNA', 'vaultRNA']
    protein_coding = ["protein_coding", "polymorphic_pseudogene", "TR_V_gene", "IG_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene",
                      "TR_C_gene", "IG_D_gene", "TR_D_gene", "IG_LV_gene", "IG_M_gene", "IG_Z_gene", "nonsense_mediated_decay",
                      "non_stop_decay", "nontranslating_CDS", "TR_gene"]
    pseudogene = ["IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "polymorphic_pseudogene", "processed_pseudogene",
                  "pseudogene", "disrupted_domain", "IG_pseudogene", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "unitary_pseudogene",
                  "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "unprocessed_pseudogene", "transcribed_unitary_pseudogene"]
    lncrna = ["3prime_overlapping_ncrna", "antisense", "lincRNA", "processed_transcript",
              "sense_intronic", "sense_overlapping", "non_coding", "retained_intron", "lncRNA"]
    
    sm_ctr = pc_ctr = ps_ctr = ln_ctr = all_ctr = 0
    for x in range(len(f)):
        if f[x] == 'NA':
        	continue 
        	 
        feature = f[x].replace("hg19_F_", "")
        feature = feature.replace("hg19_", "")
        feature = feature.split("(")[1].replace(")", "")
                
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
        else:
            cl_class = "SM-PC-PS-LN"
    else:
        cl_class = 'NO-SM'
    
    return (cl_class)
    
def determine_clusterFeatures(output_dir):
    mergedFile = "%s" % output_dir+"/tmp/mergedClustersFeature.tsv"
    f = open(mergedFile,"r")
    fout = "%s" % output_dir+"/tmp/mergedClustersFeatureProcessed.tsv"
    outFile=open(fout, "w")
    outFile.write("chrom\tstart\tend\tmerged_features\tfeatures\tannotation_class\n")
	
    for line in f:
        line = line.rstrip()
        chrom = line.split("\t")[0]
        start = line.split("\t")[1]
        end = line.split("\t")[2]
        mergedFeatures = line.split("\t")[3]
        ### determine featues class 
        FC = determine_feature_class(mergedFeatures.split(','))
        
        features = line.split("\t")[3].split(",")
        features = [x for x in features if x != 'NA']
        features = ",".join(np.unique(features))
        if features == '':
        	features = 'Unannotated'
        	FC = 'NC'
        	
        outFile.write(chrom+'\t'+start+'\t'+end+'\t'+mergedFeatures+'\t'+features+'\t'+FC+'\n')
        
    f.close()
    outFile.close()

def combine_all_merged_res(output_dir):
    f1 = pd.read_csv(output_dir+'/tmp/mergedClustersClassProcessed.tsv', sep='\t')
    f2 = pd.read_csv(output_dir+'/tmp/mergedClustersFeatureProcessed.tsv', sep='\t')
    #f3 = pd.read_csv(output_dir+'/tmp/mergedClustersAnnotClassProcessed.tsv', sep='\t')
    
    merge_res = pd.merge(f1[['chrom','start','end','cluster_class']], 
                         f2[['chrom','start','end','features', 'annotation_class']], 
                         on=['chrom','start','end'])                  
    #merge_res = pd.merge(merge_res, 
    #                     f3[['chrom','start','end','annotation_class']], 
    #                     on=['chrom','start','end'])

    annot_res= merge_res.loc[merge_res['annotation_class'].isin(['SM-ONLY']) ]
    annot_ov_pc = merge_res.loc[merge_res['annotation_class'].isin(['SM-PC','SM-LN','SM-PS','SM-PC-PS-LN']) ]
    novel_res= merge_res.loc[merge_res['cluster_class'] == 'U']
    
    merge_res.to_csv(output_dir+'/all.clusters.tsv', sep='\t', index=False)
    annot_res.to_csv(output_dir+'/annotated.smallRNAs.tsv', sep='\t', index=False)
    annot_ov_pc.to_csv(output_dir+'/smallRNAs.overlap.with.PCGs.tsv', sep='\t', index=False)
    novel_res.to_csv(output_dir+'/novel_smallRNAs.tsv', sep='\t', index=False)
                     
def main(argv):
    setDefault()
    getParameters(argv[1:])
    if sample_info_file=='':
        usage()
        exit(1);
    else:
        print ('Merging cluster files for all samples...')
        use_real_path()
        #make_dir(output_dir)
        make_dir(output_dir+'/tmp')
        combine_files(sample_info_file, output_dir)
        run_merge_clusters(output_dir)
        determine_clusterClass(output_dir)
        #determine_clusterAnnotClass(output_dir)
        determine_clusterFeatures(output_dir)
        combine_all_merged_res(output_dir)
        print ('done.')
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
