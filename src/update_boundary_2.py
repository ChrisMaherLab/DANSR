#!/usr/bin/python
import sys
import os
import math
import getopt
from subprocess import Popen, PIPE, STDOUT
from os.path import basename

def usage():
    print ("""
    update_boundary -i <bedfile>
    
    Parameters:
        -i/--input-bed-file    [string  :    path to BED with reads                                          ]
        -o/--output-dir        [string  :    path to output dir.    Default: ./                              ]
        -p/--percent           [float   :    percentage of reads at each side to consider                    ]
        -c/--cutoff            [float   :    increase cutoff to change boundary                              ]
        -k/--keep-tmp          [        :    keep the tmp folder.   Default: not keep                        ]
    
    Version:                    1.0.0
          """)
#parameters
input_bed_file = ''
output_dir = ''
#is_rm_tmp=None
percent = 0.3
cutoff = 0.33

def setDefault():
    global output_dir
    output_dir = './'
    global percent
    percent = 0.3
    global cutoff
    cutoff = 0.33
    #global is_rm_tmp
    #is_rm_tmp=True

def use_real_path():
    global output_dir
    output_dir=os.path.realpath(output_dir)

def getParameters(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:p:c:",["help",
                                                      "input-bed-file=",
                                                      "output-dir=",
                                                      "percent=",
                                                      "cutoff="])
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
            cutoff = float(arg)

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

def getNextPosition(p,sigar):
    if 'S' in sigar or 'H' in sigar: ##### not supporting clipped reads
        return -1
    else:
        s=''
        nums=[]
        num_str=''
        for x in range(len(sigar)):
            if sigar[x] in {'M','I','D','S','H'}:
                s=s+sigar[x]
                nums.append(int(num_str))
                num_str=''
            else:
                num_str=num_str+sigar[x]
        totalLen=0
        for x in range(len(s)):
            if s[x]=='M' or s[x]=='D':
                totalLen=totalLen+nums[x]
        return p+totalLen-1

#def remove_tmp():
#    if is_rm_tmp:
#        cmd = 'rm -rf ' + output_dir +'/tmp'

long_reads_sorted = [] # long means add columns of start and length for sorting

def sort_reads(reads):
    long_reads = []
    for x in range(len(reads)):
        tmp=reads[x].split(",")
        start=tmp[3][1:]     ### read start
        cigar=tmp[4]         ### read length
        one_long = [] #[0: full read string, 1: start position, 2: length of reference spanned by read, 3: count after cut, 4: keep this read? (0/1)]
        one_long.append(reads[x])
        one_long.append(int(start))
        one_long.append(getNextPosition(int(start),cigar)-int(start)+1)
        one_long.append(0.0)# for count after cut
        one_long.append(1)  # for wheather to keep
        long_reads.append(one_long)
    global long_reads_sorted
    long_reads_sorted=sorted(long_reads, key = lambda x: (x[1], x[2]))
    #print "long_reads_sorted length=",len(long_reads_sorted)


counter = []

#get array of per-base read depths for this cluster
def count_reads(cs,ce):
    global counter
    counter = [0] * (ce-cs+1)
    #print "cs,ce=",cs,ce
    global long_reads_sorted
    for x in range(len(long_reads_sorted)):
        #print long_reads_sorted[x][0]
        lr=long_reads_sorted[x]
        if lr[4]==0:
            continue
        be=lr[1]
        en=lr[1]+lr[2]-1
        #print "be,en=",be,en
        for y in range(en-be+1):
            #print "be-cs+y=",be-cs+y
            counter[be-cs+y]=counter[be-cs+y]+1

#Get average read depth across each read in cluster:
def cal_density(cs,ce):
    for x in range(len(long_reads_sorted)):
        lr=long_reads_sorted[x]
        be=lr[1]
        en=lr[1]+lr[2]-1
        nc=0
        for y in range(en-be+1):
            nc=nc+counter[be-cs+y]
        density=nc*1.0/(en-be+1)
        long_reads_sorted[x][3]=density

def cut_long_reads(num_rm,cs,ce):
    count_reads(cs,ce)  ### count reads per base 
    cal_density(cs,ce)  ### density for the whole read to represent the weight 
    global long_reads_sorted
    long_reads_sorted = sorted(long_reads_sorted, key = lambda x: (-x[3], x[1], x[2])) # sorted by weight + coord
    for x in range(num_rm):
        long_reads_sorted[len(long_reads_sorted)-1-x][4]=0   # asssign 0 to the low weight reads 
    long_reads_sorted = sorted(long_reads_sorted, key = lambda x: (x[1], x[2])) ### sort again by the coordinates n
    #Remove cut reads (MJI 03/25):
    long_reads_sorted = [r for r in long_reads_sorted if r[4]]
    #print "marked not included",num_rm 

begins = []
ends = []

def get_sub_clusters(cs,ce):
    global begins
    begins=[]
    global ends
    ends=[]
    count_reads(cs,ce)
    #print "counter:"
    #print counter
    status = 0 # for search next begin, 1 for search next end
    if counter[0]==0:
        status=0
    else:
        status=1
        begins.append(cs)
    for x in range(ce-cs+1):
        if x>0:
            if status == 0 and counter[x]>0:
                status = 1
                begins.append(x+cs)
                continue
            if status == 1 and counter[x]==0:
                status = 0
                ends.append(x+cs)
                continue
    if len(begins)==len(ends)+1:
        ends.append(ce)
    #print "begins,ends=",begins,ends

begins0=[]
ends0=[]

### assign mid points between consecutive subclusters. 
def get_begins_ends0():
    global begins0,ends0
    begins0 = []
    ends0 = []
    for x in range(len(begins)):
        s0=long_reads_sorted[0][1]
        if x > 0:
            s0=ends[x-1]+(begins[x]-ends[x-1])/2  
        e0=long_reads_sorted[len(long_reads_sorted)-1][1]+long_reads_sorted[len(long_reads_sorted)-1][2]-1
        if x < len(begins)-1:
            e0=ends[x]+(begins[x+1]-ends[x])/2
        begins0.append(s0)
        ends0.append(e0)
    #print "begins0,ends0=",begins0,ends0
 
be0=[] #positions in long_reads_sorted
en0=[]
be1=[]
en1=[]

def get_be_en_values():
    global be0,en0,be1,en1
    be0=[-1]*len(begins) 
    en0=[-1]*len(begins)
    be1=[-1]*len(begins)
    en1=[-1]*len(begins)
    for x in range(len(long_reads_sorted)): 
        #print "x=",x
        s=long_reads_sorted[x][1]
        e=long_reads_sorted[x][1]+long_reads_sorted[x][2]-1
        #print "s,e=",s,e
        for y in range(len(begins0)):
            if s>=begins0[y] and be0[y]==-1:
                be0[y]=x
                #print "be0[",y,"]=",x
                break
        for y in range(len(begins)):
            if s>=begins[y] and be1[y]==-1:
                be1[y]=x
                #print "be1[",y,"]=",x
                break
    for x in reversed(range(len(long_reads_sorted))):
        #print "x=",x
        s=long_reads_sorted[x][1]
        e=long_reads_sorted[x][1]+long_reads_sorted[x][2]-1
        #print "s,e=",s,e
        for y in reversed(range(len(ends0))):
            if e<=ends0[y] and en0[y]==-1:
                en0[y]=x
                #print "en0[",y,"]=",x
                break
        for y in reversed(range(len(ends))):
            if e<=ends[y] and en1[y]==-1:
                en1[y]=x
                #print "en1[",y,"]=",x
                break
    #print "be0,en0,be1,en1=",be0,en0,be1,en1 

#Average coverage depth of subcluster
def cal_bpl(ps,pe):#i in begins
    long_reads=long_reads_sorted[ps:pe+1]
    s=long_reads[0][1]
    e=long_reads[len(long_reads)-1][1]+long_reads[len(long_reads)-1][2]-1
    t=0
    for x in range(len(long_reads)):
        t=t+long_reads[x][2]
    return t*1.0/(e-s+1)
 
def evaluate_sub_clusters():
    better_be = []
    better_en = []
    better_be_p = []
    better_en_p = []
    better = [0]*len(begins)
    for x in range(len(better)):
        cal0=cal_bpl(be0[x],en0[x])   #### compute the weight over length for the sub cluster with status 0 cluster including 1+ 
        cal1=cal_bpl(be1[x],en1[x])   ###  compute the weight over length for the sub cluster with status 1  
        #print cal1, cal0
        #print 1.0+cutoff
        if(cal1>cal0*(1.0+cutoff)):
            better[x]=1
    merged = 0
    for x in range(len(better)):
        if (x==len(better)-1 and better[x]==0) or (x+1<len(better) and better[x]==0 and better[x+1]==1):
            better_be.append(begins0[x-merged])
            better_en.append(ends0[x])
            better_be_p.append(be0[x-merged])
            better_en_p.append(en0[x])
            merged=0
            continue
        if x+1<len(better) and better[x]==0 and better[x+1]==0:
            merged=merged+1
            continue
        if better[x]==1:
            better_be.append(begins[x])
            better_en.append(ends[x])
            better_be_p.append(be1[x])
            better_en_p.append(en1[x])
            continue
    #print "better_be,better_en,better_be_p,better_en_p=",better_be,better_en,better_be_p,better_en_p
    return    better_be,better_en,better_be_p,better_en_p

def get_reads_str(ps,pe):
    long_reads=long_reads_sorted[ps:pe+1]
    reads_str = ""
    for x in range(len(long_reads)):
        reads_str=reads_str+long_reads[x][0]
        if x!=len(long_reads)-1:
            reads_str=reads_str+";"
    return reads_str

def get_update_line(line):
    tmp=line.split("\t")
    chr0 = tmp[0]
    start0 = int(tmp[1])
    end0 = int(tmp[2])
    reads = tmp[3].rstrip('\n').split(";")
    num_reads=len(reads)
    num_rm = int(num_reads*percent)
    if num_rm > 0 and num_reads-num_rm > 0:
        sort_reads(reads)
        cut_long_reads(num_rm,start0+1,end0)
        get_sub_clusters(start0+1,end0) #Find starts and ends of subclusters (begins[], ends[])
        get_begins_ends0() #Set boundary b/w each subcluster to midpoint b/w end of one and start of next (begins0[],ends0[])
        get_be_en_values() #Index of the first or last read that falls in each subcluster (be0, en0 -> begins0, ends0), (be1, en1 -> begins, ends)
        be,en,bep,enp=evaluate_sub_clusters()
        lines = ""
        for x in range(len(be)):
            reads_str = get_reads_str(bep[x],enp[x])
            lines = lines+chr0+"\t"+str(be[x]-1)+"\t"+str(en[x])+"\t"+reads_str+"\n"
        #print "updated:"
        #print lines
        return lines
    #print "no update:"
    #print line 
    return line

def run_update():
    infile = "%s" % input_bed_file
    f=open(infile,"r")
    outfile2 = "%s" % output_dir+'/tmp/'+getOutputName2(input_bed_file,"update.bed")
    f2=open(outfile2, "w")
    while True:
        line=f.readline()
        if line=="":
            break
        else:
            f2.write(get_update_line(line))
    f.close()
    f2.close()

def main(argv):
    setDefault()
    getParameters(argv[1:])
    if input_bed_file=='':
        usage()
        exit(1);
    else:
        use_real_path()
        #make_dir(output_dir)
        #make_dir(output_dir+'/tmp')
        run_update()
        #remove_tmp()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
