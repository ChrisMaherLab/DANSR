# DANSR: A tool for the detection of annotated and novel small RNAs implemented in Python and C++. 

## 
DANSR is developed at [Christopher Maher Lab](http://www.maherlab.com/) at [Washington University in St. Louis](http://www.wustl.edu).
   
## SV-HotSpot Manual
### Prerequisites
Please make sure you have installed the following tools:

[Python3](https://www.python.org/) <br>
[panda](https://pandas.pydata.org/) <br>
[CMake](https://cmake.org/) <br>
[GCC](https://gcc.gnu.org/) <br>
[Matplotlib](http://matplotlib.org/) <br>
[samtools](https://github.com/samtools/samtools) <br>
[bedtools2](https://github.com/arq5x/bedtools2) <br>
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) 

### Installation
Clone the DANSR repository (with PATH_TO_DANSR is the local directory where you want to install DANSR):

```
cd PATH_TO_SV_DANSR
git clone https://github.com/ChrisMaherLab/DANSR.git
cd DANSR
```
As an optional step, you can add the absolute path to the "src" directoy to the PATH variable. 
 
#### Setup
Please note that you need to update the setup.small.ini available at the "src" directory. This file includes the paths to the tool DANSR uses. 

### Test installation
To test the installation, please type in your terminal the following command which shows the usage page of DANSR. 
```
python path_to/src/dansr.py
```

### Running DANSR
Assume PATH_TO_SV_DANSR is the local directory where DANSR was installed, the following command runs DANSR on the test data provided.

```
python3.7 src/small.py --begin-no-trimming \
        -i PATH_TO_SV_HOTSPOT/test_data/CRC1-P.fastq \
        -o OUTPUT_DIR --type single --cutoff 0.33 \
	-v -w -N 10 -U 0.5 -R 2 -V 0.75 -J 0.3 \
        -r PATH_TO_ANNOTATION/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        -l PATH_TO_ANNOTATION/smallRNA_library.gtf,PATH_TO_ANNOTATION/Homo_sapiens.GRCh37.75.gtf \
        -g PATH_TO_ANNOTATION/smallRNA_library.gtf \
```

### Output  
DANSR generates two main folders, one is for the alignment results and the other one is for the detection results. The main file in the aligment folder is the aligned BAM file. The detection results folder contains several files incldung:  
 
Main files:
1. ```annotated.smallRNAs.tsv```: contains all annotated small RNAs that pass all filters. . 
2. ```annotated.smallRNAs.close.proximity.tsv```: contains all annotated small RNAs passing all filters but were annotated with a close proximity to the feature.  
3. ```unannotated.smallRNAs.tsv```: contains all unannotated (novel) small RNAs that pass all filtering processes. . 

<u>Additional files:</u> 

4. ```all.clusters.tsv```: contains all clusters. 
5. ```annotated.others.tsv```: contais annotated small RNAs that didn't pass the filtering processes (e.g. overlapping with protein coding, small RNA not within the range of 17-200nt, small RNAs with loe jacard score (not the best feature assigned)). . 
6. ```unannotated.other.tsv```


## DANSR Docker Instructions
To use DANSR, a docker image has been created and tested on Linux and Mac. To run DANSR, you need to have [Docker](https://docs.docker.com/) installed on your machine. 

### Docker Installation
* Ubuntu: follow [the instructions](https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/) to get Docker CE for Ubuntu.
* Mac: follow [the instructions](https://store.docker.com/editions/community/docker-ce-desktop-mac) to install [the Stable verion of Docker CE](https://download.docker.com/mac/stable/Docker.dmg) on Mac.
<!--- 
* Windows: follow [the instructions](https://docs.docker.com/toolbox/toolbox_install_windows/) to install [Docker Toolbox](https://download.docker.com/win/stable/DockerTool    box.exe) on Windows. 
-->
 
To obtain the latest docker image, run the following in your command line:
 
```
docker pull chrismaherlab/dansr
```
To test the image, run the following command which shows the usage of this tool:
```
docker run chrismaherlab/dansr dansr
```

