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
[bedtools2](https://github.com/arq5x/bedtools2)
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) 

### Installation
Clone the DANST repository (with PATH_TO_DANSR is the local directory where you want to install DANSR):

```
cd PATH_TO_SV_DANSR
git clone https://github.com/ChrisMaherLab/DANSR.git
cd DANSR
```
You also need to add the path to "src" directory to your PATH variable in case you do not want to speificy the who. 

#### Setup
Please note that you need to update the setup.small.ini available at the "src" directory. This file include the paths to the tool DANSR uses. 

