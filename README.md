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

#### Setup
Please note that you need to update the setup.small.ini available at the "src" directory. This file includes the paths to the tool DANSR uses. 

### Test installation
To test the installation, please type in your terminal the following command which shows the usage page of DANSR. 
```
python src/dansr.py
```

### Running DANSR

