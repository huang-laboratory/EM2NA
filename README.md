# EM2NA
## Overview
EM2NA is a software for automatic nucleic-acid modeling from cryo-EM density maps. 

## Requirements
**Platform**: Linux (Mainly tested on CentOS 7).

**GPU**: A GPU with >10 GB memory is required, advanced GPU like A100 is recommended.

## Web server
The [EM2NA web server](http://huanglab.phys.hust.edu.cn/EM2NA/server.php).

## Installation
#### 0. Install conda

We use conda to manage the required packages, if it is not installed, please refer to https://docs.anaconda.net.cn/miniconda/install/ for installation.

#### 1. Download EM2NA

Download EMProt via github
```
git clone https://github.com/huang-laboratory/EM2NA.git
cd EM2NA
# if you do not have git, download the zipped file
# wget https://github.com/huang-laboratory/EM2NA/archive/refs/heads/main.zip
# unzip main.zip
# cd EM2NA-main
```

#### 2. Create conda environment
```
conda env create -f environment.yml
```
If conda fails, you could install the packages youself. Basically. you can first create an environment named `emprot` by `conda env create -n emprot python=3.10`, then install the packages listed in `environment.yml` using conda or pip.

#### 3. Install LKH
Refer to [LKH-3](http://webhotel4.ruc.dk/~keld/research/LKH-3) for installation of LKH-3.

#### 4. Set Environment Variables in the EM2NA.sh
```vi EM2NA.sh```

##### 4.1. Set "activate" to path of conda activator, for example
```activate="/path/to/your/coonda/bin/activate"```

##### 4.2. Set "EM2NA_env" to name of the python conda virtual environment, for example
```EM2NA_env="em2na"```

##### 4.3. Set "LKH_dir" to path of LKH-3, for example
```LKH_dir="/path/to/LKH-3"```

##### 4.4. Set "EM2NA_home" to path of EM2NA, for example
```EM2NA_home="/path/to/EM2NA"```

## Usage
Running EM2NA is very straight forward with one command like
```
/path/to/EM2NA/EM2NA.sh MAP.mrc OUTPUT_DIR \
    [-g GPUID] \ # GPUID, default is 0
    [--seq SEQ.fa] \ # optinal
    [--natype NA_TYPE] # DNA or RNA or AUTO
```
- Cryo-EM density map and output directory is **required**.
- Sequence(s) are **optional**.
- Input Fasta file SEQ.fa could include multiple (>= 1) sequences.
- If you launch > 1 modeling jobs, the output directory **MUST** be set differently.

The output model (named **output.cif**) will be saved in the specified output directory.

#### 1. Modeling with target sequence(s)
We provide an example for users, download it from our website
```
http://huanglab.phys.hust.edu.cn/EM2NA/6O1D.tgz
tar -zxvf 6O1D.tgz
cd 6O1D
```
Then run EM2NA
```
/path/to/EM2NA/EM2NA.sh emd_0586.map output --seq 6O1D.fasta -g 0
```
The output model is saved in directory *output*, and named *output.cif*.

#### Post refinement
It is recommended to use third-party programs to further refine the model-map fit and model geometries, e.g. using **phenix.real_space_refine**
```
phenix.real_space_refine emd_0586.map ooutput/output.cif resolution=3.4
```
This will typically take several minutes.


## Citation
Tao Li, et al. Automated detection and de novo structure modeling of nucleic acids from cryo-EM maps. *Nature Communications*. 2024
```
@article {EM2NA2024,
	title = {Automated detection and de novo structure modeling of nucleic acids from cryo-EM maps},
	author = {Tao Li, Hong Cao, Jiahua He, Sheng-You Huang},
	journal = {Nature Comminications},
	year = {2024},
	doi = {}
}
```
