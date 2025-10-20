# EM2NA
## Overview
EM2NA is a software for automatic nucleic-acid modeling from cryo-EM density maps. 

## Requirements
**Platform**: Linux (Mainly tested on CentOS 7).

**GPU**: A GPU with >10 GB memory is required, advanced GPU like A100 is recommended.

## Web server
We probide online service for EM2NA, check it [here](http://huanglab.phys.hust.edu.cn/EM2NA/server.php).

## Installation
#### 0. Install conda

We use conda to manage the required packages, if it is not installed, please refer to https://docs.anaconda.net.cn/miniconda/install/ for installation.

#### 1. Download EM2NA

Download EM2NA via github
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
If conda fails, you could install the packages youself. Basically, you can first create an environment named `em2na` by `conda env create -n em2na python=3.8.8`, then install the packages listed in `environment.yml` using conda or pip.

#### 3. Install LKH
EM2NA uses LKH for tracing of NA backbones, please refer to [LKH-3](http://webhotel4.ruc.dk/~keld/research/LKH-3) for installation of LKH-3.

#### 4. Set Environment Variables in the EM2NA.sh
Edit the EM2NA.sh to set environment variables for EM2NA
```
vi EM2NA.sh

# Set "activate" to path of conda activator, for example
activate="/path/to/your/coonda/bin/activate"

# Set "EM2NA_env" to name of the python conda virtual environment, for example
EM2NA_env="em2na"

# Set "LKH_dir" to path of LKH-3, for example
LKH_dir="/path/to/LKH-3"

# Set "EM2NA_home" to path of EM2NA, for example
EM2NA_home="/path/to/EM2NA"
```

#### 5. Compile python package in Fortran
In addition to online packages, the interpolation program `interp3d.f90` should be built as a python package `interp3d` using `f2py`. Users can either build interp3d from source or use our precompiled interp3d.

##### 5.1 Build interp3d from source
Check where your gfortran is and compile using `f2py` (`f2py` will be available once you have installed numpy in `em2na` env)
```
which gfortran
conda activate em2na
f2py -c ./interp3d.f90 -m interp3d --f90exec=/path/to/your/gfortran --f77exec=/path/to/your/gfortran
```

##### 5.2 Use our precompiled interp3d
By default, we already provided a compiled interp3d package with libgfortran.so.3 requirement in EM2NA home directory. We have compiled another 2 versions of interp3d that requires libgfortran.so.4 or libgfortran.so.5 in directory lib_interp3d/. Check the .so support information for all the 3 verions:
```
ldd lib_interp3d/libgfortran*/*
```
Pick one version (libgfortran3/4/5) that finds all .so files and no errors and copy it to EM2NA home directory. For example, if I have libgfortran4 in my system
```
cp lib_interp3d/libgfortran4/interp3d.cpython-38-x86_64-linux-gnu.so .
```


## Usage
Running EM2NA is very straight forward with one command like
```
/path/to/EM2NA/EM2NA.sh MAP.mrc OUTPUT_DIR \
    [-g GPUID] \ # GPUID, default is 0
    [--seq SEQ.fa] \ # optional
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
The output model is saved in directory `output`, and named `output.cif`.

#### Post refinement
It is recommended to use third-party programs to further refine the model-map fit and model geometries, e.g. using **phenix.real_space_refine**
```
phenix.real_space_refine emd_0586.map output/output.cif resolution=3.4
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
