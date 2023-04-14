# AutoQii2

The AutoQii2 Pipeline can be used to analysisng 16S amplicon based datasets. AutoQii2 is primarily designed for eliminating multi-step analysis involved in analyzing paired-end or single-end reads using QIIME2. AutoQii2 uses popular fastQC, cutadapt and QIIME2 platforms for performing quality check, adapter and chimera removal, operational taxonomic units (OTUs) identification and taxonomic assignments.

## AutoQii2 Pipeline (Workflow)

 ![Figure 2_R1](https://user-images.githubusercontent.com/60095100/228609360-39e00bb3-e2fc-4db2-9a4d-7bc72aa8d213.jpg)

## Clone this GitLab Repository

```
git clone https://gitlab.com/khatriabhi2319/autoqii2
cd autoqii2

```

## Requirements

### Install Dependencies

* Conda: 

Conda can be downloaded as part of the Anaconda or the Miniconda plattforms (Python 3.7). We recommend to install miniconda3. Using Linux you can get it with:

```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```
install qiime in conda environment.
All other dependencies will be automatically installed using conda environments and can be found in the corresponding environment.yaml files in the envs folder and the natrix.yaml file in the root directory of the pipeline.


### Softwares requirement


#### 1. FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```
after installation check with command:
$ fastqc --help
````

#### 2. Cutadapt (https://cutadapt.readthedocs.io/en/stable/installation.html)
```users/auth/google_oauth2
after installation check with command:
$ cutadapt --help
```

#### 3. qiime2

Qiime 2 can be installed natively or using virtual machines. For this pipeline we prefer to install with native method
(https://docs.qiime2.org/2021.4/install/native)

#### 4. Picrust 2

For Functional analysis user need to install picrust in qiime 2 envirnment
```
#Install q2-picrust2 with conda. This command will automatically install the other requirements, including PICRUSt2. 
*Note that the plugin version for qiime2-2021.2 is specified

$ conda install q2-picrust2=2021.2 -c conda-forge -c bioconda -c gavinmdouglas

*"2021.2" need to be edit according to your qiime 2 version
```
##### OR 
```
conda install -c bioconda -c conda-forge picrust2
pip install -e .
qiime dev refresh-cache
qiime picrust2 --help
```

#### 5. Greengenes File
For Taxonomic analysis you need to Download Greengenes file according to your qiime2 version.
Open the ([link](https://docs.qiime2.org/2021.4/tutorials/moving-pictures/)) and download the greengene file from Taxonomic analyis section. Please ensure to change the version of qiime2 in the website, similar to your current qiime2 version (Select greengene file from qiime2 official ([link](https://docs.qiime2.org/2021.4/tutorials/moving-pictures/)) According to your qiime version). 
 After download, move the file to/your/working/directory/16s-analysis-main



## Pipeline summary

After Successful compilation of all Reqirements, the pipeline currently performs the following:

* Sequencing quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
* Trimming of reads ([Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200))
* Taxonomical classification using DADA2 with [QIIME2](https://www.nature.com/articles/s41587-019-0209-9)
* Excludes unwanted taxa, produces absolute and relative feature/taxa count tables and plots, plots alpha rarefaction curves, computes alpha and beta diversity indices and plots thereof ([QIIME2](https://www.nature.com/articles/s41587-019-0209-9))


## Steps for Paired End Sequencing

### Step-1

Extract Single-End raw data and place it in Paired-End/raw_data directory.

```
$ tar -xvzf raw_data.tar.gz

```
### Step-2

Edit Paired end "metadata.tsv" file, according to your raw reads filename.

```
# Paired end metadata file

sample-id       forward-absolute-filepath       reverse-absolute-filepath       Type
#q2:types       categorical     categorical     categorical
BPP01_1 $PWD/Paired-End-Analysis/raw_data/HAM29_R1.fastq        $PWD/Paired-End-Analysis/raw_data/HAM29_R2.fastq        BPP 01
BPP01_2 $PWD/Paired-End-Analysis/raw_data/HAM30_R1.fastq        $PWD/Paired-End-Analysis/raw_data/HAM30_R2.fastq        BPP 01
BPG01_1 $PWD/Paired-End-Analysis/raw_data/ERR4362125_1.fastq    $PWD/Paired-End-Analysis/raw_data/ERR4362125_2.fastq    BPG 01
BPG01_2 $PWD/Paired-End-Analysis/raw_data/ERR4362126_1.fastq    $PWD/Paired-End-Analysis/raw_data/ERR4362126_2.fastq    BPG 01

```
(https://docs.qiime2.org/2019.7/tutorials/importing/)

### Step-3

Activate Qiime2 through conda environment

```
$ conda activate qiime2-2019.10

```
### Step-4

Run "16S_AutoQii2.sh" Script.

```
*You need to make the shell script executable using chmod command

$chmod +x 16S_AutoQii2.sh

$ ./16S_AutoQii2.sh

It will show GUI Interface (dialogbox) for your confirmation, ask for raw data directory and metadata file selection :


```

## Steps for Single End Sequencing

### Step-1

Extract Single-End raw data and place it in Single-End/raw_data directory.

```
tar -xvzf raw_data_SE.tar.gz

```
### Step-2

Edit Single-End "metadata.tsv" file, according to your raw reads filename.

```
sample-id       absolute-filepath       Type
#q2:types       categorical     categorical
BPP01_1 $PWD/Single-End-Analysis/raw_data/HAM29_R1.fastq        BPP 01
BPP01_2 $PWD/Single-End-Analysis/raw_data/HAM30_R1.fastq        BPP 01
BPG01_1 $PWD/Single-End-Analysis/raw_data/ERR4362125_1.fastq    BPG 01
BPG01_2 $PWD/Single-End-Analysis/raw_data/ERR4362126_1.fastq    BPG 01

```
### Step-3

Activate Qiime2 through conda environment

```
conda activate qiime2-2019.10

```
### Step4

Run "16S_AutoQii2.sh" Script.

```
*You need to make the shell script executable using chmod command

$chmod +x 16S_AutoQii2.sh

$ ./16S_AutoQii2.sh

It will show GUI Interface (Zenity dialogbox) for your confirmation, ask for raw data directory and metadata file selection :

```

## Examples, How to run AutoQii2 pipeline:

```
Download "AutoQii2" on your system using command

$git clone https://gitlab.com/khatriabhi2319/autoqii2.git

```

### What to do next and how to select parameters ? 

#####  1. A welcome Message:

![oie_BUPGzoM9vJWy](https://user-images.githubusercontent.com/60095100/191945730-ca665efd-3549-4991-835c-5cb0cddc7610.png)

#####  2. Raw data directory selection:

![raw-data](https://user-images.githubusercontent.com/60095100/192094339-b9a64dc0-8287-4e24-8cd0-f91617ced2a7.png)

#####  3. Information regarding metadata file preparation:

![metadata_preparation](https://user-images.githubusercontent.com/60095100/190382347-752a8e81-7969-445a-a034-db939daf9dd0.png)

#####  4. Select metadatafile:

![metadata](https://user-images.githubusercontent.com/60095100/192097082-8e3d0e6d-fd93-4076-9b91-19c5ea98eee5.png)

#####  5. Select Denoising parameters ?

After selecting metadatafile it will ask for enter some Denoising qiime parameters, the parameters will be selecting by demux file (How to select parameter : View next picture).

![Parameters](https://user-images.githubusercontent.com/60095100/191967572-6074fb08-6852-4c8e-8269-01b59009128b.jpg)

How to select Denoising parameters.

![denoising_parameter_graph](https://user-images.githubusercontent.com/60095100/191967642-5709a602-3671-4f45-b643-37d43e96c4b5.jpg)

#### 6. Select sampling Depth and max-depth Parameters

Entering sample_depth and max-depth for Alpha & Beta diversity and alpha rarefaction analysis, parameter would be select from table.qzv file (How to select parameter : View next picture). 

![sampling_parameters](https://user-images.githubusercontent.com/60095100/192087425-6743c2b2-9017-469e-ad63-8ed61a2a7e2f.jpg)



##### How to select sampling depth

![alpha_beta_graph_sampling_depth](https://user-images.githubusercontent.com/60095100/191969078-39df3911-c6a2-4007-8a6d-414839d7046c.jpg)

##### How to select max depth parameter

![oie_alrZ4XjugiKQ](https://user-images.githubusercontent.com/60095100/201357157-e4911eda-056b-49c8-843b-aacd9108a39a.png)

In above case maximum sampling depth is 595  




*Note : This pipeline is compatible with all qiime versions but users have to change greeen gene (gg) file according to their qiiime version 

Enjoy and send comments and complaints to khatriabhi2319@gmail.com.
