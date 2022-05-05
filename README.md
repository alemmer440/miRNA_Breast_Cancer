# Investigating the importance and mechanisms of aberrant miRNA expression and regulation in human breast cancer 

## Repository Overview

miRNA are short, non-coding RNAs that regulate gene expression at the post-transcriptional level by binding via seqeunce complementarity to mRNA targets. Usually this causes mRNA degradation, but it can also promote mRNA translation. miRNAs are known to play a role in tumorigenesis. Two possible mechanisms for this are the upregulation of miRNA that degrade tumor suppressor mRNAs, or the downregulation of miRNA the degrade oncogene mRNAs. The goal of the project is to study the role of miRNA in breast cancer, and answer the following questions: 1) Which miRNA are differentially expressed in breast cancer? 2) Can miRNA data alone indicate breast cancer? 3) What are the relationships between miRNA and mRNA in the data? What are the functional impacts of these relationships and miRNA aberrant expression? 

## Data

All data is pulled from the National Cancer Institute GDC Data Portal, which can be found at this link: https://portal.gdc.cancer.gov/. Specifically, we filtered for female breast cancer patients from The Cancer Genome Atlas BRCA project that had both mRNA and microRNA data files. We filtered separately for cancerous and normal tissue. For mRNA files, we selected the versions using FPKM normalization, and for microRNAs we selected the version without isoforms. When downloaded, the data came in a highly nested directory, so we flattened all of the files into one directory. We also downloaded the associated json file containing the metadata, from which we can determine which files belong to each patient based on the file Case IDs. In total, there are cancerous data files for 1,100 patients, and normal data files for 268 of these patients. 

Note: Due to git LFS storage limits, you may not be able to download the data directly from this repository. Please go to the link listed above to filter and download the data. The nested directory can be flattened using the following command:
find /dir1 -mindepth 2 -type f -exec mv -i '{}' /dir1 ';'
where /dir1 is the path to the root directory. 

## Folder Structure

Data: Contains two sub-directories, one with all of the data and one with just the normal sample data. Also contains the metadata file, which has information about which files come from the same patients. 

Figures: Contains all the figures produced in analysis.

Organized_Data: The original data was reorganized into many different formats for the different analyses. This folder contains all of these csv files. Note - due to git LFS storage limits, you may not be able to these directly from this repository. All of them can be reproduced using our code.

Scripts: Contains all the code used for analysis. 

## Installation

All python files and notebooks run using python3. The .R file runs using R version 4.1.3.

For python code, the following packages are required for the code to run, and were all installed using conda:
* Jupyter Notebook
* Numpy
* Json
* Seaborn
* Matplotlib
* Pandas
* Scipy
* Sci-kit learn

A list of package versions for the environment the code was created in, and the command you can use to create the same conda enviroment, in is in the file requirements.txt. Alternatively, these packages can be installed using pip. 

For R code, the following packages are required: 
* Bioconductor (DESeq2) v3.15
* ggplot2 version 3.3.6

Once the necessary packages are installed, run the following command:  
git clone https://github.com/alemmer440/HW4.git

If the Data does not download, follow the instructions above to download the data from the source, and move the files into the repository directory.

Finally, navigate to Scripts folder. The scripts should be run in the following order:
* Organize_Classify_Data.ipynb: Uses the metadata to create a dictionary that has patients' IDs as the keys, and lists of all the file names associated with the patient as the values. Also creates a list of IDs of patients that only have cancerous samples. This information gets output into csv files.
* Create_Cancer_Data_Csvs.ipynb: Generates 2 matricies, for miRNA and mRNA respectively, combing the data for each cancerous sample in the dataset. Each row will be a patient, and each column will be a miRNA/mRNA. Matricies are output into csv files.
* Generate_data_for_DESeq.py: Creates csv files from the data in the correct format to use for DEseq for differential expression analysis.
* DiffEx_analysis.R: Towards question 1 - which miRNA are differentially expressed in breast cancer - runs differential expression analysis for miRNA between normal and cancerous tissue. Generates a volcano plot.
* PCA_and_RF.py: Towards question 2 - can miRNA data alone indicate breast cancer - runs PCA and random forest classifier using miRNA data from cancerous and healthy tissue. Generates PCA plot, confusion matrix, and ROC curve.
* Correlation_Matricies.ipynb: Towards question 3 - what are the relationships between miRNA and mRNA in the data, and what are the functional impacts of these relationships and miRNA aberrant expression - runs correlation analysis between miRNA and mRNA data. Generates correlation heatmaps.
