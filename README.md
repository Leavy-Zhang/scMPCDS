#Readme
This project was used for analyzing scMPCDS sequecing data and data visualization related to the paper shown as below.

## Citation

If you use the code in this repository, please cite:

Zhengyang Zhang, Liwei Zhang, Peng An, Xu Zhang, Yi Xia, Yunlu Kang, Xiaoxia Chen, Rongrong Hua, Yinhua Zhu, Yanling Hao, Yuan Huang,Yongting Luo, Junjie Luo, Guisheng Wang. Single-cell profiling of mitochondrial phenotyping-coupled mtDNA genotyping. Proc. Natl. Acad. Sci. U.S.A. (in press)

## Raw Data

The whole sequencing data is available  at https://ngdc.cncb.ac.cn/ with accession number GVM001080, GVM001458 and HRA019133.

## Set up
*Requires anaconda installation

	# 1) Create the required conda environment: 
	    conda create -n scMPCDS bioconda::bwa==0.7.17 \
			bioconda::fastp==0.23.2 \
			bioconda::samtools==1.21 \
			bioconda::ptrimmer==1.3.1 \
			bioconda::snpeff==4.5 \
			bioconda::bcftools==1.21 \
			conda-forge::python==3.12.7 \
			conda-forge::r-base==4.3.3 \
			conda-forge::r-dplyr==1.1.4 \
			conda-forge::r-ggbeeswarm==0.7.2 \
			conda-forge::r-ggplot2==3.4.4 \
			conda-forge::r-ggpubr==0.6.0 \
			conda-forge::r-ggrepel==0.9.4 \
			conda-forge::r-ggridges==0.5.4  \
			conda-forge::numpy==2.0.2 \
			conda-forge::pandas==2.2.3 \
			conda-forge::scikit-learn==1.6.0 \
			conda-forge::scipy==1.14.1 \
			conda-forge::matplotlib==3.9.3 \
			conda-forge::seaborn==0.13.2 \
			bioconda::sigprofilermatrixgenerator==1.3.6
			
	# 3) 

## Major Steps

To activate environment: 
		conda activate scMPCDS

To obtain vcf files with high quality for each cell. Run seq2vcf/mapping2vcf.sh at the directory containing raw fastq files:

	  
