#! /bin/bash

########## Individual Component RNAseq Pipeline Script ##########
## Group: Bone
## Member: Emma Hruska
## Adapted from: Schwartz, T. (2025, January 9). Steps0-3_RNAseq_StandardPipeline_Clean-T_Map-H_Count-S.sh [Shell script]. 
## GitHub repository, https://github.com/Schwartz-Lab-at-Auburn/FunGen2025/tree/main. Accessed March 30, 2025.
## Purpose: The purpose of this script is to run the full RNAseq pipeline.
## Paramenters below were used to submit this script:
## 		queue: bigmem
##		core: 12
##		time limit (HH:MM:SS): 18:00:00 
##		memory: 240gb	
###############################################

########## Load modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1
module load multiqc
module load trimmomatic/0.39
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread

## Set the stack size to unlimited
ulimit -s unlimited

## Turn echo on so all commands are echoed in the output log
set -x

########## Define variables and make directories

## Make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubclsd0318          

WD=/scratch/$MyID/Individual_Component            
DD=$WD/RawData
RDQ=RawDataQuality
adapters=AdaptersToTrim_All.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing
CD=$WD/CleanData            			
PCQ=PostCleanQuality
REFD=$WD/BoxerRefGenome         ## This directory contains the indexed reference genome
MAPD=$WD/Map_HiSat2           		
COUNTSD=/$WD/Counts_StringTie      
RESULTSD=/home/$MyID/Individual_Component/Counts_2025     
REF=GCF_000002285.5_Dog10K_Boxer_Tasha_genomic   ## This is what the "easy name" will be for the genome reference

## Make directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there
mkdir -p ${WD}
mkdir -p ${DD}

## Move to the Data Directory
cd ${DD}

########## Download data files from NCBI
## From SRA use the SRA Toolkit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	## This downloads the SRA file and converts it to .fastq
		## Splits the files into R1 and R2 (forward reads, reverse reads)

fasterq-dump SRR22459519
fasterq-dump SRR22459520
fasterq-dump SRR22459526
fasterq-dump SRR22458622
fasterq-dump SRR22458624
fasterq-dump SRR22458632

########## FastQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results, a zipped file of results, and an .html file for each sample
mkdir ${WD}/${RDQ}
fastqc *.fastq --outdir=${WD}/${RDQ}

########## MultiQC to summarize the FastQC results!
cd ${WD}/${RDQ}
multiqc ${WD}/${RDQ}

## Tarball the directory containing the FastQC and MultiQC results to easily bring it back to your computer to evaluate
tar cvzf ${RDQ}.tar.gz  ${WD}/${RDQ}/*

########## Cleaning the data with Trimmomatic

## Make directories to hold the cleaned data files and a directory to hold the results for assessing quality of the cleaned data
mkdir ${CD}
mkdir ${WD}/${PCQ}

## Move to Raw Data Directory
cd ${DD}

## Make list of file names to run through Trimmomatic
        ## This line is a set of piped (|) commands
        ## ls means make a list, 
        ## grep means grab all file names that end in ".fastq", 
        ## cut that name into elements at each "_" and keep the first element (-f 1),
        ## sort the list, keep only the unique names, and put it into a file named "list"
ls | grep ".fastq" |cut -d "_" -f 1 | sort | uniq > list

## Copy over the list of sequencing adapters that we want Trimmomatic to look for (along with its default adapters)
cp /home/${MyID}/class_shared/AdaptersToTrim_All.fa . 

## Run a while loop to process through the file names in the list and process them with the Trimmomatic code
while read i
do

        ## Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format
	       java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar   \
					PE -threads 6 -phred33 \
        	"$i"_1.fastq "$i"_2.fastq  \
       	 ${CD}/"$i"_1_paired.fastq ${CD}/"$i"_1_unpaired.fastq  ${CD}/"$i"_2_paired.fastq ${CD}/"$i"_2_unpaired.fastq \
       	 ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
        
                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter
                ## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
                ## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across  
                ## requiredQuality: specifies the average quality required.

########## FastQC to assess quality of the cleaned sequence data
	## FastQC: run on each of the data files that have 'All' to check the quality of the data
	## The output from this analysis is a folder of results and a zipped file of results

fastqc ${CD}/"$i"_1_paired.fastq --outdir=${WD}/${PCQ}
fastqc ${CD}/"$i"_2_paired.fastq --outdir=${WD}/${PCQ}

done<list			## This is the end of the loop

########## Run MultiQC to summarize the FastQC results

## Move to the directory with the cleaned data
cd ${WD}/${PCQ}
multiqc ${WD}/${PCQ}

##########  Now compress your results files from FastQC 

## Tarball the directory containing the FastQC results so we can easily bring it back to our computer to evaluate
tar cvzf ${PCQ}.tar.gz ${WD}/${PCQ}/*

########## Mapping and counting

## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

## Prepare the Reference Index for mapping with HISAT2

cd $REFD

## Copy the reference genome (.fasta) and the annotation file (.gff3) to this REFD directory
cp /home/${MyID}/ncbi_dataset/data/GCF_000002285.5/${REF}.fna .
cp /home/${MyID}/ncbi_dataset/data/GCF_000002285.5/genomic.gff .

## Identify exons and splice sites on the reference genome
gffread genomic.gff -T -o ${REF}.gtf    ## gffread converts the annotation file from .gff3 to .gft format for HISAT2 to use
hisat2_extract_splice_sites.py ${REF}.gtf > ${REF}.ss
hisat2_extract_exons.py ${REF}.gtf > ${REF}.exon

## Create a HISAT2 index for the reference genome
hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fna Boxer_index

########## Map and count the data using HISAT2 and StringTie

## Move to the Data Directory
cd ${CD}  ## This is where the clean paired reads are located

## Create list of .fastq files to map
## Grab all .fastq files, cut on the underscore, use only the first of the cuts, sort, and use uniq to put in list
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list Example: SRR629651

## Move to the directory for mapping
cd ${MAPD}

## Move the list of unique IDs from the original files to map
mv ${CD}/list . 

## Process the samples in the list one by one using a while loop
while read i;
do
  ## HISAT2 is the mapping program
  ## -p indicates number of processors, --dta reports alignments for StringTie, --rf is the read orientation
   hisat2 -p 6 --dta --phred33       \
    -x "${REFD}"/Boxer_index       \
    -1 "${CD}"/"$i"_1_paired.fastq  -2 "${CD}"/"$i"_2_paired.fastq      \
    -S "$i".sam

    ## view: convert the SAM file into a BAM file  
    ## -bS: BAM is the binary format corresponding to the SAM text format
    ## sort: convert the BAM file to a sorted BAM file
  
  samtools view -@ 6 -bS "$i".sam > "$i".bam  

    ## This is sorting the BAM file using 6 threads and produces a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 6  "$i".bam  -o  "$i"_sorted.bam

    ## Index the BAM, get mapping statistics, and put them in a text file to look at
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt
  
########## Counting

## StringTie is the program that counts the reads that are mapped to each gene, exon, and transcript model
  ## The output from StringTie is counts folders in a directory that is ready to bring into the R program Ballgown 

	mkdir "${COUNTSD}"/"$i"
	stringtie -p 6 -e -B -G  "${REFD}"/"${REF}".gtf -o "${COUNTSD}"/"$i"/"$i".gtf -l "$i"   "${MAPD}"/"$i"_sorted.bam

done<list

########## Copy results to Home Directory

## These are your stats files from Samtools
cp *.txt ${RESULTSD}

## Move to the Counts Directory
cd ${COUNTSD}

## Run the python script prepDE.phy to prepare your data for downstream analysis (converts the files in your Ballgown folder to a count matrix)
cp /home/${MyID}/class_shared/prepDE.py3 .

 prepDE.py3 -i ${COUNTSD}

## Copy the final results files (the count matrices that are .csv) to your home directory
cp *.csv ${RESULTSD}

## Move these results files to your personal computer for downstream statistical analyses in RStudio
