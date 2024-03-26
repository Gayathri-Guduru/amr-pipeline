# AMR(Anti-Microbial resistance - silent pandemic)

Antimicrobial resistance happens when germs like bacteria and fungi develop the ability to defeat the drugs designed to kill them. That means the germs are not killed and continue to grow. Resistant infections can be difficult, and sometimes impossible, to treat. Inorder to resolve this issue we need to know AST (Anti-microbial susceptibility test) and AFST of the microrganisms. These are types of lab tests that cultures (grows) bacteria and fungi to determine how sensitive the germ is to different antibiotics and antifungals. So, now we use Machine learning models for AST prediction. 

# Pipeline 
  ![AMR-pipeline (1)](https://github.com/Gayathri-Guduru/amr-pipeline/assets/98939664/95696b34-9b9f-4490-88ab-551957640daf)

# amr-pipeline
nextflow command to run:
```
nextflow run main.nf --database salmonella_test --design 's3://zymo-filesystem/home/gguduru/design_sheet.csv' -profile awsbatch --outdir 's3://zymo-filesystem/home/gguduru/results/' -work-dir 's3://zymo-filesystem/home/gguduru/tmp/' --awsqueue 'arn:aws:batch:us-east-1:002226384833:job-queue/rnaseq'
```

## create a conda env

```
conda create --name myenv \
conda activate myenv \
conda install -c bioconda samtools #install required tools
conda install -c conda-forge tmux
```

## Extract SRA id's from _salmonella enterica_ species. 
Follow the steps given below:
1. Go to BV-BRC website and chose an organism (salmonella). 
2. Now focus on Genomes and phenotypes tabs. Click on genomes tab…a table appears(genome name, strain, genbank access etc,.)
3. Now you need to filter based on SRA accession ID’s. So there is a small + symbol on the right side. Click that and choose SRA accession under DB CROSS REFERENCE.
4. There are certains rows that are blank under SRA Accession. We need to filter them out and retain the rows that have information under SRA Accession. Clicking on SRA Accession tab arranges the ids alphabetically.
5. Now we need information of Genome ID’s(key element in genomes and phenotypes tab using for merging). Download the data in CSV format. 

Then open Rstudio.
```{r}
# Install the tidyverse package
install.packages("tidyverse")

# Load the tidyverse package for data manipulation
library(tidyverse)

# URL of the CSV file
csv_data <- "C:/Users/gguduru/Downloads/office/BVBRC_genome.csv"

# Read the CSV file into a data frame
genome_data <- read.csv(csv_data, header = TRUE, colClasses = "character")

# Convert the 'Genome.ID' column to string
#data$Genome.ID <- as.character(data$Genome.ID)

# Column to shift to the first position
#SRA_Accession <- data$SRA.Accession

# Remove the column from its original position
#data <- data[, !(names(data) %in% "SRA.Accession")]

# Insert the column at the first position
#data <- cbind(SRA_Accession, data)

# Filter out rows with blanks in the SRA_Accession column
#genome_data <- data %>% filter(!is.na(SRA_Accession) & SRA_Accession != "")

# Print the filtered data or save it to a new CSV file
# print(genome_data)
#write.csv(genome_data, "Genome_data.csv", row.names = FALSE)

# Read the CSV file into a data frame
csv_data_2 <- "C:/Users/gguduru/Downloads/office/BVBRC_phenotype_amr.csv"

# Read the CSV file into a data frame
phenotype_data <- read.csv(csv_data_2, header = TRUE, colClasses = "character")

Merging the datasets based on Genome Id's

# Merge datasets based on the 'ID' column
merged_dataset <- merge(genome_data, phenotype_data, by = "Genome.ID", all = FALSE)

# Display the merged dataset
write.csv(merged_dataset, "merged_dataset.csv", row.names = FALSE)
```

## Create a list of unique sra ids from the merged dataframe
```
# Create a list of unique SRA IDs
unique_sra_ids <- unique(merged_dataset$SRA.Accession)
sorted_sra_ids <- unique_sra_ids[order(unique_sra_ids, decreasing = TRUE)]

# Select the top 10 SRA IDs
top_10_sra_ids <- head(sorted_sra_ids, 10)

#print("Top 10 SRA IDs:")
print(top_10_sra_ids)

# Write the top 10 SRA IDs to a text file
writeLines(top_10_sra_ids, "top_sra_ids.txt")
```
o/p: We have a list of unique sra ids from the merged dataframe.

## After creating you own dataset
1. Selected ~10 sra IDs and retrieve the fastq files from the sra (I.e. with curl/wget command) to your VM
2. Retrieve the refseq reference genome for the species 
3. Upload the reference genomes and fastq files to your aws bucket
4. Amend the igenomes.config with the reference details
5. Create the design sheet
6. Run the pipeline!

## 1. Retrieve fastq files from sra
```
#!/bin/bash
# Array of 10 SRA IDs
sra_ids=("SRR3665082" "SRR2567031" "SRR2567121" "SRR3664650" "SRR2566949" "SRR2567089" "SRR2567095" "SRR2567181" "SRR2567114" "SRR3665078")

# Loop through the array and perform prefetch for each SRA ID
for sra_id in "${sra_ids[@]}"
do
  prefetch $sra_id

    # Run fastq-dump for single-end data
    fastq-dump --split-3 --gzip $sra_id
done
```
**o/p: This provides fastq files in zipped format.**

## 2. Next get a reference genome 
Go to NCBI -> select taxonomy and type species name(salmonella enterica)

![image](https://github.com/Gayathri-Guduru/amr-pipeline/assets/98939664/a752d352-bbee-4772-b6ab-cce121a4644a)
 
In this page, on the right side table click on genome. You are redirected to genome page. Click on that..(u can download to your local pc or use curl command to download on VM)
I directly downloaded the fasta file and transferred to vm using Winscp and then to s3.

**o/p: reference genome fasta file is generated.**

## 3. Upload the reference genomes and fastq files to your aws bucket.
First, create required folders on s3. (Here, ```s3://zymo-filesystem/home/gguduru/``` is my s3 bucket where all files and folders are stored)
Now, I created `fastq_files` folder to place my input fastq files and `reference_genome` folder to place my reference fasta file along with the indexed files.

To get indexed files
```bwa index "fasta file"```

```for f in *1.fastq.gz; do
    # Get the base name without the _1.fastq.gz suffix
    base=$(basename "$f" _1.fastq.gz)

    # Run bwa mem for the pair of files
    bwa mem -t 10 GCF_000006945.2_ASM694v2_genomic.fna "$base"_1.fastq.gz "$base"_2.fastq.gz > "$base".sam
done
```
**o/p Indexed files along with .sam files**

Now, send these files to s3.
```
## Transferring files to s3 from vm
aws s3 cp /home/gguduru/ s3://zymo-filesystem/tmp/gguduru/fastq_files/ --recursive --exclude "*" --include "*.fastq.gz" # to transfer fastq files
aws s3 cp /home/gguduru/ s3://zymo-filesystem/tmp/gguduru/reference_genome/ --recursive --exclude "*" --include "*GCF*" # to transfer indexed and fasta files
```

Cross-check check the alignment rate to the reference files after generating index files using bwa 
```samtools flagstat ../SRR2566949.sam```

## 4. Amend the igenomes.config with the reference details

on `gguduru` branch of `amr-pipeline` go to `conf -> igenomes.config`
Now, change the path of the index files to the file path on s3.
```
params{
    databases {
        'salmonella_test' {
           index_path = "s3://zymo-filesystem/home/gguduru/reference_genome/*"
           index_name = "./GCF_000006945.2_ASM694v2_genomic.fna"
           fasta_path = "s3://zymo-filesystem/home/gguduru/reference_genome/GCF_000006945.2_ASM694v2_genomic.fna"
      }
    }
}
```
## 5. Create the design sheet.
This is the sheet located in `gguduru` branch of `amr-pipeline`. Go to `test_data -> design_sheet.csv`
Now, replace the sample, read_1, read_2 with the paths on s3.

![image](https://github.com/Gayathri-Guduru/amr-pipeline/assets/98939664/028fa7be-112a-420f-bf93-36737520df1c)

## 6. Run the pipeline!
```
nextflow run main.nf --database salmonella_test --design 's3://zymo-filesystem/home/gguduru/design_sheet.csv' -profile awsbatch --outdir 's3://zymo-filesystem/home/gguduru/results/' -work-dir 's3://zymo-filesystem/home/gguduru/tmp/' --awsqueue 'arn:aws:batch:us-east-1:002226384833:job-queue/rnaseq'
```

## 7. Start developing scripts using Nextflow.
Quality Control using Fastq.
Genome Assembly using SPAdes.
Quality control of assembled contigs using Quast.
Genome Annotation using Prokka.
Pangenome assesment using Roary/genAPI.
