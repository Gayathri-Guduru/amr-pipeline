# AMR(Anti-Microbial resistance - silent pandemic)

Antimicrobial resistance happens when germs like bacteria and fungi develop the ability to defeat the drugs designed to kill them. That means the germs are not killed and continue to grow. Resistant infections can be difficult, and sometimes impossible, to treat. Inorder to resolve this issue we need to know AST (Anti-microbial susceptibility test) and AFST of the microrganisms. These are types of lab tests that cultures (grows) bacteria and fungi to determine how sensitive the germ is to different antibiotics and antifungals. So, now we use Machine learning models for AST prediction. 

# Pipeline 
  ![AMR-pipeline (1)](https://github.com/Gayathri-Guduru/amr-pipeline/assets/98939664/95696b34-9b9f-4490-88ab-551957640daf)

#### Use tmux to run the pipeline.(pipeline runs even if the computer is inactive)
```tmux``` 
and use ```tmux cntrl+b d``` to detach
and ```attach -t 11```

# amr-pipeline
nextflow command to run:
```
nextflow run main.nf --database Mycobacterium_tuberculosis_test --design 's3://zymo-filesystem/home/gguduru/design_sheet.csv' -profile awsbatch --outdir 's3://zymo-filesystem/home/gguduru/results/' -work-dir 's3://zymo-filesystem/home/gguduru/tmp/' --awsqueue 'arn:aws:batch:us-east-1:002226384833:job-queue/rnaseq'
```

## create a conda env
```
conda create --name myenv \
conda activate myenv \
conda install -c bioconda samtools #install required tools
conda install -c conda-forge tmux
```

## Extract SRA id's from _Mycobacterium_tuberculosis_ species. 
Follow the steps given below:
1. Go to BV-BRC website and chose an organism (Mycobacterium_tuberculosis).
2. We need to know taxonomy id inorder to filter the dataset further - go to this website to know taxon id https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1773&lvl=3&lin=f&keep=1&srchmode=1&unlock  
3. Now back to BV-BRC website, focus on Genomes and phenotypes tabs. Click on genomes tab…a table appears(genome name, strain, genbank access etc,.)
4. Now you need to filter based on SRA accession ID’s. So there is a small + symbol on the right side. Click that and choose SRA accession under DB CROSS REFERENCE.
5. There are certains rows that are blank under SRA Accession. We need to filter them out and retain the rows that have information under SRA Accession. Clicking on SRA Accession tab arranges the ids alphabetically.
6. Now we need information of Genome ID’s(key element in genomes and phenotypes tab using for merging). Download the data in excel format to retain the trailing zeros. 

Then open Rstudio.
This is a script to generate pheno.csv file.
**i/p: BV-BRC_genome.xlsx and BV-BRC_phenotype.xlsx** 
```{r}
install.packages(c("readxl", "dplyr", "tidyr", "openxlsx"))

library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)


# Load genome file from BV-BRC
genome_data <- suppressWarnings(
  read_excel("C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/plots/BVBRC_genome.xlsx"))

# Filter any blanks or NA from SRA_Accession column in the genome file.
genome_data_filtered <- suppressWarnings(
  read_excel("C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/plots/BVBRC_genome.xlsx") %>%
    rename(
      Genome_ID = `Genome ID`, 
      SRA_Accession = `SRA Accession`
    ) %>%
    filter(SRA_Accession != "")
)

# Check if there are any blanks or NA values in the SRA_Accession column of the genome_data_filtered
any(is.na(genome_data_filtered$SRA_Accession) | genome_data_filtered$SRA_Accession == "") 

# Load phenotype data and filter out rows with NA in Resistant_Phenotype into phenotype_data_filtered
phenotype_data <- suppressWarnings(
  read_excel("C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/BVBRC_phenotype.xlsx") %>%
    rename(
      Genome_ID = `Genome ID`, 
      Resistant_Phenotype = `Resistant Phenotype`
    )
)

# Create phenotype_data_filtered with rows where Resistant_Phenotype is not NA
phenotype_data_filtered <- phenotype_data %>%
  filter(!is.na(Resistant_Phenotype))

# Merge datasets based on the 'Genome_ID' column
# Only genome_data with non-blank SRA_Accession is included
merged_dataset <- merge(genome_data_filtered, phenotype_data_filtered, by = "Genome_ID", all = FALSE)
merged_dataset <- merged_dataset %>%
  select(SRA_Accession, Genome_ID, Resistant_Phenotype, everything())

# Check if there are any blanks or NA values in the SRA_Accession column of the merged dataset
any(is.na(merged_dataset$SRA_Accession) | merged_dataset$SRA_Accession == "")

# Display the merged dataset
write.xlsx(merged_dataset, "C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/merged_dataset.xlsx", rownames = FALSE)

# Read the data from the Excel file
data <- suppressWarnings(read_excel("C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/merged_dataset.xlsx"))

# Rename the column from 'Resistant.Phenotype' to 'Phenotype'
# And replace "IS" values with NA
data <- data %>%
  rename(Phenotype = Resistant_Phenotype) %>%
  mutate(Phenotype = case_when(
    Phenotype == "Resistant" ~ "R",
    Phenotype == "Susceptible" ~ "S",
    Phenotype == "Intermediate" ~ "R",
    Phenotype == "IS" ~ NA_character_,
    TRUE ~ as.character(Phenotype)
  ))

# Resolve duplicates by grouping and summarizing, take the first occurrence
data_aggregated <- data %>%
  group_by(Genome_ID, SRA_Accession, Antibiotic) %>%
  summarise(Phenotype = first(Phenotype), .groups = 'drop') %>%
  distinct(Genome_ID, SRA_Accession, Antibiotic, .keep_all = TRUE)

# Spread the data to wide format, retaining SRA_Accession
data_T <- data_aggregated %>%
  spread(key = Antibiotic, value = Phenotype, fill = NA_character_)

# Ensure each antibiotic column has at least 5 "R" and 5 "S" values, and remove columns with all NAs
valid_columns <- colSums(data_T == "R", na.rm = TRUE) >= 5 & colSums(data_T == "S", na.rm = TRUE) >= 5
data_T <- data_T %>%
  select(SRA_Accession, Genome_ID, names(valid_columns)[valid_columns])

# Remove columns that have only NA values
data_T <- data_T[, colSums(!is.na(data_T)) > 0]

# Remove duplicate SRA IDs from the entire dataframe, keeping the first occurrence
data_T <- data_T %>% distinct(SRA_Accession, .keep_all = TRUE)

# Write the transformed data with unique SRA IDs to a CSV file
write.csv(data_T, "C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/pheno.csv", row.names = FALSE)

# Additionally, save the pheno data to an Excel file
write.xlsx(data_T, "C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/pheno.xlsx", rownames = FALSE)

# Extract SRA_Accession column (already unique)
sra_ids <- data_T$SRA_Accession

# Write the unique SRA IDs to a text file, one per line
writeLines(sra_ids, "C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/sra_ids.txt")

# Additional script to format the SRA IDs in the output file
# Read the input file
input_file <- "C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/sra_ids.txt"
data <- readLines(input_file)

sra_count <- length(data)
cat("The number of unique SRA IDs in the 'sra_ids.txt' file is:", sra_count, "\n")

# Split the IDs by comma and combine them into a single string with new lines
formatted_data <- unlist(strsplit(data, ","))

# Write the output to a new file
output_file <- "C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/sra_ids_formatted_r.txt"
writeLines(formatted_data, output_file)

# Read the formatted data from the file
formatted_data <- readLines("C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/plots/plots/sra_ids_formatted_r.txt")

# Count the number of SRA IDs
sra_count <- length(formatted_data)

# Print the count
cat("The number of SRA IDs in the 'sra_ids_formatted_r.txt' file is:", sra_count, "\n")
```
**o/p: merged_dataset.xlsx, pheno.csv, sra_ids.txt**

## After creating you own dataset
1. Total no of samples are 10804. Now from those sra IDs retrieve the fastq files(script attached below).
2. Retrieve the refseq reference genome for the species 
3. Upload the reference genomes and fastq files to your aws bucket
4. Amend the igenomes.config with the reference details
5. Create the design sheet (script attached below)
6. Run the pipeline!

## 1. Retrieve fastq files from sra using the below python script.
```
import os
import subprocess
import boto3

# Initialize S3 client
s3_client = boto3.client('s3')

# Function to download and compress FASTQ files from SRA
def download_and_compress_fastq(sra_id):
    try:
        # Download the SRA file using prefetch
        command_prefetch = f"prefetch {sra_id}"
        subprocess.run(command_prefetch, shell=True, check=True)

        # Convert SRA file to FASTQ using fasterq-dump (without gzip option)
        command_fasterq = f"fasterq-dump --split-files --outdir ./ --skip-technical {sra_id}"
        subprocess.run(command_fasterq, shell=True, check=True)

        # Compress the FASTQ files using gzip
        for file in os.listdir():
            if file.startswith(sra_id) and file.endswith(".fastq"):
                command_gzip = f"gzip {file}"
                subprocess.run(command_gzip, shell=True, check=True)

    except subprocess.CalledProcessError as e:
        print(f"Error processing {sra_id}: {e}")

# Function to upload to S3
def upload_to_s3(file_name, bucket_name, s3_path):
    try:
        s3_client.upload_file(file_name, bucket_name, s3_path)
    except Exception as e:
        print(f"Failed to upload {file_name} to S3: {e}")

# Function to process a single SRA ID
def process_sra_id(sra_id, bucket_name, s3_base_path):
    print(f"Processing {sra_id}...")

    # Download and compress FASTQ files
    download_and_compress_fastq(sra_id)

    # Upload to S3 and delete local files
    for file_name in os.listdir():
        if file_name.startswith(sra_id) and file_name.endswith(".fastq.gz"):
            s3_path = os.path.join(s3_base_path, file_name)
            upload_to_s3(file_name, bucket_name, s3_path)
            os.remove(file_name)  # Delete the local file after upload

    print(f"Completed processing for {sra_id}")

# Main function to iterate over SRA IDs sequentially
def process_sra_ids(sra_ids, bucket_name, s3_base_path):
    for sra_id in sra_ids:
        process_sra_id(sra_id, bucket_name, s3_base_path)

if __name__ == "__main__":
    # Path to SRA IDs file
    sra_ids_file = "/home/gguduru/storage/subset/Mycobacterium_tuberculosis_1773/sra_ids.txt"

    # Read the SRA IDs from the file
    with open(sra_ids_file, 'r') as file:
        sra_ids = [line.strip() for line in file.readlines()]

    # S3 bucket name and path
    bucket_name = "zymo-filesystem"
    s3_base_path = "home/gguduru/Mycobacterium_tuberculosis_1773/fastq_files/"

    # Process the SRA IDs one by one
    process_sra_ids(sra_ids, bucket_name, s3_base_path)
```
**o/p: The fastq.gz files are uploaded to s3 bucket.**

## 2. Script for design_sheet.csv

```{r}
# Load necessary library
library(dplyr)

# Read the SRA ID list
sra_ids <- readLines("C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/sra_ids.txt")

# Define the S3 base path
s3_base_path <- "s3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/fastq_files/"

# Define a function to generate S3 file paths based on SRA ID
generate_s3_paths <- function(sra_id) {
  data.frame(
    sample = sra_id,
    read_1 = paste0(s3_base_path, sra_id, "_1.fastq.gz"),
    read_2 = paste0(s3_base_path, sra_id, "_2.fastq.gz")
  )
}

# Create a data frame with S3 paths for each SRA ID
s3_files_df <- do.call(rbind, lapply(sra_ids, generate_s3_paths))

# Write the design sheet to a CSV file
write.csv(s3_files_df, "C:/Users/gguduru/OneDrive - Zymo Research/Local/amr_pipeline/Mycobacterium_tuberculosis/design_sheet.csv", row.names = FALSE, quote = FALSE)

# Output the data frame to the console (optional)
print(s3_files_df)
```
- Create fastq_files folder on the s3 bucket for ex:(s3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/fastq_files/)
- Create reference_genome folder in the same path on s3 bucket (s3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/reference_genome/)
- Create results folder (s3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/results/)
- Place the ```pheno.csv``` in the same path.
- Place the ```design_sheet.csv``` in the same path.

## 3. Next get a reference genome 
Go to NCBI -> select taxonomy and type species name(salmonella enterica)

![image](https://github.com/Gayathri-Guduru/amr-pipeline/assets/98939664/a752d352-bbee-4772-b6ab-cce121a4644a)
 
In this page, on the right side table click on genome. You are redirected to genome page. Click on that..(u can download to your local pc or use curl command to download on VM)
I directly downloaded the fasta file and transferred to vm using Winscp and then to s3.

**o/p: reference genome fasta file is generated.**

## 4. Upload the reference genomes and fastq files to your aws bucket.
First, create required folders on s3. (Here, ```s3://zymo-filesystem/home/gguduru/``` is my s3 bucket where all files and folders are stored)
Now, I created `fastq_files` folder to place my input fastq files and `reference_genome` folder to place my reference fasta file along with the indexed files.

To get indexed files

```bwa index "fasta file"```

**o/p Indexed files**

Now, send these files to s3.
```
## Transferring files to s3 from vm
aws s3 cp /home/gguduru/ s3://zymo-filesystem/tmp/gguduru/fastq_files/ --recursive --exclude "*" --include "*.fastq.gz" # to transfer fastq files
aws s3 cp /home/gguduru/ s3://zymo-filesystem/tmp/gguduru/reference_genome/ --recursive --exclude "*" --include "*GCF*" # to transfer indexed and fasta files
```

## 4. Amend the igenomes.config with the reference details

on `merge` branch of `amr-pipeline` go to `conf -> igenomes.config`
Now, change the path of the index files to the file path on s3.
```
params{
    databases {
        'Mycobacterium_tuberculosis_test' {
           index_path          = "s3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/reference_genome/*"
           index_name          = "./Mycobacterium_tuberculosis_GCF_000195955.2_ASM19595v2_genomic.fna"
           fasta_path          = "s3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/reference_genome/Mycobacterium_tuberculosis_GCF_000195955.2_ASM19595v2_genomic.fna"
	   gff_path            = "s3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/reference_genome/Mycobacterium_tuberculosis_genomic.gff"
           pheno_path          = "s3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/pheno.csv"
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
nextflow run main.nf --database Mycobacterium_tuberculosis_test --design 's3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/design_sheet.csv' -profile awsbatch --outdir 's3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/results/' -work-dir 's3://zymo-filesystem/home/gguduru/Mycobacterium_tuberculosis_1773/tmp/' --awsqueue 'arn:aws:batch:us-east-1:002226384833:job-queue/rnaseq'
```

## 7. Developed scripts using Nextflow.
1. Quality Control using Fastq.
2. Genome Assembly using SPAdes.
3. Quality control of assembled contigs using Quast.
4. Genome Annotation using Prokka.
5. Pangenome assesment using Roary/genAPI.
6. Build the input train datasets (3 input files: pheno.csv, gene_presence_absence.csv, snp_input_N.transposed.csv.gz)
7. Prepare the train data
8. Train the models using combine_features_and_train_models.

## Merge both the pipelines.
1. Created a new subworkflow ```input_train_data.nf``` for ```merge_snps.nf``` and ```gene_filtering.nf``` as the outputs of these along with snp_input_N.transposed.csv.gz, pheno.csv file are provided as input for the next process for ```prep_training.nf```
2. Now, the combine_features_train_models.nf subworkflow runs prep_train_data and train_models scripts.


## Testing the pipelines on multiple datasets before merging into main branch
