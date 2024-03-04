# amr-pipeline
nextflow command to run:
```
 nextflow run main.nf --database salmonella_test  --design 's3://zymo-filesystem/home/gguduru/latest_design_sheet_2.csv' -profile awsbatch --outdir 's3://zymo-filesystem/home/gguduru/results/' -work-dir 's3://zymo-filesystem/home/gguduru/tmp/' --awsqueue 'arn:aws:batch:us-east-1:002226384833:job-queue/rnaseq'
```

## create a conda env

```
conda create --name myenv \
conda activate myenv \
conda install -c bioconda samtools #install required tools
conda install -c conda-forge tmux
```

## Extract SRA id's from _salmonella enterica_ species. Follw the steps given below:
Go to BV-BRC website and chose an organism (salmonella). 
Now focus on Genomes and phenotypes tabs. Click on genomes tab…a table appears(genome name, strain, genbank access etc,.)
Now you need to filter based on SRA accession ID’s. So there is a small + symbol on the right side. Click that and choose SRA accession under DB CROSS REFERENCE.
There are certains rows that are blank under SRA Accession. We need to filter them out and retain the rows that have information under SRA Accession. Clicking on SRA Accession tab arranges the ids alphabetically.
Now we need information of Genome ID’s(key element in genomes and phenotypes tab using for merging). Download the data in CSV format. 

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
**TESTING**

```{r}
# Subset columns from both datasets
subset1 <- subset(genome_data, select = c("Genome.ID"))
subset2 <- subset(phenotype_data, select = c("Genome.ID"))

# Merge based on common IDs
merged_dataset_2 <- merge(subset1, subset2, by = "Genome.ID", all = FALSE)
```

```{r}
# Extract common genome IDs
common_genome_ids <- intersect(genome_data$Genome.ID, phenotype_data$Genome.ID)

# Filter datasets based on common genome IDs
filtered_dataset1 <- genome_data[genome_data$Genome.ID %in% common_genome_ids, ]
filtered_dataset2 <- phenotype_data[phenotype_data$Genome.ID %in% common_genome_ids, ]

# Merge filtered datasets based on common genome IDs
merged_filtered_datasets <- merge(filtered_dataset1, filtered_dataset2, by = "Genome.ID", all = FALSE)

# Check if the values are the same for common genome IDs
are_values_same <- all.equal(filtered_dataset1, filtered_dataset2)
```

```{r}
# Number of rows in the dataframe
num_rows <- 34971

# Create a dataframe with genome IDs containing trailing zeros
test_dataframe <- data.frame(
  Genome_ID = sprintf("ID%03d", seq_len(34971)),
  Value1 = rnorm(34971),
  Value2 = runif(34971)
)

# Display the test dataframe
print(test_dataframe)
```

```{r}

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
```
ing AMR.rmd…]()


## Bwa index
```bwa index "fasta file"```

```for f in *1.fastq.gz; do
    # Get the base name without the _1.fastq.gz suffix
    base=$(basename "$f" _1.fastq.gz)

    # Run bwa mem for the pair of files
    bwa mem -t 10 GCF_000006945.2_ASM694v2_genomic.fna "$base"_1.fastq.gz "$base"_2.fastq.gz > "$base".sam
done
```
## Cross-check check the alignment rate to the reference files after generating index files using bwa 
```samtools flagstat ../SRR2566949.sam```

## Transferring files to anf from vm to scp
```
aws s3 cp s3://your-bucket-name/path/in/s3/ /path/on/local/machine/
aws s3 cp /path/on/local/machine/ s3://your-bucket-name/path/in/s3/ --recursive
```
