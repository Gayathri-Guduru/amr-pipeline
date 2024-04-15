#!/usr/bin/env Rscript

library(tibble)
library(dplyr)
library(optparse)

# Set up command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Directory containing gene_presence_absence.csv file", metavar="path"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory for training_files. Defaults to ${params.outdir}/training_files", metavar="path")
)

# Parse command line options
args <- parse_args(OptionParser(option_list=option_list))

# Ensure the input directory is provided
if (is.null(args$input)) {
  stop("Input directory must be provided. Use --help for more information.")
}

# If output directory is not provided, assume default as per modules.config
if (is.null(args$output)) {
  args$output <- file.path(args$input, "training_files") # Assuming args$input is ${params.outdir}
}

# Construct file paths
input_file_path <- file.path(args$input, "roary", "gene_presence_absence.csv")
output_file_path <- file.path(args$output, "gene_presence_absence_encoded_matrix.csv")

# Load and clean the data
data <- read.csv(input_file_path, header = TRUE, stringsAsFactors = FALSE) %>%
  select(Gene, starts_with("SRR")) %>%
  mutate(across(-Gene, na_if, "")) %>%
  mutate(across(-Gene, ~ifelse(is.na(.x), 0, 1)))

# Write the cleaned and encoded data to a new CSV file
write.csv(data, file = output_file_path, row.names = FALSE)

# Read the newly saved file, transpose it, and tidy up
transposed_data <- read.csv("C:/Users/gguduru/Documents/office/28th_march/salmonella_enterica/gene_dataset/gene_presence_absence_encoded_matrix.csv") %>%
  column_to_rownames(var = "X") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "SRA_IDs") %>%
  # Use the first row to set column names and then remove that row
  {colnames(.) <- .[1, ]; .[-1, ]}

# Optionally, rename the first column if needed
colnames(transposed_data)[1] <- "SRA_IDs"

# Display the head of the transposed and cleaned dataframe
head(transposed_data)

write.csv(transposed_data, file = "C:/Users/gguduru/Documents/office/28th_march/salmonella_enterica/gene_dataset/gene_presence_absence_filtered_matrix.csv", row.names = FALSE
