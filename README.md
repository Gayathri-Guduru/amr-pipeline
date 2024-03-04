# amr-pipeline
nextflow command to run:
```
 nextflow run main.nf --database salmonella_test  --design 's3://zymo-filesystem/home/gguduru/latest_design_sheet_2.csv' -profile awsbatch --outdir 's3://zymo-filesystem/home/gguduru/results/' -work-dir 's3://zymo-filesystem/home/gguduru/tmp/' --awsqueue 'arn:aws:batch:us-east-1:002226384833:job-queue/rnaseq'
```

# create a conda env

```
conda create --name myenv \
conda activate myenv \
conda install -c bioconda samtools #install required tools
 conda install -c conda-forge tmux
```

# bwa index
```bwa index "fasta file"```

```for f in *1.fastq.gz; do
    # Get the base name without the _1.fastq.gz suffix
    base=$(basename "$f" _1.fastq.gz)

    # Run bwa mem for the pair of files
    bwa mem -t 10 GCF_000006945.2_ASM694v2_genomic.fna "$base"_1.fastq.gz "$base"_2.fastq.gz > "$base".sam
done
```
# to check the alignment rate to the reference files after generating index files using bwa 
```samtools flagstat ../SRR2566949.sam```

# transferring files to anf from vm to scp
```
aws s3 cp s3://your-bucket-name/path/in/s3/ /path/on/local/machine/
aws s3 cp /path/on/local/machine/ s3://your-bucket-name/path/in/s3/ --recursive
```
