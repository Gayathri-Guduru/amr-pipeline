# amr-pipeline
amr-pipeline

#create a conda env

```
conda create --name myenv \
conda activate myenv \
conda install -c bioconda samtools #install required tools
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
