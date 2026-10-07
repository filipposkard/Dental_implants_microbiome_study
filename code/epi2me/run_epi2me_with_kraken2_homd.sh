#kraken2
#ehomd database
#--fastq $PROJECT_DIR/data/raw/both_runs  \
#kraken2_confidence 0.2 \#https://link.springer.com/article/10.1007/s42994-024-00178-0
PROJECT_DIR="/mnt/raid1b/philip/dental_implants"
cd $PROJECT_DIR/data/epi2me
nextflow run epi2me-labs/wf-16s \
    -profile singularity \
    --fastq /mnt/raid1b/philip/dental_implants/data/raw/both_runs  \
    --min_len 600 \
    --classifier kraken2 \
    --database /mnt/raid1b/philip/databases/HOMD_16_03_kraken2/eHOMD_kraken2_database \
    --taxonomy /mnt/raid1b/philip/databases/HOMD_16_03_kraken2/eHOMD_kraken2_database/taxonomy \
    --taxonomic_rank S \
    --bracken_threshold 10 \
    --bracken_length 600 \
    --out_dir $PROJECT_DIR/data/epi2me/kraken2_ehomd/species \
    --output_unclassified False \
    --threads 50




PROJECT_DIR="/mnt/raid1b/philip/dental_implants"
cd $PROJECT_DIR/data/epi2me
nextflow run epi2me-labs/wf-16s \
    -profile singularity \
    --fastq /mnt/raid1b/philip/dental_implants/data/raw/both_runs  \
    --min_len 600 \
    --classifier kraken2 \
    --database /mnt/raid1b/philip/databases/HOMD_16_03_kraken2/eHOMD_kraken2_database \
    --taxonomy /mnt/raid1b/philip/databases/HOMD_16_03_kraken2/eHOMD_kraken2_database/taxonomy \
    --taxonomic_rank G \
    --bracken_threshold 10 \
    --bracken_length 600 \
    --out_dir $PROJECT_DIR/data/epi2me/kraken2_ehomd/genus \
    --output_unclassified False \
    --threads 50


PROJECT_DIR="/mnt/raid1b/philip/dental_implants"
cd $PROJECT_DIR/data/epi2me
nextflow run epi2me-labs/wf-16s \
    -profile singularity \
    --fastq /mnt/raid1b/philip/dental_implants/data/raw/both_runs  \
    --min_len 600 \
    --classifier kraken2 \
    --database /mnt/raid1b/philip/databases/HOMD_16_03_kraken2/eHOMD_kraken2_database \
    --taxonomy /mnt/raid1b/philip/databases/HOMD_16_03_kraken2/eHOMD_kraken2_database/taxonomy \
    --taxonomic_rank F \
    --bracken_threshold 10 \
    --bracken_length 600 \
    --out_dir $PROJECT_DIR/data/epi2me/kraken2_ehomd/family \
    --output_unclassified False \
    --threads 50

PROJECT_DIR="/mnt/raid1b/philip/dental_implants"
cd $PROJECT_DIR/data/epi2me
nextflow run epi2me-labs/wf-16s \
    -profile singularity \
    --fastq /mnt/raid1b/philip/dental_implants/data/raw/both_runs  \
    --min_len 600 \
    --classifier kraken2 \
    --database /mnt/raid1b/philip/databases/HOMD_16_03_kraken2/eHOMD_kraken2_database \
    --taxonomy /mnt/raid1b/philip/databases/HOMD_16_03_kraken2/eHOMD_kraken2_database/taxonomy \
    --taxonomic_rank P \
    --bracken_threshold 10 \
    --bracken_length 600 \
    --out_dir $PROJECT_DIR/data/epi2me/kraken2_ehomd/phylum \
    --output_unclassified False \
    --threads 50
