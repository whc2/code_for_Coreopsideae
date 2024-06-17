# Introduction
These scripts were used in the genome project of three plants in the tribe Coreopsideae within the sunflower family (Asteraceae). All the whole-genome sequencing data that support this project have been deposited at China National Genomics Data Center (https://ngdc.cncb.ac.cn) with the Project ID PRJCA017572. The genome assemblies, gene annotations, and other resources are also available at Zenodo (https://doi.org/10.5281/zenodo.8296602).

# Contigs coverage
```
## 1. reads mapping
minimap2 -t 90 -ax map-hifi -I 16G ctgs.fa reads1.fq.gz | samtools view -@ 20 -bo reads1.fq.gz.bam -
minimap2 -t 90 -ax map-hifi -I 16G ctgs.fa reads2.fq.gz | samtools view -@ 20 -bo reads2.fq.gz.bam -

## 2. merge bam files, filter unmapped and non-primary alignments, and sort bam
samtools merge -@ 40 ccs_bam_merged.bam reads1.fq.gz.bam reads2.fq.gz.bam 
samtools view -F 0x104 -@ 40 ccs_bam_merged.bam -o ccs_bam_merged.noSecond_noUnmap.bam
samtools sort -o ccs_bam_merged.noSecond_noUnmap.sort.bam -@ 80 ccs_bam_merged.noSecond_noUnmap.bam

## 3. calculate base-level coverage
bedtools genomecov -ibam ccs_bam_merged.noSecond_noUnmap.sort.bam -bg > ccs_bam.bedGraph

## 4. calculate contig-level coverage
perl contig_coverage.pl ctg.fa.len ccs_bam.bedGraph > ccs_bam.bedGraph.ctgCov
```

# Insertion time of estimation for intact LTRs
```
## 1. generate bed file of each LTR
perl ltr_extract.pl genome.EDTA.intact.gff3

## 2. get fasta of each LTR (bedtools v2.30.0), use data parallel
bedtools getfasta -fi genome.fa -bed LTR_Copia_LTR_100000.bed -fo LTR_Copia_LTR_100000.bed.fa -s

## 3. multiple sequence alignment of each LTR, use data parallel
muscle -in LTR_Copia_LTR_100000.bed.fa -out LTR_Copia_LTR_100000.bed.fa.mucle

## 4. move different types of LTRs into different directory
mv LTR_Copia*.muscle z.Copia/

## 5. estimate insertion time
R CMD BATCH insertion_time.R insertion_time.Rout

## 6. summary insertion times
perl sum_times.pl LTR_Copia_Insert_time.csv LTR_Gypsy_Insert_time.csv LTR_unknown_Insert_time.csv > LTR_Insert_sum.tsv

## 7. draw insertion density plot
Rscript intact_LTR_density.R LTR_Insert_sum.tsv
```

# Calculatation of Transcripts Per Millions (TPMs) based on full length transrciptome
```
## 1. reads mapping
minimap2 -t 48 -ax map-hifi Bidens_alba.nucleus.chr.gene.V1.cds.gz full_length.transcriptome.ccs.fq.gz > Bidens.isoseq.minimap2.sam 2> Bidens.isoseq.minimap2.sam.err

## 2. stat mapping result
perl minimap2_sam2gff.pl Bidens.isoseq.minimap2.sam > Bidens.isoseq.minimap2.sam.tsv 2> Bidens.isoseq.minimap2.sam.tsv.stat

## 3. calculate TPM
perl tsv_to_expression.pl Bidens.isoseq.minimap2.sam.tsv.stat Bidens_alba.nucleus.chr.gene.V1.cds.len Bidens.isoseq.minimap2.sam.tsv > Bidens.isoseq.minimap2.sam.tsv.TPM

```


## Reference
> Wang H#, Xu D#, Jiang F, Wang S, Wang A, Liu H, Lei L, Q W, Fan W*. The genomes of *Dahlia pinnata*, *Cosmos bipinnatus*, and *Bidens alba* in tribe Coreopsideae provide insights into polyploid evolution and inulin biosynthesis. GigaScience. 2024 June 13:giae032. [https://doi.org/10.1093/gigascience/giae032](https://doi.org/10.1093/gigascience/giae032)
