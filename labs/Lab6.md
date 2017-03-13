Abigail Watson  
MSBI 32400 Lab 6


```sh
cd /data
mkdir /data/lab6

# Set up project directories
cd lab6
mkdir bin  
mkdir data 
mkdir doc   
mkdir results
mkdir src  

# Create our fastq file
cd /data/bds-files/ chapter-11-alignment/
samtools fastq ./NA12891_CEU_sample.bam > /data/lab6/data/NA12891_CEU_sample.fastq

# Get our chromosome genome
wget hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz
mv ~/Downloads/chr1.ga.gz /data/lab6/data
cd /data/lab6/data
gunzip chr1.fa.gz

# Need to index for bwa
bwa index -a bwtsw chr1.fa
 
bwa mem -R '@RG\tID:MSBI32400_test\tSM:NA12891_CEU_sample' chr1.fa NA12891_CEU_sample.fastq > NA12891_CEU_sample.sam  

# view, sort, and index all the things
samtools faidx chr1.fa
samtools view -bt chr1.fa.fai NA12891_CEU_sample.sam > NA12891_CEU_sample.bam
samtools sort -o NA12891_CEU_sample_sorted.bam NA12891_CEU_sample.bam
samtools index NA12891_CEU_sample_sorted.bam
samtools view -H NA12891_CEU_sample_sorted.bam

# read the fscking manual
man samtools

# Generate mpileup & run bcftools  
samtools mpileup -uf chr1.fa NA12891_CEU_sample_sorted.bam | bcftools call -mv > NA12891_CEU_sample_sorted_var.raw.vcf
bcftools filter -s LowQual -e ’%QUAL<20’ NA12891_CEU_sample_sorted_var.raw.vcf > NA12891_CEU_sample_sorted_var.flt.vcf

# Count the lines
more NA12891_CEU_sample_sorted_var.flt.vcf | grep 'chr1' | wc -l   
# 759  

more NA12891_CEU_sample_sorted_var.flt.vcf | grep 'PASS' | wc -l  
# 701  
# 700 - correct answer, there is one "PASS" in the header

# USH2A
chr1:215796236-216596738  
id=NM_206933  
Exon number: 64  
Amino acid coding number: 4698  
chr1:215844314-215844635  

# Vince's Way  
samtools mpileup -v --no-BAQ --region chr1:215622894-216423396 --fasta-ref chr1.fa NA12891_CEU_sample_sorted.bam > NA12891_CEU_sample_sorted_full_region.vcf.gz

```

The -x, -B, and -s options modify the data which is contained in each alignment.

The  bcftools  filter command marks low quality sites and sites with the read depth exceeding a
limit, which should be adjusted to about twice the average read depth (bigger read depths  usu-
ally indicate problematic regions which are often enriched for artefacts).  One may consider to
add -C50 to mpileup if mapping quality is overestimated for  reads  containing  excessive  mis-
matches. Applying this option usually helps BWA-short but may not other mappers.

