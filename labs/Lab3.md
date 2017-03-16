Abigail Watson
MSBI 32400 Lab 3

```bash
cd /data
mkdir /data/lab3

# Set up project directories
cd lab6
mkdir bin  
mkdir data
mkdir doc   
mkdir results
mkdir src   

# Download Autism Gene Panel  
# http://www.bioinfbook.org/php/C9E3k  
wget http://www.bioinfbook.org/wiley/3e/chapter9/WebDocument_9-05_101AutismPanel.bed.txt  

# Download Autism BAM File  
# http://www.bioinfbook.org/php/C9E3k  
wget http://www.bioinfbook.org/wiley/3e/chapter9/WebDocument_9-7_mysample1.bam  

# Let’s get the genes from the Gene File  
cut -f4 WebDocument_9-05_101AutismPanel.bed.txt > genelist.txt  

# Convert to BED file (using output type 'exon + 10bp')  
http://genome.ucsc.edu/cgi-bin/hgTables  

# Open IGV then add BAM + BEDs   

# Installing bedtools (as non-root)  
# http://bedtools.readthedocs.io/en/latest/content/installation.html    
wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz  
tar -zxvf bedtools-2.26.0.tar.gz  
mkdir bedtools2    
cd bedtools2  
make  

cp -rp bin/ ~  
export PATH=~/bin:$PATH  

# Let’s clean up the BED file  
bedtools sort -i 101AutismGenelistExons.bed > 101AutismGenelistExons_sort.bed    
bedtools merge -c 4 -o collapse -i 101AutismGenelistExons_sort.bed > 101AutismGenelistExons_sort_merged.bed  

# Load merged BED in IGV

cp /data/bds-files/chapter-11-alignment/NA12891_CEU_sample.vcf.gz lab3/data  
zcat lab3/data/NA12891_CEU_sample.vcf.gz | more  
gunzip lab3/data/NA12891_CEU_sample.vcf.gz  

more NA12891_CEU_sample.vcf | grep -v “^#” | wc –l

# SAMtools  
# Generate a FastQ File from a BAM file  
samtools view -H WebDocument_9-7_mysample1.bam | grep '@RG'
samtools fastq WebDocument_9-7_mysample1.bam > WebDocument_9-7_mysample1.fastq  
samtools view –h WebDocument_9-7_mysample1.bam > WebDocument_9-7_mysample1.sam  

# Create a new BAM file from SAM
samtools view -bS WebDocument_9-7_mysample1.sam | samtools sort -o WebDocument_9-7_mysample1_file_sorted.bam   
samtools index WebDocument_9-7_mysample1_file_sorted.bam   
samtools view -H WebDocument_9-7_mysample1_file_sorted.bam  

```
