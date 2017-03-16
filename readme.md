Abigail Watson  
MSBI 32400 – Final Project  

______________________________________
#### Project Setup    

```bash
mkdir msbi-32400-final-project  
cd msbi-32400-final-project  
mkdir bin  
mkdir data
mkdir doc   
mkdir results
mkdir src   
```

______________________________________
#### Assignment #1    

Call variants on Pevsner autism bam (Web Document 9.7 at http://www.bioinfbook.org/php/ then annotate with snpEff + Clinvar and upload to VEP for ExAC population frequencies.  Compare variants in the 101 target gene list with ExAC frequencies expected autism frequencies.  Use hg19.fa (see below), not Pevsner's WebDocument_9-6_101autism.fa.  Use CDC 2016 frequencies: https://www.cdc.gov/ncbddd/autism/ 

______________________________________
#### CDC Statistics on Autism Spectrum Disorder  

To understand the genetic factors contributing to Autism Spectrum Disorder, we first begin looking at phenotype expression, which we can obtain from the Centers of Disease Control.  The difficulty in doing so is that the CDC reports that the prevalence of Autism has been increasing over the past few decades.  In the Surveillance Year 2000, it was reported that 1 in 150 children who were born in 1992 were diagnosed with autism (0.006%); whereas in 2012, the numbers had increased to 1 in 68 children who were born in 2004 were diagnosed with autism (0.0147%).  This presents a few hypotheses:

A.  Autism is increasing in the general population.   
  A1. Environmental factors such as pesticides and plastics are causing an increase in Autism.  
  A2. Autism is a side effect of some other benefitial process, such as change in diet.  
B.  The reporitng process is including more autisitc children.  

We also note that parents who have a child with ASD have a 2%–18% chance of having a second child who is also affected, which suggests a recessive trait since it's a less than 1/4.  ASD tends to occur more often in people who have certain genetic or chromosomal conditions, indicating involvement in a gene network (to be expected for behavioral and psychiatric disorders).  Children born to older parents are at a higher risk for having ASD, indicating that ASD is a common symptom of mutations in the germ line.  ASD commonly co-occurs with other developmental, psychiatric, neurologic, chromosomal, and genetic diagnoses in 83% of the cases; again suggesting ASD is a common symptom of regulatory pathway disfunction.  

______________________________________  
#### Tool Installation      
Note:  These instructions are for a Mac.  Please adjust use of wget/curl and compiled binaries accordingly.  

```bash
# Install the twoBitToFa utility  
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa  
mv ~/Downloads/twoBitToFa msbi-32400-final-project/bin   
cd msbi-32400-final-project/bin      
chmod +x twoBitToFa  


# Download and install samtools
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2  
mv ~/Downloads/samtools-1.3.1.tar.bz2  msbi-32400-final-project/src  
cd msbi-32400-final-project/src  
gunzip samtools-1.3.1.tar.bz2  
tar -xzvf samtools-1.3.1.tar  
cd samtools-1.3.1  
make  
make ../../bin/samtools install  


# Download and install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2  
mv ~/Downloads/bcftools-1.3.1.tar.bz2  msbi-32400-final-project/src   
cd msbi-32400-final-project/src  
gunzip bcftools-1.3.1.tar.bz2  
tar -xzvf bcftools-1.3.1.tar   
cd bcftools-1.3.1  
make  
make ../../bin/bcftools install  


# Install the Burrows-Wheeler Aligner tool  
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2  
mv ~/Downloads/bwa-0.7.15.tar.bz2  msbi-32400-final-project/src  
cd msbi-32400-final-project/src  
gunzip bwa-0.7.15.tar.bz2  
tar -xzvf bwa-0.7.15.tar.bz2  
cd bwa-0.7.15  
make  
mv ./bwa ../../bin  

# Download and install SnpEff  
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip  
mv ~/Downloads/snpEff_latest_core.zip  msbi-32400-final-project/src    
cd msbi-32400-final-project/src  
unzip snpEff_latest_core.zip  
mv snpEff ../bin 
```

______________________________________  
#### Data Files        
```bash
# Download Autism BAM File  
# http://www.bioinfbook.org/php/C9E3k  
cd data  
wget http://www.bioinfbook.org/wiley/3e/chapter9/WebDocument_9-7_mysample1.bam   
mv ~/Downloads/WebDocument_9-7_mysample1.bam msbi-32400-final-project/data  

# Download the complete hg19 1000 Genome Project in the 2bit format
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit  
mv ~/Downloads/hg19.2bit msbi-32400-final-project/data  

# Download and install gene USH2A from Vince Buffalo Chapter 11  
git clone https://github.com/vsbuffalo/bds-files  
cp ~/Downloads/bds-files/chapter-11-alignment/NA12891_CEU_sample.bam msbi-32400-final-project/data

# Alternatively, if we had the raw GRCh37 Illumina Data, we could extract the USH2A data with the following:   
# samtools view -hb ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12891/alignment/NA12891.chrom1.ILLUMINA.bwa.CEU.high_coverage.20100517.bam 1:215622894-216423396 > NA12891_CEU_sample.bam

# Get the Feb 2017 ClinVar VCF file  
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
mv ~/Downloads/clinvar.vcf.gz* msbi-32400-final-project/data 
```

______________________________________  
#### Data Pipeline for   

```bash
# Convert the raw Human Genome data (hg19.2bit) into a .fa file  
cd msbi-32400-final-project
./bin/twoBitToFa ./data/hg19.2bit ./data/hg19.fa  

# Index the hg19 Human Genome with SAMTools; this will produce our .fai file
cd msbi-32400-final-project  
bin/samtools/bin/samtools faidx data/hg19.fa  

# Create our fastq file
cd msbi-32400-final-project  
bin/samtools/bin/samtools fastq data/NA12891_CEU_sample.bam > data/NA12891_CEU_sample.fastq

# Use the Burrows-Wheeler Aligner to generate our SAM file
cd msbi-32400-final-project  
bin/bwa index -a bwtsw data/hg19.fa
bin/bwa mem -R '@RG\tID:MSBI32400_Abigail_Watson_final_project\tSM:NA12891_CEU_sample' data/hg19.fa data/NA12891_CEU_sample.fastq > data/NA12891_CEU_sample.sam  

# view, sort, and index all the things
bin/samtools/bin/samtools view -bt data/hg19.fai data/NA12891_CEU_sample.sam > data/NA12891_CEU_sample.bam
bin/samtools/bin/samtools sort -o data/NA12891_CEU_sample_sorted.bam data/NA12891_CEU_sample.bam
bin/samtools/bin/samtools index data/NA12891_CEU_sample_sorted.bam
bin/samtools/bin/samtools view -H data/NA12891_CEU_sample_sorted.bam

# Generate mpileup & run bcftools  
bin/samtools/bin/samtools mpileup -uf data/hg19.fa data/NA12891_CEU_sample_sorted.bam | bin/bcftools call -mv > data/NA12891_CEU_sample_sorted_var.raw.vcf
bin/bcftools filter -s LowQual -e "%QUAL<20" data/NA12891_CEU_sample_sorted_var.raw.vcf > data/NA12891_CEU_sample_sorted_var.flt.vcf

# Count the lines
more data/NA12891_CEU_sample_sorted_var.flt.vcf | grep 'chr1' | wc -l     # 756  
more data/NA12891_CEU_sample_sorted_var.flt.vcf | grep 'PASS' | wc -l     # 686 - without the "PASS" in the header

# Analyze with snpEff + Clinvar Feb 2017
cd msbi-32400-final-project/data  
java -Xmx4g -jar ../bin/snpEff/snpEff.jar eff -v -d -canon -noLog hg19 NA12891_CEU_sample_sorted_var.flt.vcf > NA12891_CEU_sample_sorted_var.flt.snpEff.vcf

java -Xmx2G -jar ../bin/snpEff/SnpSift.jar annotate -noLog clinvar.vcf.gz NA12891_CEU_sample_sorted_var.flt.snpEff.vcf > NA12891_CEU_sample_sorted_var.flt.snpEff.clinvar.vcf

java -Xmx2G -jar ../bin/snpEff/SnpSift.jar extractFields -s ',' -e '.' NA12891_CEU_sample_sorted_var.flt.snpEff.clinvar.vcf CHROM POS REF ALT ID "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" CLNHGVS CLNALLE CLNACC CLNSIG CLNREVSTAT CLNDBN > NA12891_CEU_sample_sorted_var.flt.snpEff.clinvar.Extracted
```
______________________________________  
#### Pipeline Results  

Resulting files after the pipeline was completed:  
![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/DataFileSizes.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/DataFileSizes.png)   

______________________________________  
#### ANNOVAR Analysis    

After running the data analysis pipeline, we have have a VCF file that we can upload into ANNOVAR and VEP.  Annovar has consistently provided more reliable in my experience, so analysis begun there.  Of particular interest is that the exome summary results are narrowed down to only 9 single necleotide polymorphisms.  Furthermore, sorting by ExAC Frequency, and we see that 8 of the 9 SNPs are present in 13% or more of the population, with three of them being present in over half the human population.  However, one of them, rs45549044, is present in only 0.0043% of the population, which is very close the observed prevelance of 0.006% that the CDC recorded in the year 2000.  And while it's noted as being non-pathogenic, it is the only variant that warrents a COSMIC ID, which happens to be COSM1338810.

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/wannovar-rs45549044.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/wannovar-rs45549044.png)  
______________________________________  
#### Variant Effect Predictor  

Running the VCF file through the Varient Effect Predictor failed the first half-dozen tries.  Eventually, a result was achieved that identified six variants of interest on three SNPs: 

rs11120616  
rs35309576  
rs45549044  

Furthermore, not only only did VEP also identify rs45549044 as having moderate impact, it's PolyPhen score was 0.934, indicating that it is a probably damaging varient.  

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-graphs-rs45549044.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-graphs-rs45549044.png)  

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-stats-rs45549044.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-stats-rs45549044.png)  


______________________________________  
#### rs45549044  

```
GAATCTATAAAAGATGTTGAGCTTC[C/T]GTTATAGATTAGGACTGGATTGGAT  
Chromosome: 1:215671031  
Gene:USH2A (GeneView)   
Functional Consequence: missense  
Clinical significance: Benign  
Validated: by 1000G,by cluster,by frequency  
Global MAF:T=0.0022/11  
HGVS: NC_000001.10:g.215844373C>T, NC_000001.11:g.215671031C>T, NG_009497.1:g.757366G>A, NM_206933.2:c.14074G>A, NP_996816.2:p.Gly4692Arg  
```

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/chromsome1-rs45549044.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/chromsome1-rs45549044.png)  

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/rs45549044.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/rs45549044.png)   

______________________________________  
#### Autism Screening Panel  

We use Pevsner's Austism Panel and the NCBI Ideogram drawer to visualize the 101 genome markers recommended for an Autism Screening Panel.  Such visualizations are important in personalized medicine endeavors, to allow clinicians and patients to navigate the genome and understand how a person's genetic profile compares to known biomarkers for clinical conditions. 

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/AutismPanel.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/AutismPanel.png)    

______________________________________  
#### Cancer Risks    

Interestingly, mutations on the rs45549044 SNP are [associated with carcinomas of the large intestine](http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=1338810).  Which means that people diagnosed with ASD may want to undergo the Autism Screening Panel to screen for colon cancer.  

______________________________________  
#### Discussion  

Returning to the original hypothesis regarding observed CDC rates...

______________________________________  
#### References   
https://www.cdc.gov/ncbddd/autism/data.html   
https://github.com/vsbuffalo/bds-files   
http://www.htslib.org/  
http://hgdownload.cse.ucsc.edu/downloads.html  
http://hgdownload.soe.ucsc.edu/admin/exe/  
https://github.com/eweitz/ideogram  

#### Other Research   
https://wikis.utexas.edu/display/bioiteam/Variant+calling+using+SAMtools  
http://bio-bwa.sourceforge.net/bwa.shtml  
https://www.ebi.ac.uk/sites/ebi.ac.uk/files/content.ebi.ac.uk/materials/2014/140217_AgriOmics/dan_bolser_snp_calling.pdf  
