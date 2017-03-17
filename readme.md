Abigail Watson  
MSBI 32400 – Final Project  


______________________________________
## Assignment   

Call variants on Pevsner autism bam (Web Document 9.7 at http://www.bioinfbook.org/php/ then annotate with snpEff + Clinvar and upload to VEP for ExAC population frequencies.  Compare variants in the 101 target gene list with ExAC frequencies expected autism frequencies.  Use hg19.fa (see below), not Pevsner's WebDocument_9-6_101autism.fa.  Use CDC 2016 frequencies: https://www.cdc.gov/ncbddd/autism/ 

______________________________________
#### CDC Statistics on Autism Spectrum Disorder  

To understand the genetic factors contributing to Autism Spectrum Disorder, we first begin looking at phenotype expression, which we can obtain from the Centers of Disease Control.  The difficulty in doing so is that the CDC reports that the prevalence of Autism has been increasing over the past few decades.  In the Surveillance Year 2000, it was reported that 1 in 150 children who were born in 1992 were diagnosed with autism (0.6%); whereas in 2012, the numbers had increased to 1 in 68 children who were born in 2004 were diagnosed with autism (1.47%).  This presents a few hypotheses:

A.  Autism is increasing in the general population.   
  A1. Environmental factors such as pesticides and plastics are causing an increase in Autism.  
  A2. Autism is a side effect of some other benefitial process, such as change in diet.  
B.  The reporitng process is including more autisitc children.  
  B1. The screening process for diagnosing autistic children is becoming better at detecting autism.
  B2. The definition of autism has changed and become more general.

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

The following data files were used for this assignment.  

```bash
# Download the complete hg19 Human Genome in the 2bit format
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit  
mv ~/Downloads/hg19.2bit msbi-32400-final-project/data  

# Download and install gene USH2A from Vince Buffalo Chapter 11  
git clone https://github.com/vsbuffalo/bds-files  
cp ~/Downloads/bds-files/chapter-11-alignment/NA12891_CEU_sample.bam msbi-32400-final-project/data

# Get the Feb 2017 ClinVar VCF file  
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
mv ~/Downloads/clinvar.vcf.gz* msbi-32400-final-project/data 

# Download Autism BAM File  
# http://www.bioinfbook.org/php/C9E3k  
cd data  
wget http://www.bioinfbook.org/wiley/3e/chapter9/WebDocument_9-7_mysample1.bam   
mv ~/Downloads/WebDocument_9-7_mysample1.bam msbi-32400-final-project/data  
```

______________________________________  
#### Data Preparaton  

The following steps are necessary for all data pipelines (single gene, screening panel, personalized medicine).  

```bash
# Convert the raw Human Genome data (hg19.2bit) into a .fa file  
cd msbi-32400-final-project
./bin/twoBitToFa ./data/hg19.2bit ./data/hg19.fa  

# Index the hg19 Human Genome with SAMTools; this will produce our .fai file
cd msbi-32400-final-project  
bin/samtools/bin/samtools faidx data/hg19.fa  
```     

______________________________________  
#### Initial Data Pipeline (Single Gene - USH2A)  

```bash
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
#### Pipeline Results (Single Gene - USH2A)    

Resulting files after the pipeline was completed:  
![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/DataFileSizes.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/DataFileSizes.png)   

______________________________________  
#### ANNOVAR Analysis (Single Gene - USH2A)      

After running the data analysis pipeline, we have have a VCF file that we can upload into ANNOVAR and VEP.  Annovar has consistently provided more reliable in my experience, so analysis begun there.  Of particular interest is that the exome summary results are narrowed down to only 9 single necleotide polymorphisms.  Furthermore, sorting by ExAC Frequency, and we see that 8 of the 9 SNPs are present in 13% or more of the population, with three of them being present in over half the human population.  However, one of them, rs45549044, is present in only 0.0043% of the population, which is very close the observed prevelance of 0.006% that the CDC recorded in the year 2000.  And while it's noted as being non-pathogenic, it is the only variant that warrents a COSMIC ID, which happens to be COSM1338810.

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/wannovar-rs45549044.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/wannovar-rs45549044.png)  
______________________________________  
#### Variant Effect Predictor (Single Gene - USH2A)    

Running the VCF file through the Varient Effect Predictor failed the first half-dozen tries.  Eventually, a result was achieved that identified six variants of interest on three SNPs: 

rs11120616  
rs35309576  
rs45549044  

Furthermore, not only only did VEP also identify rs45549044 as having moderate impact, it's PolyPhen score was 0.934, indicating that it is a probably damaging varient.  

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-graphs-rs45549044.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-graphs-rs45549044.png)  

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-stats-rs45549044.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-stats-rs45549044.png)  


______________________________________  
#### rs45549044 - USH2A  

Using the example BAM from chapter 11 of the Vince Buffalo book, we confirm the VEP pipeline identifies genes and SNPs of interest.  In this case, we identify USH2A, which is located on Chromsome 1 at 1q41, encodes the protein Usherin, and which is involved in hearing and vision.   


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
#### Cancer Risks    

Interestingly, mutations on the rs45549044 SNP associated with USH2A are [associated with carcinomas of the large intestine](http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=1338810).  Which means that people diagnosed with ASD may want to undergo the Autism Screening Panel to screen for colon cancer.  

______________________________________  
#### Autism Screening Panel  

We use Pevsner's Austism Panel and the NCBI Ideogram drawer to visualize the 101 genome markers recommended for an Autism Screening Panel.  Such visualizations are important in personalized medicine endeavors, to allow clinicians and patients to navigate the genome and understand how a person's genetic profile compares to known biomarkers for clinical conditions. 

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/AutismPanel.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/AutismPanel.png)    


______________________________________  
#### Revised Pipeline (Autism Screening Panel)  

After running through this assignment, an article from Pychology Today entitled [Harvard Study Finds Genetic 'Toggle Switch' for Sociability](https://www.psychologytoday.com/blog/the-athletes-way/201703/harvard-study-finds-genetic-toggle-switch-sociability) ran across my Facebook feed; which detailed the role of UBE3A in downregulating a glutamate based synapse organizer called CBLN1.  I had just ran Pevsner's Autism Panel through the NCBI Ideogram tool; so I took a peek at the screening panel to see if UBE3A was in there.  Sure enough, there it was.  So then it occurred to me that maybe I should try running the entire Autism Panel against the the hg19 Human Genome.  Because that would be pretty cool.  

But where to get the BAM file for the entire Autism Panel?  And that's when I remembered that I had used the NA12891_CEU_sample file from Chapter 11, but I hadn't used the WebDocument_9-7_mysample1.bam file.  

```bash
# Create our fastq file
cd msbi-32400-final-project  
bin/samtools/bin/samtools fastq data/WebDocument_9-7_mysample1.bam > data/WebDocument_9-7_mysample1.fastq

# Use the Burrows-Wheeler Aligner to generate our SAM file
cd msbi-32400-final-project  
bin/bwa index -a bwtsw data/hg19.fa
bin/bwa mem -R '@RG\tID:MSBI32400_Abigail_Watson_final_project\tSM:WebDocument_9-7_mysample1' data/hg19.fa data/WebDocument_9-7_mysample1.fastq > data/WebDocument_9-7_mysample1.sam  

# view, sort, and index all the things
bin/samtools/bin/samtools view -bt data/hg19.fai data/WebDocument_9-7_mysample1.sam > data/WebDocument_9-7_mysample1.bam
bin/samtools/bin/samtools sort -o data/WebDocument_9-7_mysample1_sorted.bam data/WebDocument_9-7_mysample1.bam
bin/samtools/bin/samtools index data/WebDocument_9-7_mysample1_sorted.bam
bin/samtools/bin/samtools view -H data/WebDocument_9-7_mysample1_sorted.bam

# Generate mpileup & run bcftools  
bin/samtools/bin/samtools mpileup -uf data/hg19.fa data/WebDocument_9-7_mysample1_sorted.bam | bin/bcftools call -mv > data/WebDocument_9-7_mysample1_sorted_var.raw.vcf
bin/bcftools filter -s LowQual -e "%QUAL<20" data/WebDocument_9-7_mysample1_sorted_var.raw.vcf > data/WebDocument_9-7_mysample1_sorted_var.flt.vcf

# Count the lines
more data/WebDocument_9-7_mysample1_sorted_var.flt.vcf | grep 'chr1' | wc -l     # 11179 !
more data/WebDocument_9-7_mysample1_sorted_var.flt.vcf | grep 'PASS' | wc -l     # 14525 !

# Analyze with snpEff + Clinvar Feb 2017
cd msbi-32400-final-project/data  
java -Xmx4g -jar ../bin/snpEff/snpEff.jar eff -v -d -canon -noLog hg19 WebDocument_9-7_mysample1_sorted_var.flt.vcf > WebDocument_9-7_mysample1_sorted_var.flt.snpEff.vcf

java -Xmx2G -jar ../bin/snpEff/SnpSift.jar annotate -noLog clinvar.vcf WebDocument_9-7_mysample1_sorted_var.flt.snpEff.vcf > WebDocument_9-7_mysample1_sorted_var.flt.snpEff.clinvar.vcf

java -Xmx2G -jar ../bin/snpEff/SnpSift.jar extractFields -s ',' -e '.' WebDocument_9-7_mysample1_sorted_var.flt.snpEff.clinvar.vcf CHROM POS REF ALT ID "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" CLNHGVS CLNALLE CLNACC CLNSIG CLNREVSTAT CLNDBN > WebDocument_9-7_mysample1_sorted_var.flt.snpEff.clinvar.Extracted
```
______________________________________  
#### Pipeline Results (Autism Screening Panel)   

For the second run, we now have the the AutismPanel aligned, indexed, sorted, piledup, and analyzed.  

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/DataFileSizes-PevsnerAutismFile.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/DataFileSizes-PevsnerAutismFile.png)  


______________________________________  
#### WANNOVAR Analysis (Autism Screening Panel)   

Once again, the WANNOVAR tool seems to have a better auto-configuration, and does a better job with returning usable results with minimal trial-and-error.  It succesfully returned thousands of results from the Pevsner Autism File, enough that I was able to specify a threshold of 0.0147 corresponding to the CDC's observed prevelance of Autism, and drill down to multiple results for the UBE3A gene.  This confirms that the VCF file is correctly parsed, that we're analyzing the complete 101 gene autism screening panel, and we're correctly reporting the ExAC frequency numbers.

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/wannovar-PevsnerAutismFile-filtered.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/wannovar-PevsnerAutismFile-filtered.png)  



______________________________________  
#### Variant Effect Predictor (Autism Screening Panel)     

After numerous runs, the best result I could get had 23 results, covering only 4 genes from the original list of 101 Autism Screening Panel.  These genes were LAMC3, WFS1, SCN1A, and LAMC3.  


![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-graphs-pevsner.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-graphs-pevsner.png)   

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-stats-pevsner.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-stats-pevsner.png)  







______________________________________  
#### Discussion  

Returning to the original hypothesis regarding observed CDC rates...


Personally, I find the Diagnostic and Statistical Manual of Mental Disorders to have a lot of shoddy science in it, and consider psychiatry to be something of a pseudoscience.  Years from now, we'll look back at it and compare it to the study of Phylogeny.  Yet, Phylogeny did lead to the field of Cladistics.  Similarly, psychiatry and the DSM is helping us tease apart genetics and neurology.

Having said that, after reviewing the genetics of Autism, here's my general opinion on the situation:  I think the definition of Autism has changed over time.  At one point in time, it was a diagnostic description of a set of behaviors that may have been generally associated with a recessive trait, such as Usher's Syndrome, that was associated with the regulatory pathways that control the communication organs (speach, hearing, sight).  The general thinking is that there are a large number of SNPs that we can rule out; because they're present in larger portions of the population than the disease under study.  Of the SNPs that have frequencies less than the prevelence of Autism, some of them are known to be autosomnal recessive; and to have mild forms when a person is a carrier, but not doubly recessive.  

Here's the thing though:  the definition of Autism has been shifting away from an underlying genetic determinant that's correlated with a regulatory pathway that controls communication and cognition.  Instead, it's shifting towards a metabolic or neurological mechanism that describes that regulatory pathway; and it's now including any genetic determinants that may cause a chance in that metabolic or neurological mechanism.







______________________________________  
#### References   
[DSM5 - Diagnostic Criteria - Autism Spectrum Disorder](https://www.autismspeaks.org/what-autism/diagnosis/dsm-5-diagnostic-criteria)  
[Centers for Disease Control and Prevention: Autism Spectrum Disorder](https://www.cdc.gov/ncbddd/autism/data.html)     
[The Supplementary Material Repository for Bioinformatics Data Skills](https://github.com/vsbuffalo/bds-files)   
[High Throughput Sequencing Libraries](http://www.htslib.org/)  
[UCSC Genome Bioinformatics: Sequnce and Annotation Downloads](http://hgdownload.cse.ucsc.edu/downloads.html)    
[Chromosome Visualization with D3.js](https://github.com/eweitz/ideogram)    
[HUGO Gene Nomenclature Committee](http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=12601)  
[Burrows-Wheeler Alignment Tool](http://bio-bwa.sourceforge.net/bwa.shtml)    
[Varient Effect Predictor Results - USH2A](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=nnArmLFBLwGQbZB7-2936516)    
[Varient Effect Predictor Results - Autism Panel](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Ticket?tl=pD6v05iE5wtpiykx)  
[Micro-RNA Binding Site Polymorphisms in the WFS1 Gene Are Risk Factors of Diabetes Mellitus](http://eds.a.ebscohost.com.proxy.uchicago.edu/eds/detail/detail?sid=3c31a350-68d8-431a-a872-3ed327fa02ce%40sessionmgr4009&vid=0&hid=4211&bdata=JnNpdGU9ZWRzLWxpdmUmc2NvcGU9c2l0ZQ%3d%3d#AN=26426397&db=mnh)  
[Children with Usher syndrome: mental and behavioral disorders](http://behavioralandbrainfunctions.biomedcentral.com/articles/10.1186/1744-9081-8-16)  
[An eQTL mapping approach reveals that rare variants in the SEMA5A regulatory network impact autism risk](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3690972/)    
[Mutations in both gene copies more common in autism](https://spectrumnews.org/news/mutations-in-both-gene-copies-more-common-in-autism/)  
 
[Detailed investiation of the role of common and low-frequency WFS1 variants in type 2 diabetes risk](https://www.ncbi.nlm.nih.gov/pubmed/?term=20028947%5Buid%5D)  

[Genetic variants and susceptibility to neurological complications following West Nile virus infection](Genetic variants and susceptibility to neurological complications following West Nile virus infection.)  



