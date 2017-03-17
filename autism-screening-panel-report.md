Patient:  GRCh37/hg19  
Autism Screening Panel  
March 17th, 2017  

______________________________________  
#### Autism Screening Panel  

We use Pevsner's Austism Panel and the NCBI Ideogram drawer to visualize the 101 genome markers recommended for an Autism Screening Panel.  Such visualizations are important in personalized medicine endeavors, to allow clinicians and patients to navigate the genome and understand how a person's genetic profile compares to known biomarkers for clinical conditions. 

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/AutismPanel.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/AutismPanel.png)    

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
#### Variant Effect Predictor (Autism Screening Panel)     

After numerous runs, the best result I could get had two dozen results, covering only 3 genes from the original list of 101 Autism Screening Panel.  

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-graphs-pevsner.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-graphs-pevsner.png)   

![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-stats-pevsner.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/vep-stats-pevsner.png)  

These genes were LAMC3, WFS1, and SCN1A, which are associated with the following proteins and macromolecules:    

- [WFS1 - wolframin ER transmembrane glycoprotein](https://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt=Graphics&list_uids=7466) 
- [LAMC3 - laminin subunit gamma 3](https://www.ncbi.nlm.nih.gov/gene/10319)    
- [SCN1A - sodium voltage-gated channel alpha subunit 1](https://www.ncbi.nlm.nih.gov/gene/6323)    

______________________________________  
#### rs1801208 - WFS1  

Most of the listed SNPs are benign; although rs1801208 has a high PolyPhen score, indicative of a known change in phenotypic expression.  

AGCACCCATGCAGAGCCCTACACGC[A/G]CAGGGCCCTGGCCACCGAGGTCACC   
Chromosome: 4:6301162  
Gene:WFS1 (GeneView)   
Functional Consequence: missense  
Clinical significance: Likely benign  
Validated: by 1000G,by cluster,by frequency,by hapmap  
Global MAF:A=0.0603/302   
HGVS: NC_000004.11:g.6302889G>A, NC_000004.12:g.6301162G>A, NG_011700.1:g.36313G>A, NM_001145853.1:c.1367G>A, NM_006005.2:c.1367G>A, NM_006005.3:c.1367G>A, NP_001139325.1:p.Arg456His, NP_005996.1:p.Arg456His, NP_005996.2:p.Arg456His, XM_017008586.1:c.1376G>A, XP_016864075.1:p.Arg459His  


![https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/chromosome4-rs1801208.png](https://raw.githubusercontent.com/awatson1978/msbi-32400-final-project/master/screenshots/chromosome4-rs1801208.png)  


______________________________________  
#### Analysis      

During analysis, we set the threshold of ExAC allele observations at the CDC's observed rate of prevelence of Autism.  The general thinking is that there are a large number of SNPs that we can rule out; because they're present in larger portions of the population than the disease under study.  Of the SNPs that have frequencies less than the prevelence of Autism, some of them are known to be autosomnal recessive; and to have mild forms when a person is a carrier, but not doubly recessive.  

______________________________________  
#### Clinical Diagnosis    

In interpreting these results, we have to think about what the GRCh37/hg19 Genome Assembly represents.  Unlike the 1000 Genome Projects, which provides the possibility of _averaged_ data; the GRCh37 and GRCh38 assemblies represent a normalized reference genome.  That is, it was assmbled piecemeal from thirteen anonymous volunteers from Buffalo, New York.  As such, if the genome was evenly divided between the 13 people (in actuality, it wasn't), and there is a 1/68 chance of having autism, then would be a 19% chance of the reference genome having genes from somebody having autism.  

So, what these three genes from our Autism Screening Panel actually represent is this:  of the thirteen anonymous volunteers who contributed to the GRCh37 reference genome, they collectively had 3 of the 101 genes from Pevsner's Autism Screening Panel.  

______________________________________  
#### Mechanism of Action        

Broadly speaking, the genes identified by the screening panel are responsible for producing the proteins that form the sodium-calcium ion pumps in laminin and transmembrane glycoproteins don't function properly.  And laminin and transmembrane glycoproteins are involved in the development and maintenance of the occular lens, and management of diabetes.  

______________________________________  
#### Risk Scores      

If we had been running a 23andMe genome through this pipeline, at this point we might calculate a risk score based on the percentage of people with those genes who are later diagnosed with autism.  Niavely, having variant markers on only 3/101 genes in the panel would indicate a 3% risk of having autism.  In reality, the gene effects are likely distributed on a logarithmic or power distribution; and the results would either be a smoking gun, or somewhere on a long tail.  

Additionally, we might run the genomes from the 1000 Genome Project through the Autism Screening Panel, to determine the specificity and sensitivity of the screening panel.  

______________________________________  
#### Familial Heritability Risks  

While it's unlikely that any of the thirteen anonymous volunteers involved in the GRCh37 assembly had autism themselves; it appears that they may be carrying a few recessive genes for autism.  If they had children with somebody else who had variants on LAMC3, WFS1, and SCN1A, their children may be at risk of having Ushers Syndrome, Wolfram Syndrom, Diabetes, and vision/hearing loss.  





______________________________________  
#### References   
[Human Genome Browser - hg19 assembly](https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19)  
[National Center for Biotechnology Information - GRCh37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)  
[DSM5 - Diagnostic Criteria - Autism Spectrum Disorder](https://www.autismspeaks.org/what-autism/diagnosis/dsm-5-diagnostic-criteria)  
[Centers for Disease Control and Prevention: Autism Spectrum Disorder](https://www.cdc.gov/ncbddd/autism/data.html)     
[The Supplementary Material Repository for Bioinformatics Data Skills](https://github.com/vsbuffalo/bds-files)   
[High Throughput Sequencing Libraries](http://www.htslib.org/)  
[1000 Genomes Browser](https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/)  
[UCSC Genome Bioinformatics: Sequnce and Annotation Downloads](http://hgdownload.cse.ucsc.edu/downloads.html)    
[Chromosome Visualization with D3.js](https://github.com/eweitz/ideogram)    
[HUGO Gene Nomenclature Committee](http://www.genenames.org/)  
[Burrows-Wheeler Alignment Tool](http://bio-bwa.sourceforge.net/bwa.shtml)    
[Varient Effect Predictor Results - USH2A](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=nnArmLFBLwGQbZB7-2936516)
[Varient Effect Predictor Results - Autism Panel](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Ticket?tl=pD6v05iE5wtpiykx)  
[Micro-RNA Binding Site Polymorphisms in the WFS1 Gene Are Risk Factors of Diabetes Mellitus](http://eds.a.ebscohost.com.proxy.uchicago.edu/eds/detail/detail?sid=3c31a350-68d8-431a-a872-3ed327fa02ce%40sessionmgr4009&vid=0&hid=4211&bdata=JnNpdGU9ZWRzLWxpdmUmc2NvcGU9c2l0ZQ%3d%3d#AN=26426397&db=mnh)  
[Children with Usher syndrome: mental and behavioral disorders](http://behavioralandbrainfunctions.biomedcentral.com/articles/10.1186/1744-9081-8-16)  
[An eQTL mapping approach reveals that rare variants in the SEMA5A regulatory network impact autism risk](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3690972/)  
[Mutations in both gene copies more common in autism](https://spectrumnews.org/news/mutations-in-both-gene-copies-more-common-in-autism/)  
[Detailed investiation of the role of common and low-frequency WFS1 variants in type 2 diabetes risk](https://www.ncbi.nlm.nih.gov/pubmed/?term=20028947%5Buid%5D)  
[Genetic variants and susceptibility to neurological complications following West Nile virus infection](https://www.ncbi.nlm.nih.gov/pubmed/21881118)  
[Characterization and expression of the laminin gamma3 chain: a novel, non-basement membrane-associated, laminin chain.](https://www.ncbi.nlm.nih.gov/pubmed/10225960)  
[Linkage of the gene for Wolfram syndrome to markers on the short arm of chromosome 4.](https://www.ncbi.nlm.nih.gov/pubmed/7987399)  
[A gene encoding a transmembrane protein is mutated in patients with diabetes mellitus and optic atrophy (Wolfram syndrome).](https://www.ncbi.nlm.nih.gov/pubmed/9771706)  
[Localization of a putative human brain sodium channel gene (SCN1A) to chromosome band 2q24.](https://www.ncbi.nlm.nih.gov/pubmed/8062593)  
[Autosomal dominant epilepsy with febrile seizures plus with missense mutations of the (Na+)-channel alpha 1 subunit gene, SCN1A.](https://www.ncbi.nlm.nih.gov/pubmed/11823106)  
[International Union of Pharmacology. XLVII. Nomenclature and structure-function relationships of voltage-gated sodium channels.](https://www.ncbi.nlm.nih.gov/pubmed/16382098)  
