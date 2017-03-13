Abigail Watson  
MSBI 32400 – Lab #8  

```sh

# Set up lab files
mkdir Lab8
cd Lab8
mkdir bin
mkdir src
mkdir results
mkdir docs
mkdir data

# Annotate to find interesting regions  
cd results  
java -Xmx2G -jar /data/snpEff/snpEff.jar eff - canon -noLog hg19 /data/variantCaller_out.170/ IonXpress_003/TSVC_variants.vcf > TSVC_variants.snpEff.vcf  
java -Xmx2G -jar /data/snpEff/SnpSift.jar annotate -noLog /data/snpEff/data/hg19/clinvar/ clinvar_20170130.vcf.gz TSVC_variants.snpEff.vcf > TSVC_variants.snpEff.clinvar.vcf  

# Do some quick filtering  
# These create chord diagram chromosome data blobs of gene networks
grep -v "^#" TSVC_variants.snpEff.clinvar.vcf | grep -v '0/0’
grep -v "^#" TSVC_variants.snpEff.clinvar.vcf | grep -v '0/0' | grep stop



```



### Select COSMIC in hotspot bed  

**Gene Name**  
TP53  

**Mutation Id**   
COSM10996  

**AA Mutation**  
p.E171* (Substitution - Nonsense)  

**CDS Mutation**   
c.511G>T (Substitution, position 511, G➞T)  

**GRCh38**  
17:7675101..7675101, view Ensembl Contig  

**COSMIC Genome Browser**  
17:7675101..7675101, view in COSMIC JBrowse  

**Ever confirmed somatic**  
Yes   

**FATHMM prediction**  
Pathogenic (score 0.99)  


#### Filter IGV view with gene list

COSM6223  


Gene Name:
EGFRMutation

Id:
COSM6223AA

Mutation:
p.E746_A750delELREA (Deletion - In frame)

CDS Mutation:
c.2235_2249del15 (Deletion) GRCh38:7:55174772..55174786, view Ensembl Contig

COSMIC Genome Browser:7:55174772..55174786, view in COSMIC JBrowse Ever confirmed somatic:Yes


FATHMM prediction:none (score 0.00)





-----------
#### Download EGFR hotspot sample  

chr7:55242464 - 55242482  
KEGG_NON_SMALL_CELL_LUNG_CANCER  

TP53 - chr17:7571720-7590868
TP53 - chr17:7571720-7578811
