Abigial Watson  
MSBI 32400 – Lab 7   


```bash
cd Lab7
mkdir bin
mkdir data
mkdir docs
mkdir results
mkdir src


# Update Clinvar
cd data
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi

# Get the testdata
cd data
wget http://raw.githubusercontent.com/personalis/hgvslib/master/example_test_set/hgvs_test_cases.vcf --no-check-certificate

# Analyze with snpEff + Clinvar Feb 2017
cd ..
java -Xmx2G -jar ../Lab7/bin/snpEff/snpEff.jar eff -canon -noLog hg19 ./data/hgvs_test_cases.vcf > ./results/hgvs_test_cases_snpEff.vcf
java -Xmx2G -jar ../Lab7/bin/snpEff/SnpSift.jar annotate -noLog ../Lab7/data/clinvar_20170207.vcf.gz ./results/hgvs_test_cases_snpEff.vcf > ./results/hgvs_test_cases_snpEff.clinvar.vcf

java -Xmx2G -jar ../Lab7/bin/snpEff/SnpSift.jar extractFields -s ',' -e '.' ./results/hgvs_test_cases_snpEff.clinvar.vcf CHROM POS REF ALT ID "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" CLNHGVS CLNALLE CLNACC CLNSIG CLNREVSTAT CLNDBN > ./results/hgvs_test_cases_snpEff.clinvar.Extracted
```  

### Compare Results

#### VEP  
http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=SW1cOdujfhXtg7pg-2788628  

#### wANNOVAR   
http://wannovar.wglab.org/done/92239/WYSHHkfhaxkURwF3/index.html  

#### SeattleSeq  
http://snp.gs.washington.edu/SeattleSeqAnnotation141/  

### Compare results   
I chose FIG4 for my SNPs, because I like figs with cheese and marmalade.   

According to [GeneCards.org](http://www.genecards.org/cgi-bin/carddisp.pl?gene=FIG4), the protein encoded by this FIG4 "belongs to the SAC domain-containing protein gene family. The SAC domain, approximately 400 amino acids in length and consisting of seven conserved motifs, has been shown to possess phosphoinositide phosphatase activity. The yeast homolog, Sac1p, is involved in the regulation of various phosphoinositides, and affects diverse cellular functions such as actin cytoskeleton organization, Golgi function, and maintenance of vacuole morphology. Membrane-bound phosphoinositides function as signaling molecules and play a key role in vesicle trafficking in eukaryotic cells. Mutations in this gene have been associated with Charcot-Marie-Tooth disease, type 4J." [provided by RefSeq, Jul 2008]

```bash
# method
cut -f9,14-16 hgvs_test_cases_snpEff.clinvar.Extracted > cut.extracted
more cut.extracted | grep "FIG4"
```

#### snpEff + SnpSift Result  
FIG4	c.124_126delGAT	p.Asp42del	.  
FIG4	c.123_124insCAG	p.Ile41_Asp42insGln	.  

#### wANNOVAR Result  
FIG4:NM_014845:exon2:c.123_125del:p.41_42del  
FIG4:NM_014845:exon2:c.123_124insCAG:p.I41delinsIQ  

#### VEP Result  
Not Found (?!)

#### SeattleSeq Result  
Site Crashed

### Sequence Variant Nomenclature Comparison  
There seems to be some consensus between snpEff, SnpSift and wANNOVAR regarding c.123_124insCAG, indicating that there was a CAG (Lipid acid) inserted between nucleotides 123 and 124.  snpEff+SnpSift seems to think that a GAT (Aspartic Acid) nucleotide was removed; while wANNOVAR is unclear what was removed; only that something was removed between nucleotides 123 and 125.  

### Other Links of Interest  
http://varnomen.hgvs.org/    
http://www.genecards.org  
