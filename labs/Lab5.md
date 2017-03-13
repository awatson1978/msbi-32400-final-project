Abigail Watson  
MSBI 32400 - Lab 5  


```sh
# set up the lab directory structure  
cd Lab5
mkdir src
mkdir data
mkdir results
mkdir source
mkdir docs
touch docs/README_rawatson.md

# take a peek at the 23andme file
less 5657.23andme.4141
# rsid  
# chromosome      
# position        
# genotype

# oops; we don't have perl.  download and install it
curl -L http://xrl.us/installperlosx | bash

# lets convert our 23andme file into a vcf file
cd 23andme2vcf/
perl 23andme2vcf.pl ../../data/5657.23andme.4141 ../../results/5657.23andme.4141.vcf 4 ;

# 8030 sites were not included; these unmatched references can be found in sites_not_in_reference.txt.Try running again, but specify the other reference version:
# ./23andme2vcf.pl ../../data/5657.23andme.4141 ../../results/5657.23andme.4141.vcf 3
```


#### How many lines were excluded?  
8030  

#### How many lines are there in the VCF (not counting the header)?  
577020  


#### List a few stop codons

```
chr1	12776218	rs3000860	A	C	.	.	ANN=C|stop_lost|HIGH|AADACL3|AADACL3|transcript|NM_001103170.2|protein_coding|1/4|c.39A>C|p.Ter13Cysext*?|101/4049|39/1224|13/407||	GT	0/1  

chr1	21053559	i6019527	G	A	.	.	ANN=A|stop_gained|HIGH|SH2D5|SH2D5|transcript|NM_001103161.1|protein_coding|4/10|c.178C>T|p.Arg60*|680/3834|178/1272|60/423||;LOF=(SH2D5|SH2D5|1|1.00);NMD=(SH2D5|SH2D5|1|1.00)	GT	1/1  

chr1	47080679	rs6671527	G	A	.	.	ANN=A|stop_gained|HIGH|MOB3C|MOB3C|transcript|NM_145279.4|protein_coding|1/4|c.70C>T|p.Arg24*|127/2804|70/807|24/268||;LOF=(MOB3C|MOB3C|1|1.00);NMD=(MOB3C|MOB3C|1|1.00)	GT	0/1  
```

#### How many were classified as “stop_lost” and how many as “stop_gained”?  

*stop_gained*	 	
47  
0.015%   

*stop_lost*  	 	  
15	 
0.005%    


#### Annotate snpEff.vcf with SnpSift
```
java -Xmx2G -jar ./SnpSift.jar annotate -noLog ../../data/clinvar_20170104.vcf.gz ../../results/5657.23andme.4141.snpEff.vcf > ../../results/5657.23andme.4141.clinvar.snpEff.vcf  
```

```
chr1	98348885	rs1801265	G	A	.	.	ANN=A|missense_variant|MODERATE|DPYD|DPYD|transcript|NM_000110.3|protein_coding|2/23|c.85C>T|p.Arg29Cys|222/4447|85/3078|29/1025||;ASP;CAF=0.2602,0.7398;CLNACC=RCV000000464.2,RCV000086506.1;CLNALLE=0,1;CLNDBN=Dihydropyrimidine_dehydrogenase_deficiency,not_provided;CLNDSDB=MedGen:OMIM,MedGen;CLNDSDBID=C2720286:274270,CN221809;CLNHGVS=NC_000001.10:g.98348885G\x3d,NC_000001.10:g.98348885G>A;CLNORIGIN=1,0;CLNREVSTAT=no_criteria,no_assertion;CLNSIG=5,1;CLNSRC=OMIM_Allelic_Variant|UniProtKB_(protein),.;CLNSRCID=612779.0004|Q12882#VAR_005173,.;COMMON=1;G5;GENEINFO=DPYD:1806;GNO;HD;INT;KGPhase1;KGPhase3;LSD;NSM;OM;PM;PMC;REF;RS=1801265;RSPOS=98348885;RV;SAO=1;SLO;SSR=0;TPA;VC=SNV;VLD;VP=0x050178080a0515053f110100;WGT=1;dbSNPBuildID=89	GT	1/1  

chr1	100672060	rs12021720	T	C	.	.	ANN=C|missense_variant|MODERATE|DBT|DBT|transcript|NM_001918.3|protein_coding|9/11|c.1150A>G|p.Ser384Gly|1183/10815|1150/1449|384/482||;ASP;CAF=0.1082,0.8918;CLNACC=RCV000012727.23,RCV000116865.5;CLNALLE=0,1;CLNDBN=Intermediate_maple_syrup_urine_disease_type_2,not_specified;CLNDSDB=MedGen,MedGen;CLNDSDBID=CN069615,CN169374;CLNHGVS=NC_000001.10:g.100672060T\x3d,NC_000001.10:g.100672060T>C;CLNORIGIN=1,1;CLNREVSTAT=no_criteria,no_criteria;CLNSIG=5,255;CLNSRC=OMIM_Allelic_Variant,.;CLNSRCID=248610.0008,.;COMMON=1;G5;GENEINFO=DBT:1629;GNO;HD;KGPhase1;KGPhase3;LSD;NSM;OM;PM;PMC;REF;RS=12021720;RSPOS=100672060;SAO=1;SLO;SSR=0;VC=SNV;VLD;VP=0x050168000a0515053f110100;WGT=1;dbSNPBuildID=120	GT	1/1  

chr1	115236057	rs17602729	G	A	.	.	ANN=A|stop_gained&splice_region_variant|HIGH|AMPD1|AMPD1|transcript|NM_000036.2|protein_coding|2/16|c.133C>T|p.Gln45*|181/2406|133/2343|45/780||;LOF=(AMPD1|AMPD1|1|1.00);NMD=(AMPD1|AMPD1|1|1.00);ASP;CAF=0.9619,0.03814,.;CLNACC=RCV000019933.29;CLNALLE=1;CLNDBN=Muscle_AMP_deaminase_deficiency;CLNDSDB=MedGen:OMIM:SNOMED_CT;CLNDSDBID=C0268123:615511:9105005;CLNHGVS=NC_000001.10:g.115236057G>A;CLNORIGIN=1;CLNREVSTAT=no_criteria;CLNSIG=5;CLNSRC=OMIM_Allelic_Variant;CLNSRCID=102770.0001;COMMON=1;G5;GENEINFO=AMPD1:270;GNO;HD;INT;KGPhase1;KGPhase3;LSD;NSN;OM;PM;PMC;REF;RS=17602729;RSPOS=115236057;RV;SAO=1;SLO;SSR=0;VC=SNV;VLD;VP=0x05016808060515053f110100;WGT=1;dbSNPBuildID=123	GT	0/1  
```

#### Pharmacogenomics  

```sh
grep 'rs1800462' 5657.23andme.4141  # CC  
grep 'rs1142345' 5657.23andme.4141  # TT  
grep 'rs1800460' 5657.23andme.4141  # CC  
grep 'rs1800584' 5657.23andme.4141  # CC  
```

Star allelle *1   

#### Pharmacogenomics - TPMT Gene  

Patients with the CC genotype may have a decreased risk for toxicity with thiopurine drugs and purine analogues as compared to patients with the CT or TT genotype. Patients with the CC genotype may still be at risk for toxicity when taking thiopurine drugs and purine analogues based on their genotype. Other genetic and clinical factors may also influence a patient's risk for toxicity.

Patients with the CC genotype: 1) may have decreased metabolism of thiopurines 2) may have a decreased risk for toxicity with thiopurine drugs as compared to patients with the CG or GG genotype. Patients with the CC genotype may still be at risk for toxicity when taking thiopurine drugs based on their genotype. Other genetic and clinical factors may also influence a patient's risk for toxicity.
