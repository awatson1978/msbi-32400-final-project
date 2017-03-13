## MSBI 32400 Lab 4
#### Abigail Watson
#### Jan 25th, 2017


```sh
# We’ll use wget to retrieve the human reference proteome from UniProt (~7 MB file)
cd ~/data/ncbi-blast-2.6.0+/
mkdir db
cd db
curl -O http://ftp.gnu.org/gnu/wget/wget-1.16.3.tar.xz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606.fasta.gz
gunzip UP000005640_9606.fasta.gz

# grep the new fasta file for “^>”
# (lines that start with the greater than symbol, meaning they have FASTA annotation)
# pipe that to the word count program wc with the -l (just show the number of lines)
more UP000005640_9606.fasta
grep "^>" UP000005640_9606.fasta           
grep "^>" UP000005640_9606.fasta | wc       
grep "^>" UP000005640_9606.fasta | wc -l       # 21,031 lines

# Tell your server where to find the BLAST binaries
PATH=$PATH:~/data/ncbi-blast-2.6.0+/bin
export PATH

# Need a hidden file for BLAST to work correctly
echo "; Start the selection for BLAST configuration" > ~/.ncbirc
echo "[BLAST]" >> ~/.ncbirc
echo "; Specifies the path where BLAST databases are installed" >> ~/.ncbirc
echo "BLASTDB=/data/ncbi-blast-2.6.0+/db" >> ~/.ncbirc
echo ~/.ncbirc
more ~/.ncbirc


# Make this file a BLAST database
cd ~/data/ncbi-blast-2.6.0+/db
makeblastdb
makeblastdb -in UP000005640_9606.fasta - parse_seqids -dbtype prot
ls

# Let’s download some files to test
# from the protein database
~/edirect/esearch -db protein -query "NP_000509" | ~/edirect/efetch -format fasta > hbb.fasta
~/edirect/esearch -db protein -query "NP_976312" | ~/edirect/efetch -format fasta > myog.fasta
mv *.fasta ~/Code/IntroToBioinformatics/Lab\ 4/data/

# from Canvas
cd ~/Downloads/
mv *.fasta ~/Code/IntroToBioinformatics/Lab\ 4/data/

```

#### Let’s BLAST some proteins

**Run 1 - hbb.fasta**  
```sh
cd ~/Code/IntroToBioinformatics/Lab\ 4/data/
blastp -query hbb.fasta -db ~/data/ncbi-blast-2.6.0+/db/UP000005640_9606.fasta
blastp -query hbb.fasta -db ~/data/ncbi-blast-2.6.0+/db/UP000005640_9606.fasta -out mysearch.html -html
```
P68871  Hemoglobin subunit beta OS=Homo sapiens GN=HBB PE=1 SV=2      301     **4e-107**


**Run 2 - unknown.fasta**  
```sh
blastp -query unknown.fasta -db ~/data/ncbi-blast-2.6.0+/db/UP000005640_9606.fasta
blastp -query unknown.fasta -db ~/data/ncbi-blast-2.6.0+/db/UP000005640_9606.fasta -out mysearch2.html -html
```

P69905  Hemoglobin subunit alpha OS=Homo sapiens GN=HBA1 PE=1 SV=2    286     **1e-101**  

**Run 3 - unknown2.fasta**
```sh
blastp -query unknown2.fasta -db ~/data/ncbi-blast-2.6.0+/db/UP000005640_9606.fasta
blastp -query unknown2.fasta -db ~/data/ncbi-blast-2.6.0+/db/UP000005640_9606.fasta -out mysearch3.html -html
```

P68871  Hemoglobin subunit beta OS=Homo sapiens GN=HBB PE=1 SV=2      298     **3e-106**  

MVHLTPVEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK    
Query had a V at position 6   
Reference sequence has a E at position 6  

**Run 4 - unknown_fragment.fasta**  
BRCT_assoc	Serine-rich domain associated with BRCT; This domain is found on BRCA1 proteins.  
E-Value:  	1.07e-96  
