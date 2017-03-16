#### Setting up my project

- cd /data
- mkdir myproject
- mkdir myproject/data
- mkdir myproject/doc
- mkdir myproject/results
- mkdir myproject/src
- mkdir myproject/bin
- touch myproject/doc/README_rawatson.md
- tree myproject


```bash
# Installing NCBI Command Line Tools  
# Try a few of the examples in the book  
~/edirect/esearch –db pubmed –query “pevsner j AND gnaq” | ~/edirect/efetch –format pubmed > doc/ example1.txt  
~/edirect/esearch -db pubmed -query "bioinformatics [MAJR] AND software [TIAB]" | ~/edirect/efetch -format xml | xtract -pattern PubmedArticle -block Author -sep " " -tab "\n" -element LastName,Initials | sort-uniq-count-rank > doc/bioinformatics_authors.txt  
~/edirect/esearch -db protein -query 'NP_000509.1' | ~/ edirect/efetch -format fasta > doc/hbb.fasta  
```
