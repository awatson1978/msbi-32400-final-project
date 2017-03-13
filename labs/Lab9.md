Abigail Watson  
MSB32400 - Lab #9  


## Bash exercises  
```
mkdir Lab9  
cd Lab9  
mkdir bin  
mkdir docs  
mkdir results  
mkdir src  

cd /data/bds-files/chapter-12-pipelines  
cp template.sh /data/Lab9/bin  
cd /data/Lab9/bin  

# Use find + xargs to find the WebDocument files and show their size (slide #12)  
find /data -name "WebDocument*" -print | xargs ls -lh
# -rw-rw-r--. 1 student student 3.1K Jan 18 19:18 /data/Lab3/data/WebDocument_9-05_101AutismPanel.bed.txt
# -rw-rw-r--. 1 student student 3.0M Jan 18 19:18 /data/Lab3/data/WebDocument_9-7_mysample1.bai
# -rw-rw-r--. 1 student student 326M Jan 18 19:55 /data/Lab3/data/WebDocument_9-7_mysample1.bam
# -rw-rw-r--. 1 student student 1.1G Jan 18 20:40 /data/Lab3/data/WebDocument_9-7_mysample1.fastq
# -rw-rw-r--. 1 student student 327M Jan 23 12:15 /data/Lab3/data/WebDocument_9-7_mysample1_file_sorted.bam
# -rw-rw-r--. 1 student student 3.0M Jan 23 12:18 /data/Lab3/data/WebDocument_9-7_mysample1_file_sorted.bam.bai
# -rw-rw-r--. 1 student student 1.6G Jan 18 20:40 /data/Lab3/data/WebDocument_9-7_mysample1.sam


# Use find + xargs to create a bash file to find all .fastq files and delete them   
echo "find /data -name "*.fastq" -print | xargs rm" >> /data/Lab9/bin/template.sh  

# Find all data directories beneath /data  
find /data -name "*data*" -print  

# /data
# /data/Lab8/data
# /data/Lab3/data
# /data/snpEff/data
# /data/snpEff/galaxy/tool-data
# /data/bds-files/chapter-09-working-with-range-data
# /data/bds-files/chapter-06-bioinformatics-data
# /data/bds-files/chapter-08-r/motif-example/data
# /data/bds-files/chapter-07-unix-data-tools
# /data/Lab2/data
# /data/Lab9/data
# /data/lab6/data


# Use find to locate files modified within the past 14 days then use xargs to find their size
find /data -ctime 14 -print | xargs ls -lh
```


