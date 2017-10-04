# imrep-10x
10x


# One sample 
```
python ../barcode.py 03T3.txt 03T3_S48_L001_R2_001.fastq full_cdr3_03T3_S48_L001_R1_001.fastq.txt  test
python ../perReads2perbarcode.py test_raw.txt test2
``` 

# Multiple samples

```
while read line; do sample=$(echo $line | awk -F "_" '{print $1}'); echo $sample;python ~/code/imrep-10x/imrep-10x_step1.py barcodes/${sample}.txt R2/${line}_R2_001.fastq cdr3_per_read/full_cdr3_${line}_R1_001.fastq.txt  analyis_per_barcode/${line}.csv ;done<samples.txt
while read line; do cat ../raw/${line}* >${line}_perRead.csv;done<../samples.txt
while read line; do python ../../imrep-10x_step2.py ${line}_perRead.csv ${line}_perB.csv;done<../samples.txt
```


