# Usage

## 1.Download genome reference sequences and refGene.txt

First download the genome reference sequences and refGene.txt table, which will be used in the following script to get the promoter seuqnces based the gene annotation. The file will be store in the **data** folder, and the reference sequences and refGene.txt will be named as **hg19_genome.fa** and **hg19_refGene.txt** individually.

```
downloadData.pl -build hg19
```

## 2. Prepare the promoter sequences

Once we downloaded the genome reference and refGene.txt, we will prepare the promoter sequences in fasta format. For genes having multiple isoforms on the same chromosome and same strand, we will apply the leftmost position on the 5' as their TSS. To get the promoters, simply typing the following,(the result file is named as **hg19_promoter_upstream2000_downstream500.fa**)

```
preparePromoter.pl -build hg19 -up 2000 -down 500
``` 
