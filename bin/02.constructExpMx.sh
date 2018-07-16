#!/bin/bash

path0=/Share/home/younglee/projects/exRNA/hcc_lulab

## gene ##
# merge expression matrix for each RNA type
for i in `ls $path0/04.counts`;
do j=${i%.*};
echo $j;
for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
sed '1d' $path0/04.counts/${j}/sequentialMap/${j}.${k}.homer.ct | cut -f 9 | sed -e "1i lulab_${j}" | sed 's/-/_/g' > $path0/tmp/${j}.${k}.homer.ct.tmp
head -n -5 $path0/04.counts/${j}/sequentialMap/${j}.${k}.htseq.ct | cut -f 2 | sed -e "1i lulab_${j}" | sed 's/-/_/g' > $path0/tmp/${j}.${k}.htseq.ct.tmp
sed '1,2d' $path0/04.counts/${j}/sequentialMap/${j}.${k}.featureCounts.ct | cut -f 7 | sed -e "1i lulab_${j}" | sed 's/-/_/g' > $path0/tmp/${j}.${k}.featureCounts.ct.tmp
done;
done;

for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
sed '1d' $path0/04.counts/Normal-1/sequentialMap/Normal-1.${k}.homer.ct | cut -f 1 | sed -e "s/^/${k}_/g" | sed -e "1i geneID" > $path0/tmp/$k.homer.header
head -n -5 $path0/04.counts/Normal-1/sequentialMap/Normal-1.${k}.htseq.ct | cut -f 1 | sed -e "s/^/${k}_/g" | sed -e "1i geneID" > $path0/tmp/$k.htseq.header
sed '1,2d' $path0/04.counts/Normal-1/sequentialMap/Normal-1.${k}.featureCounts.ct | cut -f 1 | sed -e "s/^/${k}_/g" | sed -e "1i geneID" > $path0/tmp/$k.featureCounts.header
done;

for i in homer htseq featureCounts;
do echo $i;
for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
paste $path0/tmp/$k.$i.header $path0/tmp/*.${k}.$i.ct.tmp > $path0/05.matrix/hcc_lulab.sequentialMap.$i.$k.mx
sed -i "1s/......\t//1" $path0/05.matrix/hcc_lulab.sequentialMap.$i.$k.mx
done;
done;


## merge all RNA-type
for i in `ls $path0/04.counts`;
do j=${i%.*};
echo $j;
for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
sed '1d' $path0/04.counts/${j}/sequentialMap/${j}.${k}.homer.ct | cut -f 9 > $path0/tmp/${j}.${k}.homer.ct.tmp
head -n -5 $path0/04.counts/${j}/sequentialMap/${j}.${k}.htseq.ct | cut -f 2 > $path0/tmp/${j}.${k}.htseq.ct.tmp
sed '1,2d' $path0/04.counts/${j}/sequentialMap/${j}.${k}.featureCounts.ct | cut -f 7 > $path0/tmp/${j}.${k}.featureCounts.ct.tmp
done;
for m in homer htseq featureCounts;
do echo $m;
cat $path0/tmp/${j}.miRNA.${m}.ct.tmp $path0/tmp/${j}.piRNA.${m}.ct.tmp $path0/tmp/${j}.Y_RNA.${m}.ct.tmp $path0/tmp/${j}.snRNA.${m}.ct.tmp $path0/tmp/${j}.snoRNA.${m}.ct.tmp $path0/tmp/${j}.srpRNA.${m}.ct.tmp $path0/tmp/${j}.tRNA.${m}.ct.tmp $path0/tmp/${j}.lncRNA.${m}.ct.tmp $path0/tmp/${j}.mRNA.${m}.ct.tmp | sed -e "1i lulab_${j}" | sed 's/-/_/g' > $path0/tmp/${j}.merged.${m}.ct.tmp
cat $path0/tmp/miRNA.${m}.header $path0/tmp/piRNA.${m}.header $path0/tmp/Y_RNA.${m}.header $path0/tmp/snRNA.${m}.header $path0/tmp/snoRNA.${m}.header $path0/tmp/srpRNA.${m}.header $path0/tmp/tRNA.${m}.header $path0/tmp/lncRNA.${m}.header $path0/tmp/mRNA.${m}.header | grep -v "geneID" | sed -e "1i geneID" > $path0/tmp/merged.${m}.header
done;
done;

for m in homer htseq featureCounts;
do echo $m;
paste $path0/tmp/merged.${m}.header $path0/tmp/*.merged.${m}.ct.tmp > $path0/05.matrix/hcc_lulab.sequentialMap.${m}.merged.mx
sed -i "1s/......\t//1" $path0/05.matrix/hcc_lulab.sequentialMap.${m}.merged.mx
done;


## binned ##
# merge expression matrix for each RNA type
for i in `ls $path0/04.counts`;
do j=${i%.*};
echo $j;
for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
sed '1d' $path0/04.counts/${j}/sequentialMap/${j}.${k}.binned.homer.ct | cut -f 9 | sed -e "1i lulab_${j}" | sed 's/-/_/g' > $path0/tmp/${j}.${k}.binned.homer.ct.tmp
head -n -5 $path0/04.counts/${j}/sequentialMap/${j}.${k}.binned.htseq.ct | cut -f 2 | sed -e "1i lulab_${j}" | sed 's/-/_/g' > $path0/tmp/${j}.${k}.binned.htseq.ct.tmp
sed '1,2d' $path0/04.counts/${j}/sequentialMap/${j}.${k}.binned.featureCounts.ct | cut -f 7 | sed -e "1i lulab_${j}" | sed 's/-/_/g' > $path0/tmp/${j}.${k}.binned.featureCounts.ct.tmp
done;
done;

for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
sed '1d' $path0/04.counts/Normal-1/sequentialMap/Normal-1.${k}.binned.homer.ct | cut -f 1 | sed -e "s/^/${k}_/g" | sed -e "1i geneID" > $path0/tmp/$k.binned.homer.header
head -n -5 $path0/04.counts/Normal-1/sequentialMap/Normal-1.${k}.binned.htseq.ct | cut -f 1 | sed -e "s/^/${k}_/g" | sed -e "1i geneID" > $path0/tmp/$k.binned.htseq.header
sed '1,2d' $path0/04.counts/Normal-1/sequentialMap/Normal-1.${k}.binned.featureCounts.ct | cut -f 1 | sed -e "s/^/${k}_/g" | sed -e "1i geneID" > $path0/tmp/$k.binned.featureCounts.header
done;

for i in homer;
#htseq featureCounts;
do echo $i;
for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA;
do echo $k;
paste $path0/tmp/$k.binned.$i.header $path0/tmp/*.${k}.binned.$i.ct.tmp > $path0/05.matrix/hcc_lulab.sequentialMap.$i.$k.binned.mx
sed -i "1s/......\t//1" $path0/05.matrix/hcc_lulab.sequentialMap.$i.$k.binned.mx
done;
done;


for i in htseq featureCounts;
do echo $i;
for k in mRNA;
do echo $k;
paste $path0/tmp/$k.binned.$i.header $path0/tmp/*.${k}.binned.$i.ct.tmp > $path0/05.matrix/hcc_lulab.sequentialMap.$i.$k.binned.mx
sed -i "1s/......\t//1" $path0/05.matrix/hcc_lulab.sequentialMap.$i.$k.binned.mx
done;
done;



## merge all RNA-type
for i in `ls $path0/04.counts`;
do j=${i%.*};
echo $j;
for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
sed '1d' $path0/04.counts/${j}/sequentialMap/${j}.${k}.binned.homer.ct | cut -f 9 > $path0/tmp/${j}.${k}.binned.homer.ct.tmp
head -n -5 $path0/04.counts/${j}/sequentialMap/${j}.${k}.binned.htseq.ct | cut -f 2 > $path0/tmp/${j}.${k}.binned.htseq.ct.tmp
sed '1,2d' $path0/04.counts/${j}/sequentialMap/${j}.${k}.binned.featureCounts.ct | cut -f 7 > $path0/tmp/${j}.${k}.binned.featureCounts.ct.tmp
done;


for m in homer;
# htseq featureCounts;
do echo $m;
cat $path0/tmp/${j}.miRNA.binned.${m}.ct.tmp $path0/tmp/${j}.piRNA.binned.${m}.ct.tmp $path0/tmp/${j}.Y_RNA.binned.${m}.ct.tmp $path0/tmp/${j}.snRNA.binned.${m}.ct.tmp $path0/tmp/${j}.snoRNA.binned.${m}.ct.tmp $path0/tmp/${j}.srpRNA.binned.${m}.ct.tmp $path0/tmp/${j}.tRNA.binned.${m}.ct.tmp $path0/tmp/${j}.lncRNA.binned.${m}.ct.tmp $path0/tmp/${j}.mRNA.binned.htseq.ct.tmp | sed -e "1i lulab_${j}" | sed 's/-/_/g' > $path0/tmp/${j}.merged.binned.${m}.ct.tmp
cat $path0/tmp/miRNA.binned.${m}.header $path0/tmp/piRNA.binned.${m}.header $path0/tmp/Y_RNA.binned.${m}.header $path0/tmp/snRNA.binned.${m}.header $path0/tmp/snoRNA.binned.${m}.header $path0/tmp/srpRNA.binned.${m}.header $path0/tmp/tRNA.binned.${m}.header $path0/tmp/lncRNA.binned.${m}.header $path0/tmp/mRNA.binned.htseq.header | grep -v "geneID" | sed -e "1i geneID" > $path0/tmp/merged.binned.${m}.header
done;
done;

for m in homer;
#htseq featureCounts;
do echo $m;
paste $path0/tmp/merged.binned.${m}.header $path0/tmp/*.merged.binned.${m}.ct.tmp > $path0/05.matrix/hcc_lulab.sequentialMap.${m}.merged.binned.mx
sed -i "1s/......\t//1" $path0/05.matrix/hcc_lulab.sequentialMap.${m}.merged.binned.mx
done;



