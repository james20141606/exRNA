#!/bin/bash

path0=/BioII/lulab_b/younglee/exRNA/hcc_lulab

for i in rRNA miRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $i;
gffread -w $path0/src/$i.fa -g /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p10.genome.fa /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/$i.gencode27.gtf
rsem-prepare-reference --gtf /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/$i.gencode27.gtf --bowtie2 /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p10.genome.fa $path0/src/${i}
done;

gffread -w $path0/src/piRNA.fa -g /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p10.genome.fa /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/piRNA.piRBase.hg38.gtf
rsem-prepare-reference --gtf /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/piRNA.piRBase.hg38.gtf --bowtie2 /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p10.genome.fa $path0/src/piRNA

for i in `ls 02.mapping/`;
do echo $i;
cleanN=`cat 01.preprocess/$i.cutadapt.fastq | wc -l`
collapseN=`cat 01.preprocess/$i.cutadapt.collapsed.fastq | wc -l`
rRNAN=`samtools view 02.mapping/$i/remove_rRNA/$i.rRNA.clean.bam | wc -l`
keptN=`cat 02.mapping/$i/remove_rRNA/$i.rRNA.unmapped.fastq | wc -l`
map2hg38N=`samtools view 02.mapping/$i/hg38/$i.hg38.bam | wc -l`
dedupN=`samtools view 02.mapping/$i/hg38/$i.hg38.dedup.bam | wc -l`
markDuplicatesN=`samtools view 02.mapping/$i/hg38/$i.hg38.markDuplicates.bam | wc -l`
echo -e "$i\t$cleanN\t$collapseN\t$rRNAN\t$keptN\t$map2hg38N\t$dedupN\t$markDuplicatesN" | awk 'BEGIN{FS=OFS="\t"}{print $1,$2/4,$3/2,$4,$5/4,$6,$7,$8}' >> hcc_lulab.readsN.stat.txt
done;

sed -i "1i sample\tcleanN\tcollapseN\trRNAN\tkeptN\tmap2hg38N\tdedupN\tmarkDuplicatesN" hcc_lulab.readsN.stat.txt


filename=`grep -v ">" 385247-A.cutadapt.collapsed.fastq`
filename="foo"
while read -r line
do
    name="$line"
    num=`grep -w $line 385247-A.cutadapt.fastq | wc -l`
    if [$num >= 2]
    then 
        grep -w -A 2 -B 1 "$name" 385247-A.cutadapt.fastq > foo2
    fi
#    echo "Name read from file - $name"
done < "$filename"

for i in `cat | grep -v ">"`; do echo $i >> foo; done


## merge annotations
# cat piRNA.piRBase.hg38.gtf Y_RNA.gencode27.gtf snRNA.gencode27.gtf snRNA.gencode27.gtf srpRNA.gencode27.gtf tRNA.gencode27.gtf > canonical_ncRNA.gtf
# cat piRNA.piRBase.hg38.gff Y_RNA.gencode27.gff snRNA.gencode27.gff snRNA.gencode27.gff srpRNA.gencode27.gff tRNA.gencode27.gff > canonical_ncRNA.gff

gffread -w $path0/src/canonical_ncRNA.fa -g /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p10.genome.fa /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/piRNA.piRBase.hg38.gtf
rsem-prepare-reference --gtf /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/canonical_ncRNA.gtf --bowtie2 /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p10.genome.fa $path0/src/canonical_ncRNA


##----------------------------------------
## summary the results and statistics
path0=/Share/home/younglee/projects/exRNA/hcc_lulab
mkdir $path0/stat
## number of mapped readsfor each RNA types
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# total library size
libSizeN=`echo $(cat $path0/00.rawdata/$j.fastq | wc -l)/4 | bc`
echo -e "$j\tpreprocess\tlibSizeN\t$libSizeN" >> $path0/stat/$j.readsN.stat.tsv
# too short
tooShortN=`echo $(cat $path0/01.preprocess/$j.tooShort.fastq | wc -l)/4 | bc`
echo -e "$j\tpreprocess\ttooShortN\t$tooShortN" >> $path0/stat/$j.readsN.stat.tsv
# clean reads
cleanN=`echo $(cat $path0/01.preprocess/$j.cutadapt.fastq | wc -l)/4 | bc`
echo -e "$j\tpreprocess\tcleanN\t$cleanN" >> $path0/stat/$j.readsN.stat.tsv
# rRNA mapped reads
rRNA_N=`samtools flagstat $path0/02.mapping/$j/remove_rRNA/$j.rRNA.clean.bam | awk 'NR==5{print $1}'`
echo -e "$j\tpreprocess\trRNA_N\t$rRNA_N" >> $path0/stat/$j.readsN.stat.tsv
# kept reads
keptN=`echo $(cat $path0/02.mapping/$j/remove_rRNA/$j.rRNA.unmapped.fastq | wc -l)/4 | bc`
echo -e "$j\tpreprocess\tkeptN\t$keptN" >> $path0/stat/$j.readsN.stat.tsv
# map to hg38
hg38_N=`samtools flagstat $path0/02.mapping/$j/hg38/$j.hg38.bam | awk 'NR==5{print $1}'`
echo -e "$j\tmap2hg38\thg38\t$hg38_N" >> $path0/stat/$j.readsN.stat.tsv
# map to different RNA types
for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
RNAtypes_N=`samtools flagstat $path0/02.mapping/$j/sequentialMap/$j.$k.rsem.clean.bam | awk 'NR==5{print $1}'`
echo -e "$j\tsequentialMap\t$k\t$RNAtypes_N" >> $path0/stat/$j.readsN.stat.tsv
done;
# map to hg38 other region
hg38other_N=`samtools flagstat $path0/02.mapping/$j/sequentialMap/$j.hg38other.bam | awk 'NR==5{print $1}'`
echo -e "$j\tmap2hg38other\thg38other\t$hg38other_N" >> $path0/stat/$j.readsN.stat.tsv
# non-human
nonHuman_N=`echo $(cat $path0/02.mapping/$j/sequentialMap/$j.hg38other.unmapped.fastq | wc -l)/4 | bc`
echo -e "$j\tmap2hg38other\tnonHuman_N\t$nonHuman_N" >> $path0/stat/$j.readsN.stat.tsv
done;


cut -f 2,3 $path0/stat/LY.readsN.stat.tsv | sed '1i sample' > $path0/stat/readsN.stat.header
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
cut -f 4 $path0/stat/$j.readsN.stat.tsv | sed -e "1i $j" > $path0/stat/$j.readsN.stat.foo
done;
paste $path0/stat/readsN.stat.header $path0/stat/*.readsN.stat.foo > $path0/stat/hcc_lulab.readsN.stat.tsv
rm -rf $path0/stat/*.readsN.stat.foo





## length of mapped reads for each RNA types
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
for k in miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA;
do echo $k;
samtools view $path0/02.mapping/$j/sequentialMap/$j.$k.rsem.clean.bam | awk 'BEGIN{FS=OFS="\t"}{print length($10)}' | sort -n | uniq -c | awk 'BEGIN{FS=" "; OFS="\t"}{print $2,$1}' | sort -nk1,1 | sed -e "s/^/$j\t$k\t/g" >> $path0/stat/$j.lengthN.stat.tsv
done;
sed -i -e "1i sample\ttype\tlen\tnum" $path0/stat/$j.lengthN.stat.tsv
done;

# plot the histogram
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
Rscript $path0/bin/plot_readsLens.R $path0/stat/$j.lengthN.stat.tsv $path0/stat/$j.lengthN.stat.pdf
done;



