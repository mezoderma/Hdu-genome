#!/bin/bash 

## kraken2 

TAG=illumina
R1=$TAG.R1.fastq.gz
R2=$TAG.R2.fastq.gz
DB1=k2_pluspf_20210517
DB2=kraken2_nt_20200127

kraken2 --threads 64 --gzip-compressed --db $DB1 --output $TAG.kraken2_pluspf.out --report $TAG.kraken2_pluspf.report --paired $R1 $R2
kraken2 --threads 64 --gzip-compressed --db $DB2 --output $TAG.kraken2_nt.out --report $TAG.kraken2_nt.report --paired $R1 $R2

## bbduk v38.07

bbduk.sh in1=${TAG}.R1.fastq.gz in2=${TAG}.R2.fastq.gz out1=$TAG.bbduk.R1.fastq out2=$TAG.bbduk.R2.fastq ref=/bbmap/resources/adapters.fa ktrim=rl k=23 mink=11 hdist=1 tpe tbo &> $TAG.bbduk.log 

## filtlong v0.2.0

filtlong -1 illumina.bbduk.R1.fastq.gz -2 illumina.bbduk.R2.fastq.gz --min_length 1000 --keep_percent 90 --target_bases 500000000 --trim --split 500 nanopore.fastq.gz > Fl90_pass.fastq

## flye v2.9-b1768

FASTQ=Fl90_pass.fastq
flye --nano-raw $FASTQ --out-dir new_flye.out --threads 64 &> new_flye.log
flye --nano-raw $FASTQ --out-dir new_flye_8k.out --threads 64 --min-overlap 8000 &> new_flye_8k.log
flye --nano-raw $FASTQ --out-dir new_flye_10k.out --threads 64 --min-overlap 10000 &> new_flye_10k.log

## canu v2.1.1

FASTQ=Fl90_pass.fastq
canu -d hadu.canu.out -p hadu genomeSize=250m maxMemory=400G -nanopore $FASTQ &> hadu.canu.log 

## shasta v0.7.0

FASTQ=Fl90_pass.fastq
CFG=shasta/conf/Nanopore-Sep2020.conf
GPY=shasta/conf/SimpleBayesianConsensusCaller-10.csv
shasta-Linux-0.7.0 --threads 64 --input $FASTQ --config $CFG --Reads.minReadLength 2000 --Assembly.consensusCaller Bayesian:$GPY &> shasta.log # resulted in assembly_shasta.fa.gz
shasta-Linux-0.7.0 --threads 64 --input $FASTQ --config $CFG --Reads.minReadLength 2000 --Assembly.consensusCaller Bayesian:$GPY --assemblyDirectory Shasta_sep2020_rl2k.out &> Shasta_sep2020_rl2k.log # resulted in assembly_shasta_rl2k.fa.gz
shasta-Linux-0.7.0 --threads 64 --input $FASTQ --config $CFG --Reads.minReadLength 5000 --Assembly.consensusCaller Bayesian:$GPY --assemblyDirectory Shasta_sep2020_rl5k.out &> Shasta_sep2020_rl5k.log # resulted in assembly_shasta_r5k.fa.gz, choosen for polishing
shasta-Linux-0.7.0 --threads 64 --input $FASTQ --config $CFG --Reads.minReadLength 8000 --Assembly.consensusCaller Bayesian:$GPY --assemblyDirectory Shasta_sep2020_rl8k.out &> Shasta_sep2020_rl8k.log # resulted in assembly_shasta_r8k.fa.gz 

## raven v1.5.1

FASTQ=Fl90_pass.fastq
raven -t32 $FASTQ &> raven.log

## wtdbg2 v2.5

FASTQ=Fl90_pass.fastq
wtdbg2.pl -t32 -x ont -g 250m -o Fl_90_pass $FASTQ &> wtdbg2.log 

## racon

cat illumina.bbduk.R1.fastq illumina.bbduk.R2.fastq > illumina.combined.fastq
FA=assembly_shasta_rl5k.fa.gz
ILL=illumina.combined.fastq
FNAME=`basename $FA`
TAG=${FNAME%%.fa.gz}
cp $FA $TAG.ill.0.fa
for i in `seq 1 5`
do
  j=$((i-1))
  echo "["`date +%H:%M:%S`"] Iteration $i: Mapping reads to fasta file $TAG.ill.$j.fa"
  minimap2 -x sr -t64 -a $TAG.ill.$j.fa $ILL > $TAG.mini.$i.sam 2> $TAG.mini.$i.log
  echo "["`date +%H:%M:%S`"] Iteration $i: Racon correction to generate new consensus $TAG.ill.$i.fa"
  racon -w 100 -t64 $ILL $TAG.mini.$i.sam $TAG.ill.$j.fa > $TAG.ill.$i.fa 2> $TAG.racon.$i.log
  rm $TAG.mini.$i.sam
done

## polca

FA=assembly_shasta_rl5k.fa.gz
TAG=${FA%%.fa.gz}
FA=`readlink -f $FA`
R1=illumina.bbduk.R1.fastq.gz
R2=illumina.bbduk.R2.fastq.gz
mkdir $TAG.polca.out
cd $TAG.polca.out
polca.sh -t 16 -a $FA -r "$R1 $R2" &> $TAG.polca.log

## pilon

FA=assembly_shasta_rl5k.fa.gz
TAG=${FA%%.fa.gz}
FA=`readlink -f $FA`
R1=illumina.bbduk.R1.fastq.gz
R2=illumina.bbduk.R2.fastq.gz
bowtie2-build $FA $TAG
bowtie2 -x $TAG -1 $R1 -2 $R2 -S $TAG.sam
samtools sort -@ 20 -o $TAG.bam $TAG.sam
samtools view -@ 20 -bS -o $TAG.sorted.bam $TAG.bam
samtools index $TAG.sorted.bam
java -Xmx64G -jar pilon.jar

## busco v5.2.2
pilon --genome $FA --frags $TAG.sorted.bam



