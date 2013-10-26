#!/bin/bash
#
# initial example of a pipeline script

bwa aln -I -t 8 reference.fa s_1.txt > out.sai 
bwa samse reference.fa out.sai s_1.txt > out.sam 

samtools view -bSu out.sam  | samtools sort -  out.sorted

java -Xmx1g -jar /apps/PICARD/1.95/MarkDuplicates.jar \
                            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
                            METRICS_FILE=out.metrics \
                            REMOVE_DUPLICATES=true \
                            ASSUME_SORTED=true  \
                            VALIDATION_STRINGENCY=LENIENT \
                            INPUT=out.sorted.bam \
                            OUTPUT=out.dedupe.bam 

samtools index out.dedupe.bam 

samtools mpileup -uf reference.fa out.dedupe.bam | /apps/SAMTOOLS/0.1.19/bin/bcftools view -bvcg - > out.bcf

