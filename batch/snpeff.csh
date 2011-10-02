#!/bin/csh
#PBS -q small 
#PBS -N sift
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -o /home/aheiberg/andrew/batch/output/snpeff_out.txt
#PBS -e /home/aheiberg/andrew/batch/output/snpeff_err.txt
#PBS -V
#PBS -M andrew.heiberg@gmail.com
#PBS -m abe
#PBS -A gleeson-lab

cd /projects/gleeson-lab/bin/snpEff_v1_9_5

java -Xmx20g -jar snpEff.jar -f hg37 /projects/gleeson-lab/BroadData/PlateII/raw_data/SNPS_SUBSET.vcf

mv /home/aheiberg/andrew/batch/output/snpeff_out.txt /projects/gleeson-lab/BroadData/PlateII/intermediate_data/snpeff/output/SNPS_SUBSET.csv

