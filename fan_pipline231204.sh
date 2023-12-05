#!/bin/bash
#conda activate bulk_ATAC
bn=111
#QC before mapping
mkdir 1_QC1
cd ./1_QC1
dir_fastp=/home/fancong/anaconda3/envs/bulk_ATAC/bin/
dir_fastq=/home/yanglab_data3/user/zhaolab_data/huada_ZP-RV-ATAC_1st_118_lack-15-exp/F22FTSSCKF10179_WAIbzcjR/result
function run_fastp(){
    f1=$dir_fastq/${bn}/${bn}_L1_1.fq.gz
    f2=$dir_fastq/${bn}/${bn}_L1_2.fq.gz
    $dir_fastp/fastp -i $f1 -o ./${bn}_fastp_1.fq.gz -I $f2 -O ./${bn}_fastp_2.fq.gz
}
run_fastp $bn
cd ..

#assignment of read to reference genome.
mkdir 2_MAP
cd 2_MAP
fn_ref=/home/fancong/database_fan/ref_genome/bwa_index/mouse/mm10.fa
bwa-mem2 mem -t 8 $fn_ref ../1_QC1/${bn}_fastp_1.fq.gz ../1_QC1/${bn}_fastp_2.fq.gz > ${bn}.sam
cd ..

#QC filter after mapping
mkdir 3_pQC
cd ./3_pQC
dir_run=/home/fancong/anaconda3/envs/samtools/bin
$dir_run/sambamba view  -f bam -h -S -t 8 -F '[XS] <= 20 and not unmapped and not duplicate and mapping_quality >=50' ../2_MAP/${bn}.sam > ${bn}_mapped.bam #for check code
$dir_run/sambamba sort -t 8 ${bn}_mapped.bam -o ${bn}_sorted.bam
$dir_run/sambamba markdup -r -t 8 ${bn}_sorted.bam ${bn}_PCR.bam
$dir_run/samtools index ${bn}_PCR.bam
$dir_run/samtools idxstats ${bn}_PCR.bam > ${bn}_mitoch.stat
$dir_run/samtools view -h ${bn}_PCR.bam |grep -v \'chrM\'|$dir_run/samtools view -bS -o ${bn}_clean.bam
$dir_run/sambamba index -n 8 ${bn}_clean.bam
ATACseqQC
Rscript ./run_ATACseqQC.R $bn
Rscript ./debug_ATACseqQC.R $bn
cd ..

#peak calling
mkdir 4_peak
cd 4_peak
run_HMMRATAC=/home/fancong/anaconda3/envs/bulk_ATAC/bin/HMMRATAC
fn_bam=../3_pQC/splited/${bn}_shift.bam
fn_bam_index=../3_pQC/splited/${bn}_shift.bam.bai
genome_info=/home/yanglab_data3/user/fancong/fan_database/ATACseq_geneinfo/mm10_genometable.txt
black_list=/home/yanglab_data3/user/fancong/fan_database/ENCODE_blacklist/mm10.blacklist.bed
$run_HMMRATAC -b $fn_bam -i $fn_bam_index -g $genome_info -e $black_list
