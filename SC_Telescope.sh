#! /bin/bash -l


cd /athena/nixonlab/scratch/lpr4001/Howard_Fine/
NAME=$(head -n  ${SLURM_ARRAY_TASK_ID} SC_samples.txt| tail -n 1)
# head -n 3 SC_samples.txt
# Sample_PM1005_2D.tar
# Sample_PM1005_GLICO.tar
# Sample_PM1005_TO.tar~done
NAME=$(echo $NAME | cut -f 2,3 -d '_')
SAMP=$(echo $NAME | cut -f 1 -d '.')

NPROC=32
NPROC2=16
NPROC3=10


#loads:
# samtools
# Bowtie2
# pysam
# R
# UMItools

# conda activate ~/miniconda3/envs/UMItools



tar -xvf Sample_${SAMP}.tar


cat ./Sample_${SAMP}/${SAMP}_*_R1_001.fastq.gz > ./Sample_${SAMP}/${SAMP}_R1.fastq.gz
rm ./Sample_${SAMP}/${SAMP}_*_R1_001.fastq.gz
rm ./Sample_${SAMP}/${SAMP}_*_I1_001.fastq.gz
cat ./Sample_${SAMP}/${SAMP}_*_R2_001.fastq.gz > ./Sample_${SAMP}/${SAMP}_R2.fastq.gz
rm ./Sample_${SAMP}/${SAMP}_*_R2_001.fastq.gz

#Gives you the list of significative cells barcodes to analyze.
umi_tools whitelist --method=reads --knee-method=distance \
     -I ./Sample_${SAMP}/${SAMP}_R1.fastq.gz \
     --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
     -S ./Sample_${SAMP}/${SAMP}_R1_WHITELIST.txt

#Exctract all reads with the cell barcodes found in the previous step.
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
     --stdin=./Sample_${SAMP}/${SAMP}_R1.fastq.gz \
     --stdout=./Sample_${SAMP}/${SAMP}_R1_NEW.fastq.gz \
     --read2-in=./Sample_${SAMP}/${SAMP}_R2.fastq.gz \
     --read2-out=./Sample_${SAMP}/${SAMP}_R2_NEW.fastq.gz \
     --filter-cell-barcode \
     --whitelist=./Sample_${SAMP}/${SAMP}_R1_WHITELIST.txt

rm ./Sample_${SAMP}/${SAMP}_R1_NEW.fastq.gz

conda deactivate
conda activate ~/miniconda3/envs/teletest

bowtie2 --no-unal --score-min L,0,1.6 -p ${NPROC} -k 100 \
     --very-sensitive-local \
     -x /athena/nixonlab/scratch/lpr4001/auxfiles/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index \
     -U ./Sample_${SAMP}/${SAMP}_R2_NEW.fastq.gz -S ./Sample_${SAMP}/${SAMP}_bk100.sam \
     --un-gz ./Sample_${SAMP}/${SAMP}_unaligned.gz

rm ./Sample_${SAMP}/${SAMP}_R2_NEW.fastq.gz

samtools view -@ ${NPROC} -b ./Sample_${SAMP}/${SAMP}_bk100.sam > ./Sample_${SAMP}/${SAMP}_bk100.bam
rm ./Sample_${SAMP}/${SAMP}_bk100.sam
samtools sort -@ ${NPROC} ./Sample_${SAMP}/${SAMP}_bk100.bam -o ./Sample_${SAMP}/${SAMP}_bk100_sorted.bam
rm ./Sample_${SAMP}/${SAMP}_bk100.bam
samtools index ./Sample_${SAMP}/${SAMP}_bk100_sorted.bam

conda deactivate

#some python problems, but I modify a script from the scATACutils (https://github.com/crazyhottommy/scATACutils)
#this script takes the bam file and split the mapping reads based on the cell barcode
conda activate ~/miniconda3/envs/scATACutils

python split_SC_bam_by_cell_LPI.py -outdir ./Sample_${SAMP}/temp_bams ./Sample_${SAMP}/${SAMP}_bk100_sorted.bam

conda deactivate
conda activate ~/miniconda3/envs/UMItools

ls ./Sample_${SAMP}/temp_bams > ./Sample_${SAMP}/list_BAMS.txt

gtfHERV=/athena/nixonlab/scratch/lpr4001/auxfiles/transcripts_HERV.gtf
gtfL1=/athena/nixonlab/scratch/lpr4001/auxfiles/transcripts_L1.gtf

mkdir -p ./Telescope/HERV/${SAMP}
mkdir -p ./Telescope/L1/${SAMP}

tot=$(wc -l ./Sample_${SAMP}/list_BAMS.txt| cut -f 1 -d ' ')
#parallizing Telescope, maybe you have better ways or options, this was just to solve the problem and use CPUs

for j in $(seq ${NPROC2} ${NPROC2} ${tot})
do
    for i in $(head -n ${j} ./Sample_${SAMP}/list_BAMS.txt | tail -n ${NPROC2})
    do
        # Reads mapped a the same position with the same UMI are considered PCR duplicates and are removed with UMItools
        samtools index ./Sample_${SAMP}/temp_bams/${i}
        samp=$(echo $i | cut -f 1 -d '.' )
        umi_tools dedup --per-cell \
             --stdin=./Sample_${SAMP}/temp_bams/${i} \
             --extract-umi-method=read_id \
             --umi-separator=_ \
             --method=unique \
             --stdout=./Sample_${SAMP}/temp_bams/${samp}_dedup.bam \
             --no-sort-output --verbose 0

        # Resort the reads for Telescope input
        samtools sort -n ./Sample_${SAMP}/temp_bams/${samp}_dedup.bam | samtools view -F 256 -b -U ./Sample_${SAMP}/temp_bams/${samp}_256.bam -- > ./Sample_${SAMP}/temp_bams/${samp}_016.bam
        samtools merge -n ./Sample_${SAMP}/temp_bams/${samp}_sorted.bam ./Sample_${SAMP}/temp_bams/${samp}_016.bam ./Sample_${SAMP}/temp_bams/${samp}_256.bam

        telescope assign --tempdir ./ --quiet \
            --exp_tag ${samp} --outdir ./Telescope/HERV/${SAMP} \
            --theta_prior 200000 --max_iter 200 \
            --updated_sam ./Sample_${SAMP}/temp_bams/${samp}_sorted.bam ${gtfHERV} &
         telescope assign --tempdir ./ --quiet \
            --exp_tag ${samp} --outdir ./Telescope/L1/${SAMP} \
            --theta_prior 200000 --max_iter 200 \
            --updated_sam ./Sample_${SAMP}/temp_bams/${samp}_sorted.bam ${gtfL1} &
    done
    wait
    #remove temporal files
    for i in $(head -n ${j} ./Sample_${SAMP}/list_BAMS.txt | tail -n ${NPROC2})
    do
        samp=$(echo $i | cut -f 1 -d '.' )
        rm ./Sample_${SAMP}/temp_bams/${samp}*
        rm ./Telescope/HERV/${SAMP}/${samp}-other.bam ./Telescope/HERV/${SAMP}/${samp}-tmp_tele.bam ./Telescope/HERV/${SAMP}/${samp}-checkpoint.npz
        rm ./Telescope/L1/${SAMP}/${samp}-other.bam ./Telescope/L1/${SAMP}/${samp}-tmp_tele.bam ./Telescope/L1/${SAMP}/${samp}-checkpoint.npz
    done
done

# same as above, but with my script some samples were not analyzed. This can be improved for sure!
for i in $(ls ./Sample_${SAMP}/temp_bams)
do
    samtools index ./Sample_${SAMP}/temp_bams/${i}
    samp=$(echo $i | cut -f 1 -d '.' )
    umi_tools dedup --per-cell \
         --stdin=./Sample_${SAMP}/temp_bams/${i} \
         --extract-umi-method=read_id \
         --umi-separator=_ \
         --method=unique \
         --stdout=./Sample_${SAMP}/temp_bams/${samp}_dedup.bam \
         --no-sort-output --verbose 0


    samtools sort -n ./Sample_${SAMP}/temp_bams/${samp}_dedup.bam | samtools view -F 256 -b -U ./Sample_${SAMP}/temp_bams/${samp}_256.bam -- > ./Sample_${SAMP}/temp_bams/${samp}_016.bam
    samtools merge -n ./Sample_${SAMP}/temp_bams/${samp}_sorted.bam ./Sample_${SAMP}/temp_bams/${samp}_016.bam ./Sample_${SAMP}/temp_bams/${samp}_256.bam
    telescope assign --tempdir ./ --quiet \
        --exp_tag ${samp} --outdir ./Telescope/HERV/${SAMP} \
        --theta_prior 200000 --max_iter 200 \
        --updated_sam ./Sample_${SAMP}/temp_bams/${samp}_sorted.bam ${gtfHERV} &
     telescope assign --tempdir ./ --quiet \
        --exp_tag ${samp} --outdir ./Telescope/L1/${SAMP} \
        --theta_prior 200000 --max_iter 200 \
        --updated_sam ./Sample_${SAMP}/temp_bams/${samp}_sorted.bam ${gtfL1} &
done
wait
for i in $(ls ./Sample_${SAMP}/temp_bams)
do
    samp=$(echo $i | cut -f 1 -d '.' )
    rm ./Sample_${SAMP}/temp_bams/${samp}*
    rm ./Telescope/HERV/${SAMP}/${samp}-other.bam ./Telescope/HERV/${SAMP}/${samp}-tmp_tele.bam ./Telescope/HERV/${SAMP}/${samp}-checkpoint.npz
    rm ./Telescope/L1/${SAMP}/${samp}-other.bam ./Telescope/L1/${SAMP}/${samp}-tmp_tele.bam ./Telescope/L1/${SAMP}/${samp}-checkpoint.npz
done




conda deactivate

#Alevin is a subprogram of SALMON desinged for single cell RNA-seq (https://salmon.readthedocs.io/en/latest/alevin.html)

conda activate ~/miniconda3/envs/SALMON

#Too many cells causes a lot of memory consuption, therfore I splited the file.

cut -f 1 ./Sample_${SAMP}/${SAMP}_R1_WHITELIST.txt | split -l 1000 - ./Sample_${SAMP}/${SAMP}_R1_WHITELIS_splited_

mkdir -p ./Alevin/${SAMP}/

for i in $(ls ./Sample_${SAMP}/${SAMP}_R1_WHITELIS_splited_*)
do
    salmon alevin -l ISR -1 ./Sample_${SAMP}/${SAMP}_R1.fastq.gz -2 ./Sample_${SAMP}/${SAMP}_R2.fastq.gz \
        --chromium -i /athena/nixonlab/scratch/lpr4001/auxfiles/gencode_V31_Salmon_index \
        --tgMap /athena/nixonlab/scratch/lpr4001/auxfiles/gencode_V31_Salmon_index/txp2gen.tsv \
        –-dumpMtx \
        -p ${NPROC3} -o ./Alevin/${SAMP}/${i} --whitelist ${i}
done

#this is for joining tables of Alevin and Telescope outputs (HERVs and L1)

Rscript ./joinTelescope.Rscript ${SAMP}

rm ./Sample_${SAMP}/${SAMP}_R1.fastq.gz ./Sample_${SAMP}/${SAMP}_R2.fastq.gz
rm ./Sample_${SAMP}/${SAMP}_R1_WHITELIS_splited_*
#rm -rf ./Alevin/${SAMP}/Sample_${SAMP}

conda deactivate
