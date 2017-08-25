#!/bin/zsh
# First argument is the name of the reference sequence
# Second argument is the name of the fw query files
# Third argument is the name of the rc query files
# Fourth argument is the seed length

# First, we build the index
bowtie2-build $1 index
# run bowtie on the fw reads
bowtie2 -x index -U $2 \
	-N 0 \
	-L $4 \
	-a \
	--no-1mm-upfront \
	-i L,0,1 \
	--local \
	--score-min L,$4,0 \
	--mp 50,50 \
	--ma 1 \
	--norc \
	--ignore-quals \
	-S output_bt_fw.sam
# run bowtie on the rc reads
bowtie2 -x index -U $3 \
	-N 0 \
	-L $4 \
	-a \
	--no-1mm-upfront \
	-i L,0,1 \
	--local \
	--score-min L,$4,0 \
	--mp 50,50 \
	--ma 1 \
	--nofw \
	--ignore-quals \
	-S output_bt_rc.sam


# use samtools to convert the sam to bam
samtools view -bS output_bt_fw.sam > output_bt_fw.bam
samtools sort output_bt_fw.bam -o output_bt_fw.bam -O BAM
samtools index output_bt_fw.bam output_bt_fw.bam.bai
samtools view -bS output_bt_rc.sam > output_bt_rc.bam
samtools sort output_bt_rc.bam -o output_bt_rc.bam -O BAM
samtools index output_bt_rc.bam output_bt_rc.bam.bai
