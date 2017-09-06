#!/bin/zsh
# First argument is the name of the reference sequence
# Second argument is the prefix for the query files
# Third argument is the seed length
# Fourth argument is the prefix for the output files

# First, we build the index
bowtie2-build $1 index
# run bowtie on the forward strand to the right of the seed
bowtie2 -x index -U $2_fw.fq \
	-N 0 \
	-L $3 \
	-a \
	--no-1mm-upfront \
	-i L,0,2 \
	--local \
	--score-min L,$3,0 \
	--mp 50,50 \
	--ma 1 \
	--norc \
	--ignore-quals \
	-S $4_fw_right.sam
# run bowtie on the fw strand to the left of the seed
bowtie2 -x index -U $2_rc.fq \
	-N 0 \
	-L $3 \
	-a \
	--no-1mm-upfront \
	-i L,0,2 \
	--local \
	--score-min L,$3,0 \
	--mp 50,50 \
	--ma 1 \
	--nofw \
	--ignore-quals \
	-S $4_fw_left.sam
# run bowtie on the rc strand to the right of the seed
bowtie2 -x index -U $2_rc.fq \
	-N 0 \
	-L $3 \
	-a \
	--no-1mm-upfront \
	-i L,0,2 \
	--local \
	--score-min L,$3,0 \
	--mp 50,50 \
	--ma 1 \
	--norc \
	--ignore-quals \
	-S $4_rc_right.sam
# run bowtie on the rc strand to the left of the seed
bowtie2 -x index -U $2_fw.fq \
	-N 0 \
	-L $3 \
	-a \
	--no-1mm-upfront \
	-i L,0,2 \
	--local \
	--score-min L,$3,0 \
	--mp 50,50 \
	--ma 1 \
	--nofw \
	--ignore-quals \
	-S $4_rc_left.sam

for suffix in fw_right fw_left rc_right rc_left; do
    samtools view -bS $4_$suffix.sam > $4_$suffix.bam
    samtools sort $4_$suffix.bam -o $4_$suffix.bam -O BAM
    samtools index $4_$suffix.bam $4_$suffix.bam.bai
done
