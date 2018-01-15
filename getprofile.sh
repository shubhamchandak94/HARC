#/bin/bash

### Base-Calling Profile ###

dirname=assembly_data/datasets
for dataset in ERR204808 SRR065390 SRR823377 SRR543736
do
	prefix=$dirname"/"$dataset
	ref_fasta=$prefix".genome.fasta"
	# Or Using BWA
	#./../../variantcalling/bwa-0.7.16a/bwa mem -t 16 -M $ref_fasta $prefix"_1.fastq" $prefix"_2.fastq" > $prefix".sam" 

	./../pirs/baseCalling_Matrix_calculator -r $ref_fasta -o $prefix".baseCallProfile" -i $prefix".sam"
	#./baseCalling_Matrix_calculator -b -r ./ref/human.fa -i bam/GA0230.sort.bam -s vcf/GA0230.SNPs.filter.vcf -o matrix/GA0230 2>matrix/GA0230.matrix.err
done
