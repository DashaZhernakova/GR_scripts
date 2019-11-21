bam_list=$1
cat $bam_list | \
	while read line
	do
		tmp=${line##*/}
		sample_id=${tmp%.bam}
		out_dir="/mnt/gr2/HLA_genotyping/results/${sample_id}/"
		mkdir $out_dir
		echo "running HLA-LA for $sample_id, bam file: $line"
		echo "output folder: $out_dir" 
		/home/genomerussia/tools/HLA-LA/src/HLA-LA.pl --BAM $line --graph PRG_MHC_GRCh38_withIMGT --sampleID $sample_id --maxThreads 50 > ${out_dir}/run.log 2>&1
	
	done

echo "Finished all"