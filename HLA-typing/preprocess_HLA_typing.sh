bam_path=/home/genomerussia/bakeoff/macrogen/bam/
out_path=/home/dzhernakova/bakeoff/macrogen/HLA_typing/results_pairend/

for f in ${bam_path}*
do
id=${f##*/}
echo $id
mkdir ${out_path}/${id}

/home/genomerussia/tools/samtools-1.3/samtools view -hb ${f}/${id}.sorted.bam 6:28000000-34000000 > ${out_path}/${id}/${id}.subset.bam
/home/genomerussia/tools/samtools-1.3/samtools sort -n ${out_path}/${id}/${id}.subset.bam | \
/home/genomerussia/tools/samtools-1.3/samtools fixmate - ${out_path}/${id}/${id}.subset.qsort.fixed.bam

/home/genomerussia/tools/samtools-1.3/samtools fastq -f 3 -1 ${out_path}/${id}/${id}.HLA_1.fq -2 ${out_path}/${id}/${id}.HLA_2.fq -0 ${out_path}/${id}/${id}.HLA_0.fq  ${out_path}/${id}/${id}.subset.qsort.fixed.bam 

rm ${out_path}/${id}/${id}.subset.bam
rm ${out_path}/${id}/${id}.subset.qsort.fixed.bam

/home/genomerussia/tools/samtools-1.3/samtools view -f4 -h ${f}/${id}.sorted.bam | /home/genomerussia/tools/samtools-1.3/samtools fastq - | gzip - > ${out_path}/${id}/${id}.unnmapped.fq.gz

echo "#"
echo "#"
echo "#"
echo "#"

wc ${out_path}/${id}/${id}.HLA_1.fq
wc ${out_path}/${id}/${id}.HLA_2.fq

echo "#"
echo "#"
echo "#"
echo "#"

/home/genomerussia/tools/bwa.kit-0.7.15/bwa mem \
-t 24 \
~/tools/Athlates_2014_04_26/db/ref/hla.clean.fasta \
${out_path}/${id}/${id}.HLA_1.fq \
${out_path}/${id}/${id}.HLA_2.fq \
| /home/genomerussia/tools/samtools-1.3/samtools view -bSh | \
/home/genomerussia/tools/samtools-1.3/samtools sort -@10 - > ${out_path}/${id}/${id}.HLA_bwa.bam

/home/genomerussia/tools/bowtie2-2.2.8/bowtie2 \
-p 24 \
-X 800 \
-x ~/resources/MHC_typing/hla.clean \
-1 ${out_path}/${id}/${id}.HLA_1.fq \
-2 ${out_path}/${id}/${id}.HLA_2.fq \
-U ${out_path}/${id}/${id}.HLA_0.fq \
-U ${out_path}/${id}/${id}.unnmapped.fq.gz \
| /home/genomerussia/tools/samtools-1.3/samtools view -bhS \
| /home/genomerussia/tools/samtools-1.3/samtools sort -@10 - > ${out_path}/${id}/${id}.HLA_bowtie.bam
GR0001.HLA+unmapped_1.fq

echo "Extract genes"
genes=("A" "B" "C" "DQA1" "DQB1" "DRB1" "DRB2" "DRB3" "DRB4" "DRB5" "DRB6" "DRB7" "DRB8" "DRB9")
for gene in ${genes[@]}
do
echo "$gene"
/home/genomerussia/tools/samtools-1.3/samtools view -bh -F4 -L ~/tools/Athlates_2014_04_26/db/bed/hla.${gene}.bed \
${out_path}/${id}/${id}.HLA_bowtie.bam \
| /home/genomerussia/tools/samtools-1.3/samtools sort -n - \
> ${out_path}/${id}/${id}.HLA_bowtie.bam.hla.${gene}.qsorted.bam

/home/genomerussia/tools/samtools-1.3/samtools view -bh -F4 -L ~/tools/Athlates_2014_04_26/db/bed/hla.non-${gene}.bed \
${out_path}/${id}/${id}.HLA_bowtie.bam \
| /home/genomerussia/tools/samtools-1.3/samtools sort -n - \
> ${out_path}/${id}/${id}.HLA_bowtie.bam.hla.non${gene}.qsorted.bam

/home/genomerussia/tools/samtools-1.3/samtools view -bh -F4 -L ~/tools/Athlates_2014_04_26/db/bed/hla.${gene}.bed \
${out_path}/${id}/${id}.HLA_bwa.bam \
| /home/genomerussia/tools/samtools-1.3/samtools sort -n - \
> ${out_path}/${id}/${id}.HLA_bwa.bam.hla.${gene}.qsorted.bam

/home/genomerussia/tools/samtools-1.3/samtools view -bh -F4 -L ~/tools/Athlates_2014_04_26/db/bed/hla.non-${gene}.bed \
${out_path}/${id}/${id}.HLA_bwa.bam \
| /home/genomerussia/tools/samtools-1.3/samtools sort -n - \
> ${out_path}/${id}/${id}.HLA_bwa.bam.hla.non${gene}.qsorted.bam

done
done