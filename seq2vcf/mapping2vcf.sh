
conda activate Bioinfo

threads=6
ref_base="/pub/ref_base"
fasta_file="Homo_Ref_Numts.fa"
bwa index -a bwtsw "$ref_base/$fasta_file"

mkdir vcf_file 
cd sample_dir
mkdir clean_data
for i in `ls  `; do
    fastp -i 01.RawData/"$i"/"$i"_1.fq.gz \
        -I 01.RawData/"$i"/"$i"_2.fq.gz \
        -o clean_data/"$i"/"$i".R1.clean.fastq.gz \
        -O clean_data/"$i"/"$i".R2.clean.fastq.gz \
        -w 16 -h clean_data/"$i"/"$i".html \
        -j clean_data/"$i"/"$i".json \
        2> clean_data/"$i"/"$i".fastp.log

    pTrimmer-1.3.1 -s single -a hmp_mtDNA.txt \
        -f clean_data/"$i"/"$i".R1.clean.fastq.gz -r clean_data/"$i"/"$i".R2.clean.fastq.gz \
        -d clean_data/"$i"/"$i".R1.no_primer.fastq.gz -e clean_data/"$i"/"$i".R2.no_primer.fastq.gz -z -m 2


    bwa mem "$ref_base"/"$fasta_file" \
        clean_data/"$i"/"$i".R1.no_primer.fastq.gz \
        clean_data/"$i"/"$i".R2.no_primer.fastq.gz -t "$$threads" 2> clean_data/"$i"/"$i".mapping.log  | samtools view -@ "$$threads" -bS > clean_data/"$i"/"$i".bam

    samtools view -b clean_data/"$i"/"$i".bam "chrM" > clean_data/"$i"/"$i".mito.bam
    samtools sort clean_data/"$i"/"$i".mito.bam -o clean_data/"$i"/"$i".mito.sort.bam
    samtools index clean_data/"$i"/"$i".mito.sort.bam
    samtools mpileup -B -q 20 -Q 30 -d 1000000 -f "$ref_base"/Homo_chrM_Ref_Numts.fa clean_data/"$i"/"$i".mito.sort.bam |\
        varscan mpileup2snp --output-vcf 1 --min-var-freq 0.01 --strand-filter 0 --min-coverage 5 > ../"$vcf_file"/"$i".vcf"
done

cd ../vcf_file
for i in `ls | grep vcf$`;do
    java -Xmx100G -jar ~/tools/snpEff/snpEff.jar -c ~/tools/snpEff/snpEff.config -v -o gatk GRCh38.99 "$" > "${i%vcf}anno.vcf.gz";tabix -p vcf "$i"/"${j%vcf}anno.vcf"
done




