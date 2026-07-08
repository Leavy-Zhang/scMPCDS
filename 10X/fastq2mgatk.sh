
./cellranger-atac-2.2.0/bin/cellranger-atac count --id 10X-mito-atac \
        --reference ref_base/refdata-cellranger-arc-GRCh38-2024-A \
        --project 10X-mito-atac --fastqs 10X_rawdata \
        --sample 201647G_293T --localcores 50 --localmem 500

samtools view -@ 30 -b -h 10X-mito-atac/outs/possorted_bam.bam chrM > 10X-mito-atac/outs/chrM.possorted_bam.bam
samtools index 10X-mito-atac/outs/chrM.possorted_bam.bam
mgatk tenx -i 10X-mito-atac/outs/chrM.possorted_bam.bam -g hg38 -o major_results/mgatk -c 40 -bt CB -b 10X-mito-atac/outs/filtered_peak_bc_matrix/barcodes.tsv
