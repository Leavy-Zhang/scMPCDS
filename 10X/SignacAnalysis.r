library(Signac)
library(Seurat)

mito.data <- ReadMGATK(dir = "major_result/mgatk/final/")
valid_cells <- rownames(mito.data$depth)[mito.data$depth$mito.depth >= 20]
mito_assay <- CreateAssayObject(counts = mito.data$counts)



counts <- Read10X_h5(filename = "major_result/10X-atac/outs/filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "major_result/10X-atac/outs/fragments2.tsv.gz",
  min.cells = 10,
  min.features = 200
)


atac_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

atac_obj <- RunTFIDF(atac_obj)
atac_obj <- FindTopFeatures(atac_obj, min.cutoff = 10)
atac_obj <- RunSVD(atac_obj)
saveRDS(atac_obj, 'major_result/10X-atac/outs/object.rds')



atac_obj[['mito']] <- mito_assay
atac_obj2 <-subset(atac_obj,cells = valid_cells)
variable.sites <- IdentifyVariants(
  object = atac_obj,
  assay = "mito",
  refallele = mito.data$refallele
)

atac_obj <- AlleleFreq(
object = atac_obj,
variants = subset(variable.sites, subset = (n_cells_conf_detected >=1))$variant,
assay = "mito"
)


pdf('figures/10X.example.pdf')
DefaultAssay(atac_obj) <- "alleles"
alleles.view <- c("2149G>A", "15259C>T", "16399A>G")
FeaturePlot(
object = atac_obj,
features = alleles.view,
order = TRUE,
cols = c("grey", "darkred"),
ncol = 4
) & NoLegend()

dev.off()
