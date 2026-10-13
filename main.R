suppressPackageStartupMessages({
  library(ORFik)
  library(GenomicFeatures)
  library(txdbmaker)
  library(AnnotationDbi)
  library(BiocParallel)
  library(Rsamtools)
  library(GenomeInfoDb)
})

# ============================================================
# FICHIERS
# ============================================================

bam_dir <- "/project/def-marie87/glaudea/riboseq_analysis/bams/50"
genome_file <- "/project/def-marie87/glaudea/riboseq_analysis/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# IMPORTANT : utiliser le GTF Ensembl SANS chr
gtf_file <- "/project/def-marie87/glaudea/riboseq_analysis/data/Homo_sapiens.GRCh38.114.gtf"

txdb_file <- "/project/def-marie87/glaudea/riboseq_analysis/data/Homo_sapiens.GRCh38.114.db"

output_dir <- "/project/def-marie87/glaudea/riboseq_analysis/results/50"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. INDEX FASTA
# ============================================================

if (!file.exists(paste0(genome_file, ".fai"))) {
  Rsamtools::indexFa(genome_file)
}


# ============================================================
# 2. SEQINFO AVEC LONGUEURS DU GÉNOME
# ============================================================

fa_index <- Rsamtools::scanFaIndex(genome_file)

chrominfo <- GenomeInfoDb::Seqinfo(
  seqnames = as.character(seqnames(fa_index)),
  seqlengths = width(fa_index)
)

cat("\nSeqinfo FASTA :\n")
print(chrominfo[seqlevels(chrominfo)[1:min(30, length(chrominfo))]])
cat("\nNombre de séquences FASTA :", length(chrominfo), "\n")

cat("\nNombre de seqlengths NA dans FASTA :",
    sum(is.na(seqlengths(chrominfo))), "\n")


# ============================================================
# 3. TxDb AVEC LONGUEURS CHROMOSOMIQUES
# ============================================================

txdb <- txdbmaker::makeTxDbFromGFF(
  gtf_file,
  format = "gtf",
  organism = "Homo sapiens",
  chrominfo = chrominfo
)


# ============================================================
# 4. VÉRIFICATION DU TxDb
# ============================================================

cat("\nSeqlevels TxDb :\n")
print(head(seqlevels(txdb), 30))

sl_txdb <- seqlengths(txdb)

cat("\nLongueurs TxDb :\n")
print(sl_txdb)

cat("\nNombre de seqlengths NA dans TxDb :",
    sum(is.na(sl_txdb)), "\n")

cat("\nNombre de seqlengths <= 0 dans TxDb :",
    sum(!is.na(sl_txdb) & sl_txdb <= 0), "\n")


# ============================================================
# 5. SAUVEGARDE TxDb
# ============================================================

if (file.exists(txdb_file)) {
  file.remove(txdb_file)
}

AnnotationDbi::saveDb(
  txdb,
  file = txdb_file
)

cat("\nTxDb sauvegardé :\n", txdb_file, "\n")

# ============================================================
# 2. EXPÉRIENCE ORFik
# ============================================================

create.experiment(
  dir = bam_dir,
  exper = "50",
  fa = genome_file,
  txdb = txdb_file,
  organism = "Homo sapiens",
  pairedEndBam = FALSE,
  rep = 1,
  condition = "wt",
  saveDir = output_dir,
  viewTemplate = FALSE
)

df.rfp <- read.experiment("50", output_dir)

print(df.rfp)


trace(
  data.table:::as.data.table.list,
  quote({
    if (
      "levels" %in% names(x) &&
      all(c("seqnames", "start", "cigar", "strand") %in% names(x)) &&
      length(x$levels) == 3
    ) {
      x$strand <- factor(
        as.character(x$strand),
        levels = c("+", "-", "*")
      )
      x$levels <- NULL
    }
  }),
  print = FALSE
)
BPPARAM <- BiocParallel::SerialParam()
ORFikQC(df.rfp)
# ============================================================
# 3. P-SHIFT
# ============================================================

shiftFootprintsByExperiment(df.rfp,accepted.lengths = 25:40)
