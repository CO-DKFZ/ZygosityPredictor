


#library(getopt)
library(dplyr)
library(readr)
library(stringr)
library(GenomicRanges)
library(ZygosityPredictor)


# spec = matrix(c(
#   "germline",          "g", 0, "character",   "Print this help message and quit.",
#   "somatic",           "s", 0, "character",   "Print this help message and quit.",
#   "cnv",               "c", 0, "character",   "Print this help message and quit.",
#   "meta",              "m", 0, "character",   "Print this help message and quit.",
#   "file_gene_model",   "b", 2, "character",   "Path to exon annotation BED file for gene model. Default is.",
#   "pre_proc_fncts",    "p", 0, "character",   "Print this help message and quit.",
#   "file_sel_fncts",    "f", 0, "character",   "Print this help message and quit.",
#   #"wf_fncts",    "w", 0, "character",   "Print this help message and quit.",
#   "run_mode",    "z", 0, "character",   "Print this help message and quit.",
#   "otp_paths",    "o", 0, "character",   "Print this help message and quit."
# ),ncol=5,byrow=TRUE)
# opt <- getopt(spec)

PID <- "WGS.1174YV.metastasis02"

input_file <- tibble(
  germ=file.path("/omics/odcf/analysis/OE0246_projects/decisions/shared/results_per_pid/", PID, "/dMextraction/extracted_smallvar_germline.tsv"),
  som=file.path("/omics/odcf/analysis/OE0246_projects/decisions/shared/results_per_pid/", PID, "/dMextraction/extracted_smallvar_somatic.tsv"),
  cnv=file.path("/omics/odcf/analysis/OE0246_projects/decisions/shared/results_per_pid/", PID, "/dMextraction/extracted_cnvs.tsv"),
  meta=file.path("/omics/odcf/analysis/OE0246_projects/decisions/shared/results_per_pid/", PID, "/dMextraction/extracted_metainfo.tsv"),
  otp=file.path("/omics/odcf/analysis/OE0246_projects/decisions/shared/results_per_pid/", PID, "/otp_paths/otp_paths.csv")
)

opt <- tibble(

  file_gene_model="/omics/odcf/project/ODCF/reference_genomes/bwa06_1KGRef_PhiX/gencode/gencode19/GencodeV19_Exons_plain.bed.gz",
  pre_proc_fncts="/omics/groups/OE0660/internal/m168r/projects/TopArtCalc/TAC_preprocessing/pre_proc_fncts.R",
  file_sel_fncts="/omics/groups/OE0660/internal/m168r/projects/TopArtCalc/TAC_file_selection/file_sel_fncts.R",
  germline=input_file$germ[1],
  somatic=input_file$som[1],
  cnv=input_file$cnv[1],
  meta=input_file$meta[1],
  otp_paths=input_file$otp,
  run_mode=Sys.getenv("ZP_RUN_MODE", "full_pid_test")

)
#print(opt)



geneModel <- read_tsv(
  opt$file_gene_model,
  col_types=cols(.default="c"),
  col_names=c("chr", "start", "end","gene", "n_exon", "strand", "cdsstart", "cdsend"),
  skip=1,
  show_col_types = FALSE ) %>%
  #filter(gene %in% c("ZNF99", "TP53", "BRCA1", "BRCA2", "PTEN", "FBXW4", "CPT1A")) %>%
  GRanges()

cat("gene model loaded")

germSmallVars <- read_tsv(opt$germline,
                          col_types=cols(.default="c")) %>%
  GenomicRanges::GRanges()
somSmallVars <- read_tsv(opt$somatic,
                         col_types=cols(.default="c")) %>%
  GenomicRanges::GRanges()
somCna <- read_tsv(opt$cnv,
                   col_types=cols(.default="c")) %>%
  mutate_at(.vars = c("TCN"),
            .funs = as.numeric) %>% select(-id_dataMASTER) %>%
  GenomicRanges::GRanges()
meta <- read_tsv(opt$meta, col_types=cols(.default="c")) %>%
  mutate_at(.vars = c("TumorCellContent", "Ploidy"),
            .funs = as.numeric)

source(opt$pre_proc_fncts)
source(opt$file_sel_fncts)


otp_paths <- read_csv(opt$otp_paths)

print(otp_paths)

bamDna <- otp_paths$bamDna

id_dataMASTER <- paste(meta$DNASeq, meta$PatientID, meta$TumorID,sep=".")
too_much_rna_samples <- c("WGS.89D53Z.tumor", "WES.1KJALE.metastasis", 
                          "WGS.5MRTGM.tumor", "WGS.BG2P4P.metastasis03", "WGS.3Q1VYU.metastasis07")



if(id_dataMASTER %in% too_much_rna_samples){
  bamRna <- NULL
  raw_vcf_file <- NULL
  phasedVcfs <- NULL
  haploBlocks <- NULL
} else {
  bamRna <- otp_paths$bamRna
  raw_vcf_file <- otp_paths$raw_snv_file
  phasing_dir <- otp_paths$phasing_dir  
  if(is.na(bamRna)){
    bamRna <- NULL
  }
  
  if(meta$DNASeq=="WES"){
    phasedVcfs <- raw_vcf_file
    haploBlocks <- NULL
  } else {
    haploblock_files <- list.files(phasing_dir, pattern="haploblock", full.names = T) 
    cat("hploblocks files detected")
    
    haploBlocks <- prepare_haploblock_object(haploblock_files)
    cat("hploblocks object created")
    phasedVcfs <- get_phased_vcfs(phasing_dir)
  }
}



# if(is.na(bamRna)){
#   bamRna <- NULL
# }
# 
# if(meta$DNASeq=="WES"){
#   phasedVcfs <- raw_vcf_file
#   haploBlocks <- NULL
# } else {
#   haploblock_files <- list.files(phasing_dir, pattern="haploblock", full.names = T) 
#   cat("hploblocks files detected")
#   
#   haploBlocks <- prepare_haploblock_object(haploblock_files)
#   cat("hploblocks object created")
#   phasedVcfs <- get_phased_vcfs(phasing_dir)
# }

start_time <- Sys.time()
zp_workers <- as.integer(Sys.getenv("ZP_WORKERS", "1"))
if(is.na(zp_workers)||zp_workers < 1){
  zp_workers <- 1
}
bp_param <- if(zp_workers > 1){
  BiocParallel::MulticoreParam(workers=zp_workers)
} else {
  BiocParallel::SerialParam()
}

full_prediction <- predict_zygosity(
  purity=meta$TumorCellContent,
  sex=meta$Sex,
  somCna = somCna,
  geneModel = geneModel,
  bamDna = bamDna,
  somSmallVars = somSmallVars,
  germSmallVars = germSmallVars,
  #germSmallVars=NULL,
  bamRna = bamRna,
  ploidy = meta$Ploidy,
  haploBlocks=haploBlocks,
  vcf=phasedVcfs,
  ## options for run
  printLog = T,
  includeIncompleteDel = T,
  includeHomoDel = T,
  verbose=F,
  debug=F,
  ## highest distance of two variants that should be phased (if reduced runtime reduces as well)
  distCutOff = 2500,
  ## if provided, detailed information will be stored there
  logDir=NULL,
  showReadDetail = F,
  ## Secondary phasing approach if read-level phasing fails... can lead to wroing results
  AllelicImbalancePhasing=F,
  byTcn=TRUE,
  BPPARAM=bp_param
)

end_time <- Sys.time()

ZP_version <- as.character(packageVersion("ZygosityPredictor"))

outpath <- file.path("/omics/groups/OE0660/internal/m168r/ZP_test", paste0(PID, "_ZP_", ZP_version))
dir.create(outpath)

store_everything(full_prediction, outpath)



write_tsv(tibble(success="TRUE", time=as.character(end_time-start_time), version=ZP_version), 
          file = 
            file.path(outpath, 
                      paste("ZygosityPredictor", ZP_version, paste0(opt$run_mode,".suc"), sep="_")))


