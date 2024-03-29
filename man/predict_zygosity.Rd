% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fncts.R, R/fncts_local.R
\name{predict_zygosity}
\alias{predict_zygosity}
\title{predicts zygosity of a set of genes of a sample}
\usage{
predict_zygosity(
  purity,
  sex,
  somCna,
  geneModel,
  bamDna,
  somSmallVars = NULL,
  germSmallVars = NULL,
  bamRna = NULL,
  ploidy = NULL,
  colnameTcn = NULL,
  colnameCnaType = NULL,
  includeHomoDel = TRUE,
  includeIncompleteDel = TRUE,
  showReadDetail = FALSE,
  printLog = FALSE,
  assumeSomCnaGaps = FALSE,
  byTcn = TRUE,
  vcf = NULL,
  distCutOff = 5000
)

predict_zygosity(
  purity,
  sex,
  somCna,
  geneModel,
  bamDna,
  somSmallVars = NULL,
  germSmallVars = NULL,
  bamRna = NULL,
  ploidy = NULL,
  colnameTcn = NULL,
  colnameCnaType = NULL,
  includeHomoDel = TRUE,
  includeIncompleteDel = TRUE,
  showReadDetail = FALSE,
  printLog = FALSE,
  assumeSomCnaGaps = FALSE,
  byTcn = TRUE,
  vcf = NULL,
  distCutOff = 5000
)
}
\arguments{
\item{purity}{purity of the sample (numeric value between 0 and 1 indicating 
the fraction of relevant sample with control/unrelevant tissue)}

\item{sex}{sex of the sample (character: "male", "female", "m", "f")}

\item{somCna}{GRanges object containing all genomic regions with annotated 
total copynumber and cna_type as metadata columns. The total-copynumber 
column should be named "tcn" but also some other commonly used names. 
It should contain numeric values or characters that can be converted to 
numeric values. The cna_type column must contain the information about 
loss of heterozygosity (LOH). Therefore the term "LOH" must be explicitely 
mentioned in the column. If a genomic region is not present in the object, 
it will be taken as heterozygous with neutral TCN of 2.}

\item{geneModel}{GRanges object containing the gene-annoattion of 
the used reference genome with metadata column of the gene name (gene)}

\item{bamDna}{path to bam-file}

\item{somSmallVars}{GRanges object containing all somatic small 
variants (SNV and INDEL).
Required metadata columns are reference base (ref/REF), 
alternative base (alt/ALT),
annotation of the gene name (gene/GENE) and the allele-frequency (af/AF). 
If the object is not provided the tool assumes there are no somatic small 
variants.}

\item{germSmallVars}{GRanges object containing all germline small 
variants (SNV and INDEL).
Required metadata columns are reference base (ref/REF), alternative 
base (alt/ALT),
annotation of the gene name (gene/GENE) and the allele-frequency (af/AF)
If the object is not provided the tool assumes there are no germline small 
variants.}

\item{bamRna}{optional; path to rna file (bam format)}

\item{ploidy}{ploidy of the sample (numeric value)}

\item{colnameTcn}{character indicating the name of the metadata containing 
the tcn information in the somCna object. If not provided the tool tries to 
detect the column according to default names}

\item{colnameCnaType}{character indicating the name of the metadata 
containing cna type information in the somCna object. 
If not provided the tool tries to detect the column according to default 
names}

\item{includeHomoDel}{default = TRUE; if FALSE homozygous deleteions are 
excluded}

\item{includeIncompleteDel}{default = TRUE; if FALSE heterzygous deleteions 
are excluded}

\item{showReadDetail}{default = FALSE; if TRUE a table is added to the 
output,}

\item{printLog}{default = FALSE; if TRUE the gene which is evaluated is 
printed in console, 
containing the query-name of each read which was used to perform 
haplotype-phasing and the info into which class it was assigned.}

\item{assumeSomCnaGaps}{(logical, default=FALSE) Only required if the somCna
object lacks copy number information for genomic segments on which small 
variants are detected. By default, variants in such regions will be excluded 
from the analysis as required information about the copy number is missing. 
These variants will be attached to the final output list in a separate 
tibble. To include them, this flag must be set TRUE and the ground ploidy 
must be given as an input. This ground ploidy will then be taken as tcn in 
the missing regions. If no ploidy is given the tool will assume the ground 
ploidy of 2 when this flag is TRUE.}

\item{byTcn}{logical, default=TRUE; optional if includeHomoDel or 
includeIncompleteDelS is TRUE. If FALSE the tool will not use tcn as a 
criterion to assign large deletions. It will use the cna_type column and 
check for indicating strings like HOMDEL/HomoDel/DEL. Some commonly used 
strings are covered. It is recommended to leave this flag TRUE}

\item{vcf}{character; path to variant call file (.vcf.gz format). 
Will be used (if provided)
for extended SNP phasing if variants on the same gene are too far away from
each other for direct haplotype phasing}

\item{distCutOff}{numeric, default=5000; if input vcf is provided and SNP
phasing is performed, this will limt the distance at which the SNP phasing
should not be tried anymore. As the probability of finding overlapping reads
at such a long distance is very low and the runtime will increase
exponentially.}
}
\value{
A list of dataframes. Those are the evaluation per variant, 
the evaluation per gene and, if performed, the info about the 
haplotype-phasing.

A list of dataframes. Those are the evaluation per variant, 
the evaluation per gene and, if performed, the info about the 
haplotype-phasing.
}
\description{
predicts zygosity of a set of genes of a sample

predicts zygosity of a set of genes of a sample
}
\examples{
cnvs  = GenomicRanges::GRanges(
  dplyr::tibble(
    chr = "chr17",
    start = c(170060, 34520990),
    end = c(34520990, 83198614),
    tcn = c(2, 1),
    cna_type = c("neutral", "LOH")
  )
)
somatic_vars = GenomicRanges::GRanges(
  dplyr::tibble(
    chr="chr17",
    start = 7675088,
    end = 7675088,
    ref = "C",
    alt = "T",
    af = 0.65,
    gene = "TP53" 
  )
)
germline_vars = GenomicRanges::GRanges(
  dplyr::tibble(
    chr="chr17",
    start = 41771694,
    end = 41771694,
    ref = "GTGT",
    alt = "G",
    af = 0.95,
    gene = "JUP" 
  )
)
reference = GenomicRanges::GRanges(
  dplyr::tibble(
    chr = "chr17",
    start = c(7661778, 41754603),
    end = c(7687538, 41786931),
    gene = c("TP53", "JUP")
  )
)
sex = "female"
purity = 0.9
bamfile <- system.file("extdata", "ZP_example.bam", 
  package = "ZygosityPredictor")
predict_zygosity(purity = purity, sex = sex, 
  somCna = cnvs,
  somSmallVars = somatic_vars,
  germSmallVars = germline_vars,
  geneModel = reference,
  bamDna = bamfile
)
cnvs  = GenomicRanges::GRanges(
  dplyr::tibble(
    chr = "chr17",
    start = c(170060, 34520990),
    end = c(34520990, 83198614),
    tcn = c(2, 1),
    cna_type = c("neutral", "LOH")
  )
)
somatic_vars = GenomicRanges::GRanges(
  dplyr::tibble(
    chr="chr17",
    start = 7675088,
    end = 7675088,
    ref = "C",
    alt = "T",
    af = 0.65,
    gene = "TP53" 
  )
)
germline_vars = GenomicRanges::GRanges(
  dplyr::tibble(
    chr="chr17",
    start = 41771694,
    end = 41771694,
    ref = "GTGT",
    alt = "G",
    af = 0.95,
    gene = "JUP" 
  )
)
reference = GenomicRanges::GRanges(
  dplyr::tibble(
    chr = "chr17",
    start = c(7661778, 41754603),
    end = c(7687538, 41786931),
    gene = c("TP53", "JUP")
  )
)
sex = "female"
purity = 0.9
bamfile <- system.file("extdata", "ZP_example.bam", 
  package = "ZygosityPredictor")
predict_zygosity(purity = purity, sex = sex, 
  somCna = cnvs,
  somSmallVars = somatic_vars,
  germSmallVars = germline_vars,
  geneModel = reference,
  bamDna = bamfile
)
}
