#' Reads data from a MAF or VCF type file.
#' This is a wrapper for read.delim that reads mutation data from a tab delimited text file such as a MAF or VCF file.
#' The input file must have at least 5 columns specifying in order 1:5 chromosome, start position, end position, reference sequence, and the mutant sequence. 
#' For these any header names will be ignored and the imported data.frame will have the columns chr, start, end, ref, and mut.
#' If there are additonal columns specifying strand or sampleID these can be imported too either by name (if there is  a header) or number.
#' other columns of data will be ignored.
#' @param file a delimited file with columns containing at least "chr", "start.position", and "end.position" (a data.frame)
#' @param the name of an optional strand column
#' @param WT name of WT sequence column
#' @param MUT name of Mutant sequence column
#' @param the name of an optional sample ID column
#' @param skip VCF or MAF files often have a few rows of comments that you will want to skip over (numeric, default=0).
#' @keywords rainfall, mutants, kataegesis
#' @examples
#' # A single example data file comes with the package
#' (exF <- system.file("extdata", "1317049.tsv", package="kat"))  
#' # a minimal example
#' PD4107a_min <- read.mutations(file=exF, sampleID = Sample.Name)

read.mutations <- function(file,  strand=strand, sampleID = sampleID, skip=0, header=T, ...){
  
  # Non Standard Evaluation: see http://adv-r.had.co.nz/Computing-on-the-language.html#capturing-expressions
  strand=deparse(substitute(strand))
  sampleID = deparse(substitute(sampleID))
  temp = read.delim(file, as.is=T, ...)
  
  # check that at least the 5 main columns exist - fail if not.
  if(ncol(temp) > 5 && nrow(temp) > 1)
  {
    TSV
    VCF <- data.frame(chr= temp[,chr], start.position= temp[,start.position], end.position= temp[,end.position])
  }
  else
  {
    stop(" A valid VCF needs at least chr, start.position, and end.position columns - check your file has these correctly named columns")
  }
  
  # using the dplyr tbl throughout this package
  library(dplyr)
  VCF=as.tbl(VCF)
  
  # add these columns if they exist
  if(strand %in% colnames(temp))
    VCF <- mutate(VCF, strand = temp[,strand])
  if(WT %in% colnames(temp))
    VCF <- mutate(VCF, WT = temp[,WT])
  if(MUT %in% colnames(temp))
    VCF <- mutate(VCF, MUT = temp[,MUT])
  if(sampleID %in% colnames(temp))
    VCF <- mutate(VCF, sampleID = temp[,sampleID])
  
  # if you want the other columns in the VCF then set "other = TRUE"
  if(other==TRUE)
  {
    other.cols<- setdiff(colnames(temp), c(chr,start.position,end.position, strand, WT, MUT, sampleID))
    VCF <- as.tbl(cbind(VCF, temp[,other]))
  }
  
  
  #short summary message and return a chr/pos sorted dplyr style tbl/data.frame
  cat(c("VCF contains: ", paste0(colnames(VCF), ", "), "\n\n"))
  VCF <- VCF %>%
    arrange(chr, start.position)
  
  class(VCF)<- c(class(VCF), "VCF")
  return(VCF)
  
}