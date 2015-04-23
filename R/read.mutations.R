#' Reads data from a MAF or VCF type file.
#' This is a wrapper for read.delim that reads mutation data from a tab delimited text file such as a MAF, VCF or other tab-delimited tsv file.
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
#' # A single example data file with multiple tumour samples comes with the package
#' (exF <- system.file("extdata", "1317049.tsv.gz", package="kater"))  
#' # a minimal example
#' PD4107a <- read.mutations(file=gzfile(exF), sampleID = Sample.Name)

read.mutations <- function(file,  strand=strand, sampleID = sampleID, pos=1:5, skip=0, header=T, ...){
  
  # Non Standard Evaluation: see http://adv-r.had.co.nz/Computing-on-the-language.html#capturing-expressions
  strand <- deparse(substitute(strand))
  sampleID <- deparse(substitute(sampleID))
  tempDF <- read.delim(file, as.is=T, ...)
 
  

  # check that at least the 5 main columns exist - fail if not.
  if(ncol(tempDF) >= 5 && nrow(tempDF) > 1)
  {
    # give standard header names
    DF = tempDF[, pos]
    colnames(DF)[1:5] <- c("chr", "start", "end", "ref", "mut")
  }
  else
  {
    stop("The data has insufficient columns or was somehow corrupted on import. Check the data format and function options")
  }
  
  # if strand or sampleID columns matchin the data were specified then add them
  if(any(colnames(tempDF) == strand))
    DF$strand <- tempDF[,strand]
  if(any(colnames(tempDF) == sampleID))
    DF$sampleID <- tempDF[, sampleID]
  
  # check the data looks of the correct type
  if(!is.integer(DF$start) | !is.integer(DF$end))
    stop("The start or end position data is not of type integer. Check the data format and function options")
  if(!is.integer(DF$start) | !is.integer(DF$end))
    stop("The start or end position data is not of type integer. Check the data format and function options")
  if(!(any(DF$chr==1) | any(tolower(DF$chr) == "chr1")))
     stop("The chromosome column looks odd, no 1's or chr1's. Check the data format and function options")
  
  
  # short message and return data as a dplyr style data.frame
  library(magrittr)
  DF=dplyr::as.tbl(DF)
  cat(c("imported mutations data.frame contains: ", paste0(colnames(DF), ", "), "\n\n"))
  DF <- DF %>%
    dplyr::arrange(chr, start)
  return(DF)
  
}