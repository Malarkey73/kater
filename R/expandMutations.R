#' Get surrounding sequence context
#' uses the Bioconductor BSGenome framework to find
#' @importFrom magrittr "%>%"
#' @param DF tbl or data.frame with columns containing at least "chr", "start", "end", "ref", and "mut" 
#' @param expand The number of upstream and downstream bases to get
#' @param genomeseq The Bioconductor genome annotation package to retireve sequence information
#' @keywords Poisson, mutants, kataegesis
#' @export
#' @examples
#' # This will add a column of the triplet context surrounding point mutations
#' expandseq(DF, expand=1)

expandMutations<- function(DF, expand=1, genome = BSgenome.Hsapiens.UCSC.hg19)
{

  #NSE
  genome=deparse(substitute(genome))
  
  if(grep("UCSC", genome))
    {
    # This is the format of the UCSC.hg19 data "chr19" NOT "19", "chrX" not "chr23"
      DF <- DF %>%
        dplyr::mutate(adj.chr = ifelse(chr %in% 1:24, yes=paste0("chr", chr), no=chr)) %>%
        dplyr::mutate(adj.chr = ifelse(chr==23, yes="chrX", no=adj.chr)) %>%
        dplyr::mutate(adj.chr = ifelse(chr==24, yes="chrY", no=adj.chr))
        
    }
    else
    {
      DF <- DF %>%
        dplyr::mutate(adj.chr = chr)
    }
  
  # select the right library
  genome <- BSgenome::getBSgenome(genome)
  message("This may take a few seconds.")
  # if that library is human then use this to access the sequence data.
  gs <- BSgenome::getSeq(genome, DF$adj.chr, start= DF$start - 1, end= DF$end + 1)
  
  DF<- DF %>% 
    dplyr::mutate(expandref = as.character(gs)) %>%
    dplyr::select(-adj.chr)
  
  return(DF)
}