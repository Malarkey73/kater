#' Count of dinucleotide mutations (e.g. CC>AA, TT>CC)#'
#' This function counts the different types of dinucleotide substitutions (e.g. CC to AA). As with the simple pointMutations function it also counts the reverse complement... so CCtoAA is really the count of both CCtoAA and GGtoTT. Most NGS technologies (not all) are not stranded so you cannot tell if a CCtoAA on positive strand was caused by a GGtoTT on the negative starnd (or vice versa).
#' @param DF a delimited file with columns ref and mut, containing values including "A", "T", "G" or "C" (a data.frame)
#' @keywords monte carlo, mutants, substitutions, point mutations
#' @examples
#' res <- dinucmut.count(DF)

dinucMutations<- function(DF)
{
  
  # filter down to second substitutions (dinucleotides, dn) same as the one before
  dn <- DF %>%
    dplyr::arrange(chr, start) %>%
    dplyr::mutate(im.lag = start-lag(start)) %>%
    dplyr::filter(im.lag ==1 & ref == lag(ref) & mut ==lag(mut)) 
  
  # this function can be used alone but here it counts the dinucs - crucially ignoring any indels         
  if(nrow(dn)>0)
  {
    res<- pointMutations(dn)
  }else
  {
    res= rep(0,6)
  }
  
  names(res)=c("CCtoAA", "CCtoGG", "CCtoTT", "TTtoAA", "TTtoCC", "TTtoGG")
  return(res)
  
}