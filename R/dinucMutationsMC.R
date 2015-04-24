#' Dinucleotide Monte Carlo Simulation'
#' In Nik-Zainal et al. 2012 paper "Mutational Processes Molding the Genomes of 21 Breast Cancers", they note a higher rate of dinucleotide mutations than expected by chance. They describe this as estimated by Monte Carlo simulations. This function is my version of such a simulation. I scramble the existing data n=100 times, assuming the same number of mutations on each chromosoem Poisson distributed across the observed range e.g. between the first to last actual mutation. The function then calculates permuted counts (Expected Counts :E) for each dinucleotide type each time. The fraction of E > O (Observed Counts: O) is a Monte Carlo estimate of P (the null hypothesis).   
#' @param DF a delimited file with columns "WT" and "MUT", containing values including "A", "T", "G" or "C" (a data.frame)
#' @param n.mc the number of Monte Carlo simulations (permutations) of the expected dincleotide rate E. 
#' @keywords monte carlo, mutants, substitutions, point mutations
#' @examples
#' res <- dinucmut.mc(DF, n=100)

dinucMutationsMC <- function(DF, n.mc=100)
{
  
  # the observed dinucleotide counts
  Observed <- dinucMutations(DF)
  
  # n=100 rows by default, a column per dinuc type, same as Observed 
  Expected=as.data.frame(matrix(NA, n.mc, 6))
  names(Expected)=c("CCtoAA", "CCtoGG", "CCtoTT", "TTtoAA", "TTtoCC", "TTtoGG")
  
  # recording the observed chromosome ranges
  DF.summ<- DF %>%
    dplyr::group_by(chr) %>%    
    dplyr::summarise(first.position=first(start), last.position=last(end), n=n())
  
  
  
  for(i in 1:n.mc)
  {
    
    # This essentially scrambles the same mutations of same types on same chromosomes but moves them to random positions
    # within the observed chromosome ranges
    
    DF.scramble<- DF %>%
      dplyr::mutate(start = unlist(apply(DF.summ, 1, function(x) floor(runif(x[4], x[2], x[3]))))) %>%
      dplyr::arrange(chr, start)
    
    Expected[i,] <- dinucMutations(DF.scramble)
    
    if((i %% 10) == 0)
      message(i, " scrambles")
  }
  
  P.mc = rowMeans(apply(Expected,1, function(x) x > Observed))
  names(P.mc) = names(Expected)
  return(P.mc)
}