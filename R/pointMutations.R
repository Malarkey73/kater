#' Count of different point mutations(e.g. C>A, T>C)#'
#' Ignoring Indels (e.g. ATGCGC > A) this function simply counts the different types of point mutation. Note that the CtoA category  contains both CA and the reverse complement - as in most sequencing experiments we cannot distinguish which strand was actually mutated.
#' @param DF a delimited file with columns "ref" and "mut", containing values including "A", "T", "G" or "C" (a data.frame)
#' @keywords monte carlo, mutations, substitutions, point mutations
#' @examples
#' res <- pointmuts(DF)
#' 
pointmuts= function(DF)
{
  
  #C>A is the same as G>T on the -strand 
  CtoA<- DF %>% summarise(CA=sum((ref=="C" & mut=="A")|(ref=="G" & mut=="T")))
  #C>G is the same as G>C on the -strand 
  CtoG<- DF %>% summarise(CG=sum((ref=="C" & mut=="G")|(ref=="G" & mut=="C")))
  #C>T is the same as G>A on the -strand 
  CtoT<- DF %>% summarise(CT=sum((ref=="C" & mut=="T")|(ref=="G" & mut=="A")))
  
  #T>A is the same as A>T on the -strand 
  TtoA<- DF %>% summarise(TA=sum((ref=="T" & mut=="A")|(ref=="A" & mut=="T")))
  #T>C is the same as A>G on the -strand 
  TtoC<- DF %>% summarise(TC=sum((ref=="T" & mut=="C")|(ref=="A" & mut=="G")))
  #T>G is the same as A>C on the -strand 
  TtoG<- DF %>% summarise(TG=sum((ref=="T" & mut=="G")|(ref=="A" & mut=="C")))
  
  message("NB A count of  e.g. CtoA mutations also includes reverse complement GtoT mutations.")
  return(cbind(CtoA,CtoG,CtoT,TtoA,TtoC,TtoG))
}