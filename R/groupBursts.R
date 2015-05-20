#' Grouping Bursts of Mutations
#'
#' This function aims to group consecutive bursts of mutations into numbered groups (0 is not a burst).
#' These groups can then be used for simple summary statistics functions describing the spatial patterns
#' of mutation within each sample. The definition of a burst is very simple: at least 3 consecutive 
#' mutations all closer together than a certain threshold (gap default 250 bp). Bursts end when a gap of 
#' more than the threshold opens.
#' @param DF a tbl with columns containing at least "chr", "start", and "end" (a data.frame)
#' @param gap below which bursts are detected 
#' @keywords rainfall, mutants, kataegesis
#' @export
#' @examples
#' groupBurst(DF, 200)


groupBursts<- function(DF, max.threshold=250)
{
  kat<- DF %>%
    group_by(chr) %>%
    # next position minus current position
    mutate(leading = lead(start)-start,
           # TRUE/FALSE less than a certain value (bmth = below max threshold)
           bmth = leading < max.threshold,
           # remove singletons
           bmth = ifelse(lag(bmth)==FALSE & lead(bmth)==FALSE , FALSE, bmth)
    )
  
  bmth<-kat$bmth
  
  burst.gr<-0
  vec=rep(NA, length(bmth))
  last=F
  
  for(i in 1:length(bmth))
  {
    if(bmth[i]==FALSE | is.na(bmth[i]))
    {
      vec[i]<-0
      last<- FALSE
    }
    
    if(bmth[i]==TRUE & last==FALSE & !is.na(bmth[i]))
    {
      burst.gr<- burst.gr +1
      last<-TRUE
      vec[i]<- burst.gr
    }
    
    if(bmth[i]==TRUE & last==TRUE & !is.na(bmth[i]))
    {
      vec[i]<- burst.gr
    }
    
  }
  return(as.tbl(data.frame(DF, bursts=as.factor(vec))))
  
}