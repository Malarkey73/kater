statBursts <- function(DF, bursts="bursts", intra.burst=T)
{
  kat<- DF %>%
    group_by(bursts, chr) %>%
    mutate(intermutation=lead(start)-start) %>%
    summarise(burst.start= min(start),
              burst.end=max(start),
              burst.length= burst.end-burst.start,
              n=n(),
              min.gap = min(intermutation, na.rm=T),
              q25.gap = quantile(intermutation, 0.25, na.rm=T),
              median.gap = median(intermutation,na.rm=T),
              mean.gap = mean(intermutation,na.rm=T),
              q75.gap = quantile(intermutation, 0.25, na.rm=T),
              max.gap = max(intermutation,na.rm=T))
  
  if(intra.burst==TRUE){
    kat <- kat %>%
      group_by(chr) %>%
      mutate(next.burst= lead(burst.start)-burst.end) %>%
      filter(bursts != 0)
    
  }
  
  return(kat)
  
}