#' Rainfall Plot
#'
#' A rainfall plot show the spatial clustering of mutations .. usually point mutations
#' @param DF a tbl with columns containing at least "chr", "start", and "end" (a data.frame)
#' @param chrN the chromosome to plot(character)
#' @param lower the lower range of the plot (numeric)
#' @param upper the upper range of the plot x-axis (numeric)
#' @param colour a column to specify point colours (character)
#' @param facets a column to specify plot faceting (character)
#' @param upper the upper range of the plot x-axis (numeric)
#' @keywords rainfall, mutants, kataegesis
#' @export
#' @examples
#' ## before using the rainfall plot you might want to filter the VCF to contain mutations of a certain type.
#' ## For instance this is filtering for APOBEC3A/B type point mutations
#' DF <- filter(DF, (ref=="C" & mut=="T")|(ref=="C" & mut=="G")|(ref=="G" & mut=="A")|(ref=="G" & mut=="C"))
#' ## You might also add a column that can be used for shape or colour, e.g.
#' DF <- mutate(DF, mutation_type = ifelse((ref=="C" & mut=="T")|(ref=="G" & mut=="A"), "transition", "transversion"))
#' ## Or guess  the strand of the mutation actual
#' DF <- mutate(DF, strand = ifelse((ref=="C" & mut=="T")|(ref=="C" & mut=="G"), "positive", "negative"))
#' rainfall(DF, chrN = 6, lower = 1.24e8, upper = 1.4e8, facets="strand", colour="mutation_type")
#'

rainfallPlot <- function(DF, chrN=1, lower, upper, colour=NULL, facets=NULL, shape=NULL, bursts=F, s=3, gamma=3, vline=NULL, ...){
  library(dplyr)
  library(ggplot2)
  theme_set(theme_bw())
  
  # filter by chr and range, sort by position, add intermutation distance
  RF <- DF %>%
    dplyr::filter(chr == chrN, start > lower, end < upper) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(intermutation=lead(start)-start) %>%
    dplyr::filter(intermutation > 1)
  
  
  QP <- qplot(x=start, y=intermutation, data=RF, geom="point", alpha=I(0.8), size =I(3))+
    scale_y_log10()+
    theme(legend.position="bottom", strip.background=element_rect(fill="white")) +
    scale_x_continuous()
  
  # add facets for e.g. strand or different samples
  if(!is.null(facets))
    QP <- QP + facet_wrap(as.formula(sprintf("~ %s", facets)))
  
  # add a colour e.g. for mutation types
  if(!is.null(colour))
    QP <- QP + aes_string(colour=colour) + scale_colour_brewer(palette="Set1") +
    
    # add a shape e.g. for transition transversion
    if(!is.null(shape))
      QP <- QP + aes_string(shape=shape)
  
  if(!is.null(vline))
    QP<- QP + geom_vline(x=vline, linetype="longdash")
  
  # add kleinberg burst detection
  if(bursts==TRUE)
  {
    library(bursts)
    kl= kleinberg(RF$start,  s=s, gamma=gamma)
    kl<- kl[kl$level>1,]
    QP <- QP + geom_segment(aes(x = start, xend = end, y=level, yend=level), colour= I(3), size=I(6), data=kl)
  }
  
  
  return(QP)
}