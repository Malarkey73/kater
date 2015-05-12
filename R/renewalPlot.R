#' A renewal plot to detect non homogenous Poisson processes
#'
#' This idea is borrowed from the analysis of spike train dat described in the STAR package - and reapplied here to a completely different problem. The basic idea is to use the
#' rank order of the intermutation distances OJ plotted against OJ+1 (the lag). This uses the plot space better than just plotting intermutation vs intermutaion lag.
#' Homogenous Poisson processes fill the square  space randomly, deviations from homogeneity leave gaps.
#' @param VCF
#' @keywords rainfall, mutants, kataegesis
#' @examples
#' # a minimal example
#' DF<-read.mutations(file="1317049.tsv", chr="Chromosome", start.position="Genome.start", end.position="Genome.stop")
#' renewal.plot(DF)

renewalPlot<-function(DF)
{
 

REN<- DF %>%
  dplyr::group_by(chr) %>%
  dplyr::mutate(intermutation=start-lag(start)) %>%
  dplyr::mutate(OJ=rank(intermutation)) %>%
  dplyr::mutate(OJ.l=lag(OJ))

library(ggplot2)
theme_set(theme_bw())
QP <- qplot(OJ, OJ.l, data=REN)+facet_wrap(~chr, scales = "free")
return(QP)
}