#' Create a table of types of point mutation
#' A type of point mutation is TAA > TGA (triplet = "TGA", mutation ="A" )
#' @param triplets E.g. TCG, ACT, AAA etc. These might be recovered using the expandSeq function on a DF file.
#' @param mutations The mutation at the centre of the triplet e.g. TCG > TGG is mutation G
#' @keywords Poisson, mutants, kataegesis
#' @export
#' @examples
#' tripletmutcount(DF, triplets = triplets, mutations = MUT)
#'

tripletSignatures= function(DF)
  
{
  
  #prep
  all.mutations <- c("A","C", "G", "T")
  all.triplets <- paste0(rep(all.mutations, times=16), rep(all.mutations, each=4), rep(all.mutations, each=16))
  
  # This is awkward but I need to make a table of all possible mutation types and merge with the actual mutation table
  # remember a mutation canot be AAA > A or CTT > T - so these are removed now.
  mutation.allposs <- dplyr::as.tbl(data.frame(expandref=rep(all.triplets,each=4), mut=all.mutations)) %>%
    dplyr::filter(substr(expandref,2,2) != mut) %>%
    dplyr::arrange(expandref, mut)
  
  # this is the actual mutation type counts
  mutation.actual <- DF %>%
    dplyr::filter(nchar(mut)==1, nchar(expandref)==3) %>%
    dplyr::group_by(expandref,mut) %>% 
    dplyr::summarise(count = n())
    
  # this is just in case there are missing mutations in the actual data, to get 0 counts for them rather than blank rows.
  mutation.join = dplyr::left_join(mutation.allposs, mutation.actual) %>%
    dplyr::mutate(count = ifelse(is.na(count), 0, count))
  
  # each mutation has a single reverse complement which (with most types of suequencing)
  # you cannot distinguish ... so really there are 96 unique point mutations NOT 192
  mutation.join <- mutation.join %>%
    dplyr::mutate(rctriplets = as.character(reverseComplement(DNAStringSet(expandref)))) %>%
    dplyr::mutate(rcmutations = as.character(reverseComplement(DNAStringSet(mut))))
  
  # this code adds the counts from the reverse complement then removes the 96 redundant rows.
  matches <- with(mutation.join, match(paste0(expandref,mut), paste0(rctriplets,rcmutations)))
  res <- mutation.join %>%
    dplyr::mutate(count.both= count + count[matches]) %>%
    dplyr::filter(substr(expandref,2,2) %in% c("C", "T"))
  
  
  return(res)
}