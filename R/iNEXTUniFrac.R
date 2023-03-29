#' function to calculate UniFrac based on dissimilarity measure
#'
#' \code{iNEXTUniFrac}: function to calculate UniFrac based on Sørensen- and Jaccard-type dissimilarity measure
#'
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a \code{matrix/data.frame} (species by assemblages), or a \code{list} of \code{matrices/data.frames}, each matrix represents species-by-assemblages abundance matrix.\cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a \code{list} (a region) with several \code{lists} (assemblages) of \code{matrices/data.frames}, each matrix represents species-by-sampling units.
#' @param q a numerical vector specifying the diversity orders. Default is c(0, 1, 2).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being \code{0} (non-detection) or \code{1} (detection).
#' @param level A numerical vector specifying the particular value of sample coverage (between 0 and 1). \code{level = 1} means complete coverage (the corresponding UniFrac represents asymptotic UniFrac).\cr
#' If \code{level = NULL}, this function computes the gamma and alpha diversity estimates up to one (for \code{q>0}) or up to the coverage of double the reference sample size (for \code{q = 0});
#' the corresponding beta diversity and UniFrac is computed up to the same maximum coverage as the alpha diversity.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter \code{0} to skip the bootstrap procedures. Default is \code{10}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param PDtree a phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param PDreftime  a numerical value specifying reference time for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).
#'
#' @import tidyverse
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import ape
#' @import phytools
#' @import phyclust
#' @import tidytree
#' @import RColorBrewer
#' @import iNEXT.3D
#' @import future.apply
#' @import ade4
#' @import tidyr
#' @import tibble
#' @import iNEXT.beta3D
#'
#' @return a list of two data frames with two types dissimilarity measure for UniFrac.
#'
#' @examples
#'
#' data("tongue_cheek")
#' data("tongue_cheek_tree")
#' output <- iNEXTUniFrac(tongue_cheek, q=c(0,1,2), level = seq(0.5, 1, 0.05), nboot = 10, conf = 0.95, PDtree = tongue_cheek_tree, PDreftime = NULL)
#'
#' @export
iNEXTUniFrac = function(data, q=c(0,1,2), datatype = "abundance", level = NULL, nboot = 10, conf = 0.95, PDtree = NULL, PDreftime = NULL){
  out = iNEXTbeta3D(data, diversity = "PD", q=q, datatype = datatype, level = level, nboot = nboot, conf = conf, PDtree = PDtree,  PDreftime = PDreftime)

  UniFrac_out = list()
  if(class(data) == "list"){
    for(i in 1:length(data)){
      UniFrac_out[[i]] = list(C = out[[i]]$C, U = out[[i]]$U)
    }
  }else{
    UniFrac_out[[1]] = list(C = out[[1]]$C, U = out[[1]]$U)
  }
  names(UniFrac_out) = names(data)
}

# tongue_cheek %>% colnames()
# model1 <- data.frame(Cheek = tongue_cheek[,1],
#                      Tongue = tongue_cheek[,31],
#                      row.names = rownames(tongue_cheek))
# model1 <- model1 %>% filter(Cheek>0 | Tongue>0)
# model1 %>% dim
# output <- iNEXTbeta3D(model1, diversity = "PD", q=c(0,1,2), level = seq(0.5, 1, 0.05), nboot = 0, conf = 0.95, PDtree = tree, PDreftime = NULL)



