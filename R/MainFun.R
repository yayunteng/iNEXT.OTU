#' function to calculate UniFrac distance based on dissimilarity measure
#'
#' \code{iNEXTUniFrac}: function to calculate UniFrac distance based on Sørensen- and Jaccard-type dissimilarity measure
#'
#' @param data OTU data can be input as a \code{matrix/data.frame} (species by assemblages), or a \code{list} of \code{matrices/data.frames}, each matrix represents species-by-assemblages abundance matrix.\cr
#' @param q a numerical vector specifying the diversity orders. Default is c(0, 1, 2).
#' @param base Sample-sized-based rarefaction and extrapolation for gamma and alpha diversity (\code{base = "size"}) or coverage-based rarefaction and extrapolation for gamma, alpha and beta diversity (\code{base = "coverage"}). Default is \code{base = "coverage"}.
#' @param level A numerical vector specifying the particular value of sample coverage (between 0 and 1 when \code{base = "coverage"}) or sample size (\code{base = "size"}). \code{level = 1} (\code{base = "coverage"}) means complete coverage (the corresponding diversity represents asymptotic diversity).\cr
#' If \code{base = "size"} and \code{level = NULL}, then this function computes the gamma and alpha diversity estimates up to double the reference sample size. \cr
#' If \code{base = "coverage"} and \code{level = NULL}, then this function computes the gamma and alpha diversity estimates up to one (for \code{q = 1, 2}) or up to the coverage of double the reference sample size (for \code{q = 0});
#' the corresponding beta diversity is computed up to the same maximum coverage as the alpha diversity.
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
#' @return If \code{base = "coverage"}, return a list of seven data frames with three coverage-based diversity (gamma, alpha, and beta) and four types dissimilarity measure for UniFrac distance. If \code{base = "size"}, return a list of two data frames with two diversity (gamma and alpha).
#'
#' @examples
#'
#' data("tongue_cheek")
#' data("tongue_cheek_tree")
#' output <- iNEXTUniFrac(tongue_cheek, q=c(0,1,2),
#'                        level = seq(0.5, 1, 0.05), nboot = 10,
#'                        conf = 0.95, PDtree = tongue_cheek_tree, PDreftime = NULL)
#'
#' @export
iNEXTUniFrac = function(data, q=c(0,1,2), base = "coverage", level = NULL, nboot = 10, conf = 0.95, PDtree = NULL, PDreftime = NULL){
  out = iNEXTbeta3D(data, diversity = "PD", q=q, datatype = "abundance", base = base, level = level, nboot = nboot, conf = conf, PDtree = PDtree,  PDreftime = PDreftime)
  out
  # UniFrac_out = list()
  # if(class(data)=="list"){
  #   for(i in 1:length(data)){
  #     UniFrac_out[[i]] = list(C = out[[i]]$C, U = out[[i]]$U)
  #   }
  #   if(is.NULL(names(data))){
  #     names(UniFrac_out) = paste0("Region_", 1:length(data))
  #   }else{
  #     names(UniFrac_out) = names(data)
  #   }
  #
  # }else{
  #   UniFrac_out[[1]] = list(C = out[[1]]$C, U = out[[1]]$U)
  #   names(UniFrac_out) = "Region_1"
  # }
  #UniFrac_out
}

# tongue_cheek %>% colnames()
# model1 <- data.frame(Cheek = tongue_cheek[,1],
#                      Tongue = tongue_cheek[,31],
#                      row.names = rownames(tongue_cheek))
# model1 <- model1 %>% filter(Cheek>0 | Tongue>0)
# model1 %>% dim
# output <- iNEXTbeta3D(model1, diversity = "PD", q=c(0,1,2), level = seq(0.5, 1, 0.05), nboot = 0, conf = 0.95, PDtree = tree, PDreftime = NULL)

#' ggplot2 extension for an iNEXT.UniFrac object
#'
#' \code{ggiNEXTUniFrac}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXTUniFrac}} object to plot coverage-based rarefaction/extrapolation curves for UniFrac.
#'
#' @param output the output from iNEXTUniFrac
#' @param type (required only when \code{base = "coverage"}), selection of plot type : \cr
#' \code{type = 'B'} for plotting the gamma, alpha, and beta diversity ;  \cr
#' \code{type = 'D'} for plotting 4 turnover dissimilarities.
#' @param scale Are scales shared across all facets (the default, \code{"fixed"}), or do they vary across rows (\code{"free_x"}), columns (\code{"free_y"}), or both rows and columns (\code{"free"})?
#' @param transp a value between 0 and 1 for controlling transparency. \code{transp = 0} is completely transparent, default is 0.4.
#'
#' @return a figure for two types of UniFrac based on dissimilarity measure.
#'
#' @examples
#'
#' data("tongue_cheek")
#' data("tongue_cheek_tree")
#' output <- iNEXTUniFrac(tongue_cheek, q=c(0,1,2), nboot = 0, PDtree = tongue_cheek_tree)
#' ggiNEXTUniFrac(output, scale = 'free', transp = 0.4)
#'
#' @export
ggiNEXTUniFrac = function(output, type = "B", scale = "fixed", transp = 0.4){
  ggiNEXTbeta3D(output, type = type, scale = scale, transp = transp)

  # cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  # ylab = "UniFrac distance"
  #
  #   C = lapply(output, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
  #   U = lapply(output, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
  #   C = C %>% filter(Method != 'Observed')
  #   U = U %>% filter(Method != 'Observed')
  #   C[C == 'Observed_alpha'] = U[U == 'Observed_alpha'] = 'Observed'
  #
  #   df = rbind(C, U)
  #   for (i in unique(C$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
  #   df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN"))
  #
  #   id_obs = which(df$Method == 'Observed')
  #
  #   for (i in 1:length(id_obs)) {
  #     new = df[id_obs[i],]
  #     new$SC = new$SC - 0.0001
  #     new$Method = 'Rarefaction'
  #
  #     newe = df[id_obs[i],]
  #     newe$SC = newe$SC + 0.0001
  #     newe$Method = 'Extrapolation'
  #
  #     df = rbind(df, new, newe)
  #   }
  #
  #   lty = c(Rarefaction = "solid", Extrapolation = "dashed")
  #   df$Method = factor(df$Method, levels = c('Rarefaction', 'Extrapolation', 'Observed'))
  #
  #   double_size = unique(df[df$Method == "Observed",]$Size)*2
  #   double_extrapolation = df %>% filter(Method == "Extrapolation" & round(Size) %in% double_size)
  #
  #   ggplot(data = df, aes(x = SC, y = Estimate, col = Region)) +
  #     geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha = transp) +
  #     geom_line(data = subset(df, Method != 'Observed'), aes(linetype = Method), size=1.1) + scale_linetype_manual(values = lty) +
  #     # geom_line(lty=2) +
  #     geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 2) +
  #     geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 2, stroke = 1.2)+
  #     geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 2) +
  #     geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 2, stroke = 1.2) +
  #     scale_colour_manual(values = cbPalette) +
  #     scale_fill_manual(values = cbPalette) +
  #     facet_grid(div_type ~ Order.q, scales = scale) +
  #     theme_bw() +
  #     theme(legend.position = "bottom", legend.title = element_blank()) +
  #     labs(x = 'Sample coverage', y = ylab)
}

