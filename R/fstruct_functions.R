# THIS FILE CONTAINS:
# structure.plot - generate a population structure plot using ggplot2
# fst_stat - generate fst, fstmax, and fst/fstmax for a Q matrix
# bootstrap - generate bootstrap Q matrices and related plots and statistics

# structure.plot -----------------------------------------------------------------
#' Plot a Q matrix using \code{\link[ggplot2]{ggplot2}}
#'
#' This function enables graphical visualization of a Q matrix, the default
#' output of population structure inference software programs such as
#' \href{https://web.stanford.edu/group/pritchardlab/structure.html}{STRUCTURE}
#' and \href{http://dalexander.github.io/admixture/index.html}{ADMIXTURE}. In
#' the output plot, each vertical bar represents a single individual's ancestry;
#' the height of each color in the bar corresponds to the individual membership
#' coefficients given by the Q matrix. Because this function produces a
#' \code{\link[ggplot2]{ggplot}} object, its output can be modified using
#' standard \code{\link[ggplot2]{ggplot2}} syntax. For a more comprehensive
#' population structure visualization program, see the program
#' \emph{\href{https://rosenberglab.stanford.edu/distruct.html}{distruct}}.
#'
#' @param Q A dataframe, matrix, or array representing a Q matrix. Each row
#'   represents an individual and the last \code{k} columns contain individual
#'   membership coefficients. The first few columns may contain information not
#'   relevant to this plot; their inclusion is optional. When restricted to the
#'   last \code{k} columns, the rows of this matrix should sum to approximately
#'   1.
#' @param k The number of ancestral clusters in the Q matrix. Each individual
#'   should have \code{k} membership coefficients. The default color scheme is
#'   "spectral" from \code{\link[RColorBrewer]{RColorBrewer}}, which tolerates
#'   \code{k} values of 11 or fewer.
#' @return A \code{ggplot} object describing a bar plot of membership
#'   coefficients from the Q matrix.
#' @examples
#' structure.plot(
#' # Make an example matrix of membership coefficients.
#' # Each row is an individual. Rows sum to 1.
#' Q = matrix(c(.4,.2,.4,
#'              .5,.3,.2,
#'              .5,.4,.1,
#'              .6,.1,.3,
#'              .6,.3,.1),
#'            nrow = 5,
#'            byrow = TRUE),
#' # How many ancestry coefficients per individual?
#' k = 3
#' ) +
#' # Below are example, optional modifications to the default plot
#'   ggtitle("Population A") +
#'   scale_fill_brewer("Blues") +
#'   xlab("Individuals")
#'
#' @export
structure.plot <- function(Q, k){
  # Check if Q matrix is in STRUCTURE/ADMIXTURE output form, or if it only contains columns of ancestry coefficients
  # If it is in STRUCTURE/ADMIXTURE output form, extract the ancestry coefficients--the last k columns.
  if(ncol(Q) > k){
    Q <- Q[(ncol(Q)-k+1):ncol(Q)]
  }

  # Generate the structure plot
  ggplot(data.frame(cbind(data.frame(Individuals =1:nrow(Q)),Q)) %>%
           pivot_longer(cols = 2:(ncol(Q)+1)),
         aes(fill=name, y=value, x=Individuals)) +
    geom_bar(position="stack", stat="identity", width=1) +
    scale_fill_brewer(palette = "Spectral") +
    theme_void() + ylab("") +
    theme(legend.position = "none")
}



#' @export
fst_stat <- function(Q){
  ## INPUT:
  # Q = matrix where each of row is a vector of ancestry coefficients
  ## OUTPUT:
  # A named list of statistics:
  # Fst = F_ST computed as if each individual were a population,
  # and each ancestral cluster were an allele
  # FstMax = The maximum value of Fst (for fixed frequency of the most frequent allele)
  # (result from Alcala, Rosenberg 2020(?))
  # ratio = the ratio of Fst/FstMax

  K = nrow(Q)
  p = colSums(Q) # yields vector of summed allele frequencies across populations
  sig1 = max(p)
  J = ceiling(1/sig1)
  sig1.frac = sig1 - floor(sig1)

  if(sig1 == K){
    FstMax = 0
    Fst = 0
    ratio = 0
  }else{if(sig1 <= 1){
    FstMax = ((K-1)*(1-sig1*(J-1)*(2-J*sig1)))/
      (K-1+sig1*(J-1)*(2-J*sig1))
  }else{
    FstMax = (K*(K-1)-2*(K-1)*sig1.frac*(1-sig1.frac)-
                floor(sig1)*(floor(sig1)-1)-2*sig1.frac*floor(sig1))/
      (K*(K-1)+2*sig1.frac*(1-sig1.frac)-floor(sig1)*(floor(sig1)-1)-2*sig1.frac*floor(sig1))
  }

    Fst =
      (sum(Q^2)/K-sum(colSums(Q/K)^2))/
      (1-sum(colSums(Q/K)^2))

    ratio = Fst/FstMax
  }

  return(list(Fst = Fst,
              FstMax = FstMax,
              ratio = ratio))
}


#' @export
bootstrap <- function(matrices, n_replicates){
  ## INPUT:
  # matrices = list of Q-matrices we seek to compare. Each of row is a vector of ancestry coefficients. Each matrix can be named (e.g. matrices = list(A=matrix(...), B=matrix(...)))
  # n_replicates = number of bootstrap replicates desired
  ## OUTPUT:
  # bootstrap_replicates = list of matrices of bootstrap replicates for each input matrix
  # statistics = dataframe; rows = bootstrap reps, columns = Matrix, Statistic, Value
  # histogram, CDF, boxplot of F_ST/F_ST max ratio bootstrap estimates

  n_matrix = length(matrices)
  # List of names of matrices
  names <- if(sum(!is.na(names(matrices)))){names(matrices)}else{1:n_matrix}

  ## For each matrix: ##
  for(m in 1:n_matrix){

    bs_list <- list()
    matrix = matrices[[m]]
    # Generate bootstrap data sets
    for(replicate in 1:n_replicates){
      bs_list[[replicate]] <- matrix[rdunif(n = nrow(matrix),
                                            a = 1, b = nrow(matrix)),]
    }

    # Compute statistics for these reps
    stats <- lapply(X = bs_list, FUN = fst_stat) %>% unlist %>%
      matrix(ncol = 3, byrow = TRUE) %>% data.frame()
    colnames(stats) <- c("Fst", "FstMax", "ratio")

    # Name this dataset, based on the name of the matrices in the list or the entry number
    assign(paste0("stats_", names[m]), stats, pos = sys.frame())
    assign(paste0("bootstrap_matrices_", names[m]), bs_list, pos = sys.frame())
  }

  # Make a long dataset with all matrices' statistics:
  all_stats <- cbind(Matrix =
                       names %>%
                       lapply(function(name)rep(name, n_replicates)) %>%
                       unlist,
                     do.call(rbind,
                             paste0("stats_", names) %>% lapply(get))) %>%
    pivot_longer(data = ., cols = 2:4, names_to = "Statistic")

  all_stats$Matrix <- factor(all_stats$Matrix, levels = unique(all_stats$Matrix))

  compare_ratio_hist <-
    ggplot(data = all_stats %>%
             filter(Statistic == "ratio")) +
    geom_histogram(aes(value, fill = Matrix)) +
    xlab("Ratio of Fst to Max. Fst")

  compare_ratio_ecdf <- ggplot(data = all_stats %>%
                                 filter(Statistic == "ratio")) +
    stat_ecdf(aes(value, color = Matrix)) +
    xlab("Ratio of Fst to Max. Fst") + ylab("Cumulative Probability")

  compare_ratio_boxplot <- ggplot(data = all_stats %>%
                                    filter(Statistic == "ratio"),
                                  aes(x = Matrix, y = value, color = Matrix)) +
    geom_boxplot() +
    ylab("Ratio of Fst to Max. Fst")

  kruskal_wallace <- kruskal.test(value ~ Matrix,
                                  data= filter(all_stats, Statistic == "ratio"))

  return(list(statistics = all_stats,
              bootstrap_replicates = lapply(paste0("bootstrap_matrices_",
                                                   names),
                                            get),
              histogram = compare_ratio_hist,
              cdf = compare_ratio_ecdf,
              boxplot = compare_ratio_boxplot,
              KW_test_result = kruskal_wallace))
}
