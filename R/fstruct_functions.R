# THIS FILE CONTAINS:
# Q_plot - generate a population structure plot using ggplot2
# Q_stat - generate fst, fstmax, and Fst/FstMax for a Q matrix
# Q_bootstrap - generate bootstrap Q matrices and related plots and statistics
# Q_simulate - simulate Q matrices using the Dirichlet distribution

# Q_plot -----------------------------------------------------------------
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
#'   represents an individual and the last \code{K} columns contain individual
#'   membership coefficients. The first few columns may contain information not
#'   relevant to this plot; their inclusion is optional. When restricted to the
#'   last \code{K} columns, the rows of this matrix should sum to approximately
#'   1.
#' @param K The number of ancestral clusters in the Q matrix. Each individual
#'   should have \code{K} membership coefficients. The default color scheme is
#'   "spectral" from \code{\link[RColorBrewer]{RColorBrewer}}, which tolerates
#'   \code{K} values of 11 or fewer.
#' @return A \code{ggplot} object describing a bar plot of membership
#'   coefficients from the Q matrix.
#' @examples
#' Q_plot(
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
#' K = 3
#' ) +
#' # Below are example, optional modifications to the default plot
#'   ggplot2::ggtitle("Population A") +
#'   ggplot2::scale_fill_brewer("Blues") +
#'   ggplot2::xlab("Individuals")
#'
#'@importFrom dplyr %>%
#' @export
Q_plot <- function(Q, K){
  # Check if Q matrix is in STRUCTURE/ADMIXTURE output form, or if it only contains columns of ancestry coefficients
  # If it is in STRUCTURE/ADMIXTURE output form, extract the ancestry coefficients--the last K columns.
  if(ncol(Q) > K){
    Q <- Q[,(ncol(Q)-K+1):ncol(Q)]
  }

  # Generate the structure plot
  ggplot2::ggplot(data.frame(cbind(data.frame(Individuals =1:nrow(Q)),Q)) %>%
           tidyr::pivot_longer(cols = 2:(ncol(Q)+1)),
           ggplot2::aes(fill=name, y=value, x=Individuals)) +
    ggplot2::geom_bar(position="stack", stat="identity", width=1) +
    ggplot2::scale_fill_brewer(palette = "Spectral") +
    ggplot2::theme_void() +
    ggplot2::ylab("") +
    ggplot2::theme(legend.position = "none")
}



# Q_stat -----------------------------------------------------------------
#' Compute Fst, FstMax, and the ratio Fst/FstMax for a Q matrix
#'
#' This function computes a statistical measure of ancestry variability, Fst/FstMax, for a Q matrix, the default output of population structure inference software programs such as \href{https://web.stanford.edu/group/pritchardlab/structure.html}{STRUCTURE} and \href{http://dalexander.github.io/admixture/index.html}{ADMIXTURE}. The function returns a named list containing the ratio Fst/FstMax as well as the values of Fst and FstMax.
#'
#'  Fst/FstMax is a statistic which takes a value of 0 when every individual in a population has identical ancestry, and a value of 1 when the ancestry is maximally variable (see *our paper* for more details). It is based on the population differentiation statistic Fst which, in its traditional application, is used to measure variability in allele frequencies
#'
#' @param Q A dataframe, matrix, or array representing a Q matrix. Each row
#'   represents an individual and the last \code{K} columns contain individual
#'   membership coefficients. The first few columns may contain information not
#'   relevant to this plot; their inclusion is optional. When restricted to the
#'   last \code{K} columns, the rows of this matrix should sum to approximately
#'   1.
#' @param K The number of ancestral clusters in the Q matrix. Each individual
#'   should have \code{K} membership coefficients.
#' @return A named list of containing the following entries:
#' \itemize{
#' \item  \code{Fst}: Fst computed as if each individual were a population, and each ancestral cluster were an allele.
#' \item \code{FstMax}: The maximum value of Fst (for fixed frequency of the most frequent allele, or, in the analogy, the membership of the most prevalent ancestral cluster).
#' \item \code{ratio}: The ratio of Fst/FstMax. We recommend that this statistic is used to quantify ancestry variability and compare the variability of two or more Q matrices.
#' }
#' @examples
#' #' Q_stat(
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
#' K = 3
#' )
#'
#' @export
Q_stat <- function(Q, K){
  # Check if Q matrix is in STRUCTURE/ADMIXTURE output form, or if it only contains columns of ancestry coefficients
  # If it is in STRUCTURE/ADMIXTURE output form, extract the ancestry coefficients--the last K columns.
  if(ncol(Q) > K){
    Q <- Q[,(ncol(Q)-K+1):ncol(Q)]
  }

  I=nrow(Q) # Here, I is the number of individuals (number of subpopulations)
  p = colSums(Q) # yields vector of summed allele frequencies across populations
  sig1 = max(p)
  J = ceiling(1/sig1)
  sig1.frac = sig1 - floor(sig1)

  if(sig1 == I){
    FstMax = 0
    Fst = 0
    ratio = 0
  }else{if(sig1 <= 1){
    FstMax = ((I-1)*(1-sig1*(J-1)*(2-J*sig1)))/
      (I-1+sig1*(J-1)*(2-J*sig1))
  }else{
    FstMax = (I*(I-1)-2*(I-1)*sig1.frac*(1-sig1.frac)-
                floor(sig1)*(floor(sig1)-1)-2*sig1.frac*floor(sig1))/
      (I*(I-1)+2*sig1.frac*(1-sig1.frac)-floor(sig1)*(floor(sig1)-1)-2*sig1.frac*floor(sig1))
  }

    Fst =
      (sum(Q^2)/I-sum(colSums(Q/I)^2))/
      (1-sum(colSums(Q/I)^2))

    ratio = Fst/FstMax
  }

  return(list(Fst = Fst,
              FstMax = FstMax,
              ratio = ratio))
}


# Q_bootstrap ----------------------------------------
#' Generate and analyze bootstrap replicates of one or more Q matrices
#'
#' Generates bootstrap replicate Q matrices, computes Fst/FstMax for each bootstrap replicate, produces several plots of the bootstrap distributions of Fst/FstMax for each provided Q matrix, and runs two statistical tests comparing these bootstrap distributions. The tests comparing bootstrap distributions of Fst/FstMax facilitate statistical comparison of the variability in each of multiple Q matrices.
#'
#' @param matrices A dataframe, matrix, or array representing a Q matrix, or a (possibly named) list of arbitrarily many such objects. For each Q matrix, matrix rows represents an individual and the last \code{K} columns contain individual membership coefficients (when restricted to the last \code{K} columns, the rows should sum to approximately 1). If the list of matrices is not named, the matrices are numbered and the numbers are used in place of names.
#' @param n_replicates The number of bootstrap replicate matrices to generate for each provided Q matrix.
#' @param K The number of ancestral clusters in the Q matrix. Each individual should have \code{K} membership coefficients.
#'
#' @return A named list of containing the following entries:
#' \itemize{
#' \item \code{bootstrap_replicates}: A named list of lists. Each element is named for a Q matrix provided in \code{matrices} and contains a list of \code{n_replicates} bootstrap replicates of the provided matrix. E.g., if \code{n_replicates = 100} and the first Q matrix in \code{matrices} is named \code{A}, then the first element of \code{bootstrap_replicates}, \code{bootstrap_replicates$bootstrap_matrices_A}, is itself a list of 100 matrices, each representing a bootstrap replicate of matrix A.
#' \item \code{statistics}: A dataframe containing the output of \code{Q_stat}: \code{Fst}, \code{FstMax}, and \code{ratio} (Fst/FstMax), computed for each bootstrap replicate matrix in \code{bootstrap_replicates}. The ratio Fst/FstMax quantifies the variability of each  Q matrix. The first column, titled \code{Matrix}, is a factor indicating which provided Q matrix the row corresponds to (the matrix name if \code{matrices} is a named list, or a number otherwise). The row names are of the form \code{stats_matrix.replicate} where \code{matrix} is the name of one of the provided Q matrices (or the entry number if the list elements were not named) and replicate is the number of bootstrap replicate (rep takes values from 1 to \code{n_replicates}).
#' \item \code{plot_boxplot}: A ggplot2 box plot depicting the bootstrap distribution of Fst/FstMax for each matrix in \code{matrices}.
#' \item \code{plot_violin}: A ggplot2 violin plot depicting the bootstrap distribution of Fst/FstMax for each matrix in \code{matrices}.
#' \item \code{plot_ecdf}: A ggplot2 empirical cumulative distribution function plot depicting the bootstrap distribution of Fst/FstMax for each matrix in \code{matrices}.
#' \item \code{test_kruskal_wallis}: Results of a Kruskal-Wallis test performed on the bootstrap distributions of Fst/FstMax. This test is a non-parametric statistical test of whether all provided bootstrap distributions are identically distributed.
#' \item \code{test_pairwise_wilcox}: Results of a Wilcoxon rank-sum test performed on the bootstrap distributions of Fst/FstMax. This test is a non-parameteric statistical test of whether \emph{each pairwise combination} of provided bootstrap distributions is identically distributed.
#' }
#' @examples
#' A = matrix(c(.4,.2,.4,
#'              .5,.4,.1,
#'              .6,.1,.3,
#'              .6,.3,.1),
#'            nrow = 4,
#'            byrow = TRUE)
#'
#' B = matrix(c(0,0,1,
#'              .5,.5,0,
#'              .5,.4,.1,
#'              .6,.2,.2,
#'              .6,.4,0),
#'            nrow = 5,
#'            byrow = TRUE)
#'
#' bs <- Q_bootstrap(matrices = list(A=A, B=B),
#'                   n_replicates = 10, K = 3)
#'
#'@importFrom dplyr %>%
#' @export
Q_bootstrap <- function(matrices, n_replicates, K){
  ## INPUT:
  # matrices = list of Q-matrices we seek to compare. Each of row is a vector of ancestry coefficients. Each matrix can be named (e.g. matrices = list(A=matrix(...), B=matrix(...)))
  # n_replicates = number of bootstrap replicates desired
  # K = # of ancestral clusters. Used to clean Q matrix.
  ## OUTPUT:
  # bootstrap_replicates = list of matrices of bootstrap replicates for each input matrix
  # statistics = dataframe; rows = bootstrap reps, columns = Matrix, Statistic, Value
  # histogram, CDF, boxplot of Fst/FstMax ratio bootstrap estimates

  # Do computations if matrices = a single matrix ---------------------------------------

  if(is.data.frame(matrices) | is.array(matrices)){

    n_matrix = 1

    names <- "Q"

    # Clean Q matrix - isolate ancestry coefficients
    if(ncol(matrices) > K){
      matrices <- matrices[,(ncol(matrices)-K+1):ncol(matrices)]
      if(sum(rowSums(matrices) != 1)){
        warning("The rows of the Q matrix (restricted to the last K columns) do not sum to 1.")
      }
    }

    bootstrap_matrices_Q <- list()
    matrix = matrices
    # Generate bootstrap data sets
    for(replicate in 1:n_replicates){
      bootstrap_matrices_Q[[replicate]] <- matrix[purrr::rdunif(n = nrow(matrix),
                                                         a = 1, b = nrow(matrix)),]
    }

    # Compute statistics for these reps
    stats_Q <- lapply(X = bootstrap_matrices_Q,
                      FUN = function(matrix) Q_stat(Q = matrix, K=ncol(matrix))) %>%
      unlist %>%
      matrix(ncol = 3, byrow = TRUE) %>% data.frame() %>%
      `colnames<-`(c("Fst", "FstMax", "ratio"))


    # Do computations if matrices = a list ---------------------------------------------

  }else if(is.list(matrices)){

    n_matrix = length(matrices)

    # List of names of matrices
    names <- if(sum(!is.na(names(matrices)))){names(matrices)}else{1:n_matrix}

    ## For each matrix: ##
    for(m in 1:n_matrix){

      bs_list <- list()
      matrix = matrices[[m]]

      # Check format of matrix
      if(ncol(matrix) > K){
        matrix <- matrix[,(ncol(matrix)-K+1):ncol(matrix)]
      }
      if(sum(round(rowSums(matrix),3) != 1)){
        warning(paste0(
          "At least one of the rows of Q matrix number ", m,
          " (restricted to the last K columns) does not sum to 1."))
      }

      # Generate bootstrap data sets
      for(replicate in 1:n_replicates){
        bs_list[[replicate]] <- matrix[purrr::rdunif(n = nrow(matrix),
                                              a = 1, b = nrow(matrix)),]
      }

      # Compute statistics for these reps
      stats <- lapply(X = bs_list, FUN = function(matrix) Q_stat(Q = matrix, K=ncol(matrix))) %>% unlist %>%
        matrix(ncol = 3, byrow = TRUE) %>% data.frame()
      colnames(stats) <- c("Fst", "FstMax", "ratio")

      # Sometimes, for values of Fst very close to 0 (i.e. order 10^-6, 10^-7), the
      # value of Fst ends up negative due to precision errors.
      # Find matrices for which this is the case, and replace them and their statistics

      while(sum(stats$ratio < 0)){
        # Which bootstrap matrices have negative values for Fst/FstMax?
        negatives <- which(stats$ratio < 0)

        # Replace those bootstrap replicates with new, random bootstrap replicates
        bs_list[negatives] <- lapply(X = 1:length(negatives),
                                     FUN = function(x){
                                       matrix[purrr::rdunif(n = nrow(matrix),
                                                     a = 1, b = nrow(matrix)),]
                                     })

        # Replace the corresponding entries of the statistics matrix
        stats[negatives,] <- lapply(X = bs_list[negatives],
                                    FUN = function(matrix){
                                      Q_stat(Q = matrix, K=ncol(matrix))
                                    }) %>%
          unlist %>%
          matrix(ncol = 3, byrow = TRUE) %>% data.frame()

      } # repeat this until there are no more errors

      # Name this dataset, based on the name of the matrices in the list or the entry number
      assign(paste0("stats_", names[m]), stats, pos = -1)
      assign(paste0("bootstrap_matrices_", names[m]), bs_list, pos = -1)
    }

  }else{
    stop("Error: The entry `matrices` must be a data frame, matrix, or array, or a list of these objects.")
  }

  # Make a dataset with all matrices' statistics:
  all_stats <- cbind(Matrix =
                       names %>%
                       lapply(function(name)rep(name, n_replicates)) %>%
                       unlist,
                     mget(paste0("stats_", names)) %>%
                       do.call(what = rbind, args = .) %>%
                       rbind)

  all_stats$Matrix <- factor(all_stats$Matrix, levels = unique(all_stats$Matrix))

  plot_ecdf <- ggplot2::ggplot(data = all_stats) +
    ggplot2::stat_ecdf(ggplot2::aes(x = ratio, color = Matrix)) +
    ggplot2::xlab(TeX('Fst/FstMax')) +
    ggplot2::ylab("Cumulative Probability") +
    ggplot2::xlim(0,1) + ggplot2::theme_bw() + ggplot2::scale_color_viridis_d()

  plot_boxplot <- ggplot2::ggplot(data = all_stats,
                                  ggplot2::aes(x = Matrix, y = ratio)) +
    ggplot2::geom_boxplot() +
    ggplot2::ylab(TeX('Fst/FstMax')) + ggplot2::xlab("") +
    ggplot2::theme_bw()

  plot_violin <- ggplot2::ggplot(data = all_stats,
                                 ggplot2::aes(x = Matrix, y = round(ratio,5))) +
    ggplot2::geom_violin(scale = "width") + ggplot2::geom_boxplot(width = 0.3) +
    ggplot2::ylab(TeX('Fst/FstMax')) + ggplot2::xlab("") +
    ggplot2::theme_bw()

  test_kruskal_wallis <- stats::kruskal.test(ratio ~ Matrix,
                                       data= all_stats)

  test_pairwise_wilcox <-  stats::pairwise.wilcox.test(x = all_stats$ratio,
                                                g = all_stats$Matrix,
                                                paired = FALSE)

  bootstrap_replicates <- mget(paste0("bootstrap_matrices_",
                                      names))

    return(list(bootstrap_replicates = bootstrap_replicates,
                statistics = all_stats,
                plot_boxplot = plot_boxplot,
                plot_violin = plot_violin,
                plot_ecdf = plot_ecdf,
                test_kruskal_wallis = test_kruskal_wallis,
                test_pairwise_wilcox = test_pairwise_wilcox))
}










# Q_simulate ----------------------------------------
#' Simulate one or more Q matrices using the Dirichlet distribution
#'
#' Simulates Q matrices by drawing vectors of membership coefficients from a Dirichlet distrubtion parameterized by two variables: \eqn{\alpha}, which controls variability, and \eqn{\lambda=(\lambda_1, \lambda_2, ...., \lambda_K)} which controls the mean of each of the K ancestry coefficients.
#'
#' @param alpha A number that sets the variability of the membership coefficients. The variance of coefficient k is Var[x_k] = \eqn{\lambda_k/(\alpha+1)}. Larger values of \eqn{\alpha} lead to lower variability.
#' @param lambda A vector that sets the mean membership of each ancestral cluster across the population. The vector should sum to 1.
#' @param rep The number of Q matrices to generate.
#' @param popsize The number of individuals to include in each Q matrix.
#' @param seed Optional; sets the random seed. Use if reproducibility of random results is desired.
#'
#' @return A data frame containing the simulated Q matrices. Each row represents a single simulated individual. The data frame has the following columns
#' \itemize{
#' \item \code{rep}: Which random Q matrix the row belongs to (a number between 1 and the parameter \code{rep})
#' \item \code{ind}: Which individual in each Q matrix the row corresponds to (a number between 1 and the parameter \code{popsize})
#' \item \code{alpha}: The alpha value used to simulate the Q matrix.
#' \item \code{Pop}: alpha_rep (where alpha and rep are the columns described above). Serves as a unique identifier for each Q matrix (useful if running simulations with many different values of \eqn{\alpha}).
#' \item \code{lambda1, lambda2, etc.}: Membership coefficients (sum to 1).
#' }
#'
#' @examples
#' # Simulate 100 random Q matrices.
#' # Each Q matrix has 100 individuals.
#' # On average these individuals have mean ancestry (1/2, 1/4, 1/4)
#' # from each of 3 ancestral clusters.
#' # The variance of each cluster x_i is Var[x_i] = lambda_i/(alpha + 1)
#' # Here lambda_1 = 1/2, lambda_2 = lambda_3 = 1/4
#'
#' Q_list <- Q_simulate(alpha = 1, lambda = c(1/2, 1/4, 1/4),
#'                     rep = 100, popsize = 50, seed = 1)
#'
#'@importFrom dplyr %>%
#' @export
Q_simulate <- function(alpha, lambda, rep, popsize, seed){

  # How many clusters are there?
  K = length(lambda)

  if(!missing(seed)){set.seed(seed)}

  # (a1, a2, ...) = (lambda1, lambda2, ...)*alpha parameterizes the Dirichlet distribution
  alpha_vec = alpha %>%
    data.frame("alpha" = .) %>%
    cbind(apply(X = ., MARGIN = 1, FUN = function(a) a*lambda) %>%
            matrix(byrow=TRUE, ncol=K)) %>% # this step necessary for K=1 case
    `colnames<-`(c("alpha", sapply(1:K, function(k) paste0("a",k))))

  ### Simulate Q matrices ###

  # Generate a data frame that uses the alpha_vec matrix above to simulate ancestry coefficients for a population of individuals for each alpha
  Q <- apply(X = alpha_vec %>%
               dplyr::select(-alpha), # Input is list of alpha vectors for each alpha
             MARGIN = 1, # call function on each row
             FUN = function(a){ # Draw ancestry vectors for a population for each alpha
               gtools::rdirichlet(n = popsize,
                          alpha = a) %>%
                 round(10)
             }) %>%
    data.frame %>% ### Restructure output so the columns are: individual, alphas, lambda1, lambda2, rep
    `colnames<-`(alpha_vec$alpha) %>%
    dplyr::mutate(ind = rep(1:popsize,K),
           lambda = (sapply(X = 1:K,
                       FUN = function(x) paste0("lambda", x)) %>%
                  lapply(function(lambda) rep(lambda, popsize)) %>%
                  unlist)) %>%
    tidyr::pivot_longer(cols = 1:length(alpha),
                 names_to = "alpha",
                 values_to = "Q") %>%
    tidyr::pivot_wider(names_from = lambda,
                values_from = Q) %>%
    cbind("rep" = 1, .)

  if(rep>1){
    # Repeat this process for iter times
    for(iter in 2:(rep)){
      Q <- Q %>%
        rbind(
          apply(X = alpha_vec %>%
                  dplyr::select(-alpha), # Input is list of alpha vectors for each alpha
                MARGIN = 1, # call function on each row
                FUN = function(a){ # Draw ancestry vectors for a population for each alpha
                  gtools::rdirichlet(n = popsize,
                             alpha = a) %>%
                    round(10)
                }) %>%
            data.frame %>% ### Restructure output so the columns are: individual, alphas, lambda1, lambda2, rep
            `colnames<-`(alpha_vec$alpha) %>%
            dplyr::mutate(ind = rep(1:popsize,K),
                   lambda = (sapply(X = 1:K,
                               FUN = function(x) paste0("lambda", x)) %>%
                          lapply(function(lambda) rep(lambda, popsize)) %>%
                          unlist)) %>%
            tidyr::pivot_longer(cols = 1:length(alpha),
                         names_to = "alpha",
                         values_to = "Q") %>%
            tidyr::pivot_wider(names_from = lambda,
                        values_from = Q) %>%
            cbind("rep" = iter, .)
        )
    }
  }

  # Give each population a unique identifier: alpha_rep
  Q <- Q %>% dplyr::mutate(ind = as.factor(ind),
                    alpha = as.numeric(alpha),
                    Pop = paste(round(alpha,3), rep, sep = "_") %>% as.factor,
                    rep = as.factor(rep),
                    .before = lambda1) %>%
    dplyr::arrange(rep, Pop, ind)
  return(Q)
}

