## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warnings = FALSE,
  message = FALSE,
  fig.width = 7,
  comment = "#>"
)

library(dplyr)

Q_example <- matrix(
  c("1", "1", "(x)", "1", ":", "1", "0", "0", 
    "2", "2", "(x)", "1", ":", "1",  "0", "0",
    "3", "3", "(x)", "2", ":", "0.0137", "0.7567", "0.2296",
    "4", "4", "(x)", "2", ":", "0", "0.7325", "0.2674",
    "5", "5", "(x)", "2", ":", "0.0224", "0.7279", "0.2496",
    "6", "6", "(x)", "3", ":", "0.0222", "0.7372", "0.2406",
    "7", "7", "(x)", "3", ":", "0", "0.7396", "0.2603",
    "8", "8", "(x)", "3", ":", "0.0231", "0.7496", "0.2272",
    "9", "9", "(x)", "3", ":", "0.046", "0.6861", "0.2678",
    "10", "10", "(x)", "3", ":", "0.0843", "0.6739", "0.2418"),
  ncol = 8, 
  byrow = TRUE
) %>% 
  data.frame %>% 
  mutate(X6 = as.numeric(X6), 
         X7 = as.numeric(X7))


## ----echo = FALSE-------------------------------------------------------------
knitr::kable(Q_example, col.names = NULL)

## ----setup--------------------------------------------------------------------
library(FSTruct)
library(dplyr)

## ---- eval = FALSE------------------------------------------------------------
#  Q_matrix_name <- data.table::fread("file path to your Q matrix.output")

## -----------------------------------------------------------------------------
A = Q_simulate(alpha = .1, lambda = c(.75, .25), rep = 1, popsize = 20, seed = 1)

B = Q_simulate(alpha = .1, lambda = c(.75, .25), rep = 1, popsize = 20, seed = 2)

C = Q_simulate(alpha = 1, lambda = c(.75, .25), rep = 1, popsize = 20, seed = 3)

D = Q_simulate(alpha = 1, lambda = c(.75, .25), rep = 1, popsize = 20, seed = 4)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(A)

## -----------------------------------------------------------------------------
# Generate and modify a plot for each Q matrix

# arrange(var) sorts the individuals in order of that variable
# scale_fill_brewer() modifies the color scheme
# Because Q_plot() outputs a ggplot, you can do lots of modifications!

plot_A <- Q_plot(Q = A %>% arrange(lambda1),
                 K=2) + 
  ggplot2::scale_fill_brewer("Blues")

plot_B <- Q_plot(Q = B %>% arrange(lambda1), 
                 K=2) + 
  ggplot2::scale_fill_brewer("Blues")

plot_C <- Q_plot(Q = C %>% arrange(lambda1), 
                 K=2) + 
  ggplot2::scale_fill_brewer("Blues")

plot_D <- Q_plot(Q = D %>% arrange(lambda1), 
                 K=2) + 
  ggplot2::scale_fill_brewer("Blues")

# Display them in a grid
cowplot::plot_grid(plot_A, plot_B, plot_C, plot_D,
                   labels = "AUTO", 
                   nrow = 1,
                   vjust = 2) 

## -----------------------------------------------------------------------------
Q_stat(Q = A, K = 2)
Q_stat(Q = B, K = 2)
Q_stat(Q = C, K = 2)
Q_stat(Q = D, K = 2)

## -----------------------------------------------------------------------------
bs <- Q_bootstrap(matrices = list(A = A, # name = matrix
                                  B = B,
                                  C = C,
                                  D = D), 
                  n_replicates = 100,
                  K = 2, 
                  seed = 1)

## ---- fig.height=4------------------------------------------------------------
cowplot::plot_grid(bs$plot_boxplot + ggplot2::ggtitle("Box Plot"), 
                   bs$plot_violin + ggplot2::ggtitle("Violin Plot"), 
                   bs$plot_ecdf + ggplot2::ggtitle("ECDF Plot"), 
                   nrow = 2, rel_widths = c(1,1,2))

## -----------------------------------------------------------------------------
bs$test_kruskal_wallis

bs$test_pairwise_wilcox

