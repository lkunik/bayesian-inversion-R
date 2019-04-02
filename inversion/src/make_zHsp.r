# calculate observed minus prior enhancements
# author: Lewis Kunik

## prerequisite scripts:
##  make_sprior.r
##  Hsplit.r
##  make_z.r
##
## output files:
##  zHsp.rds - file containing a vector of values equal to (z - Hsp)

# run dependent scripts
source("config.r")
source("src/read_sparse_h.r")


make_zHsp <- function(z,
                      prior_emiss) {

    # get num obs from length of z
    nobs <- length(z)

    # format into matrix so different time-steps are easily accessible
    sprior_mat <- matrix(prior_emiss, nrow = ntimes, byrow = T)
    ncells <- ncol(sprior_mat)

    # get Hsprior by convolving footprints
    Hsprior <- rep(0, nobs)

    path <- "H/"
    # convolve each H file and add to running total of Hsbio
    for (ii in 1:ntimes) {

        # grab H from file
        H_file <- paste0(path, "H", formatC(ii, width = 3, flag = "0"), ".rds")
        Hi <- read_sparse_h(H_file, nobs, ncells)

        # convolve footprints with prior emissions
        Hsprior <- Hsprior + Hi %*% sprior_mat[ii, ]
    }

    # subtract Hsprior from z
    zHsp <- z - Hsprior

    ret_list <- list()

    #return both vars in list
    ret_list[["Hsprior"]] <- Hsprior
    ret_list[["zHsp"]] <- zHsp

    ret_list

}
