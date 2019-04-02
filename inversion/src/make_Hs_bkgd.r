# convolve outer domain footprints with far-field anthro and bio emissions
# author: Lewis Kunik

## prerequisite scripts:
##  make_receptors.r
##  make_outer.r
##  make_sbio.r
##  Hsplit.r
##
## output files:
##  Hs_outer.rds - far-field emission contributions from outer_emiss.nc file,
##  corresponding to observations from the receptor file
##  Hs_bio.rds - biogenic flux contributions from bio_flux.nc file,
##  corresponding to observations from the receptor file

# run dependent scripts
source("config.r")
source("src/read_sparse_h.r")


make_Hs_bkgd <- function(lonlat_domain,
                         lonlat_outer,
                         bio_flux = NA,
                         outer_emiss = NA,
                         n_obs){




    ncells_outer <- nrow(lonlat_outer)

    iDomain <- apply(lonlat_domain, FUN = function(x) which((lonlat_outer[, 1] %in% x[1]) &
        (lonlat_outer[, 2] %in% x[2])), MARGIN = 1)

    # ~~~~~~~~~~~~~~~ Load emissions files ~~~~~~~~~~~~~~~#

    if (include_outer) {

        outer_mat <- matrix(outer_emiss, nrow = ntimes, byrow = T)

        # mask out the inner domain cells, so that the outer emiss contribution only
        # considers outer-domain cells
        outer_mat[, iDomain] <- 0

        #set up the convolved vector to add to
        Hs_outer <- rep(0, nobs)
    }

    if (include_bio) {

        sbio_mat <- matrix(bio_flux, nrow = ntimes, byrow = T)

        # set up the convolved vector to add to
        Hsbio <- rep(0, nobs)
    }

    # convolve each H file and add to running total of Hsbio
    for (ii in 1:ntimes) {

        # get H far-field file for time = ii
        H_file <- readRDS(paste0("H_outer/H", formatC(ii, width = 3, flag = "0"), ".rds"))
        Hi <- read_sparse_h(H_file, nobs, ncells_outer)

        # add outer outer domain contribution, if toggled
        if (include_outer)
            Hs_outer <- Hs_outer + Hi %*% outer_mat[ii, ]

        # add biospheric contributions to observations, if toggled
        if (include_bio)
            Hsbio <- Hsbio + Hi %*% sbio_mat[ii, ]

    }  #end ntimes for-loop


    ret_list <- list()

    # save the files as vectors of enhancements for each observation in obs space
    if (include_outer) {

        ret_list[["Hs_outer"]] <- Hs_outer

        #filepath <- paste0(out_path, "Hs_outer.rds")
        #saveRDS(Hs_outer, filepath)
    }

    if (include_bio) {

        ret_list[["Hs_outer"]] <- Hs_outer
        #filepath <- paste0(out_path, "Hs_bio.rds")
        #saveRDS(Hsbio, filepath)

    }

    ret_list

}
