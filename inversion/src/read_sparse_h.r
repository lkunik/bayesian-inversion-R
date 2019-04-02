

# ~~~~~~~~~~~~~ define a function to read sparse H file ~~~~~~~~~~~~~~~ #
read_sparse_h <- function(H_file, nobs, ncells) {

    Hi_tmp <- readRDS(H_file)
    # Populate the H-slice matrix (nobs x ncells) with zeros
    Hi <- array(0, dim = c(nobs, ncells))

    # checking if Hi is empty
    if (length(Hi_tmp[, 1]) > 0) {
        # for every value in H file, locate where that exists and put a value there
        iobs <- Hi_tmp[, 1]
        icell <- Hi_tmp[, 2]
        Hi[cbind(iobs, icell)] <- Hi_tmp[, 3]
    }

    # return H for this timestep
    return(Hi)
}
