#!/usr/bin/env Rscript

# launch script to run all inversion parts in order
# author: Lewis Kunik

# start the clock to keep track of how long it takes for the script to run
ptm1_all <- proc.time()

# ~~~~~~~~~~~ Clear all existing data files before proceeding ~~~~~~~~~#

source("config.r")

# get all files in the outer directory
out_files <- paste0(out_path, list.files(out_path))

if(!dir.exists("H/"))
    dir.create("H/")

if(!is.na(lonlat_outer_file) & !dir.exists("H_outer/"))
    dir.create("H_outer/")

if(!dir.exists("HQ/"))
    dir.create("HQ/")

if(!dir.exists("out/"))
    dir.create("out/")


if (clear_H) {
    # H files
    if (length(list.files("H/")) > 0) {
        out_files <- c(out_files, paste0("H/", list.files("H/")))
    }

    # H files (outer domain)
    if (!is.na(lonlat_outer_file) & length(list.files("H_outer/")) > 0) {
        out_files <- c(out_files, paste0("H_outer/", list.files("H/")))
    }
}

# HQ files
if (length(list.files("HQ/")) > 0) {
    out_files <- c(out_files, paste0("HQ/", list.files("HQ/")))
}

# determine which of these files exist (redundant but good)
iFiles <- which(file.exists(out_files))

# remove existing files to clean the directory
invisible(sapply(out_files[iFiles], FUN = function(x) system(paste("rm", x))))


# ~~~~~~~~~~~~~~~~~ run scripts in order ~~~~~~~~~~~~~~~~~#

source("src/make_receptors.r")
source("src/aggregate_receptors.r")
source("src/make_sprior.r")
source("src/make_outer.r")
source("src/make_sigma.r")
source("src/make_sbio.r")
source("src/Hsplit.r")
source("src/make_Hs_bkgd.r")
source("src/make_sp_cov.r")
source("src/make_tmp_cov.r")
source("src/make_HQ.r")
source("src/make_bg.r")
source("src/make_z.r")
source("src/make_R.r")
source("src/make_zHsp.r")
source("src/inversion.r")
source("src/make_Vshat.r")
source("src/post_to_ncdf.r")
source("src/chi_sq.r")

bio_flux <- NA
outer_emiss <- NA
lonlat_outer <- NA
Hs_outer <- NA
Hs_bio <- NA
receptors_aggr <- NA

#load any necessary inputs before running code
obs_mat <- readRDS(obs_file)
bkgd_mat <- readRDS(bg_file)
lonlat_domain <- readRDS(lonlat_domain_file)
if(!is.na(lonlat_outer_file))
    lonlat_outer <- readRDS(lonlat_outer_file)


# 1. make receptor list (list of all footprint files used for observations)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("generating receptor list")
print("~~~~~~~~~~~~~~~~~~~~~~~~")

receptor_mat <- make_receptors(obs_mat, bkgd_mat)
n_obs <- nrow(receptor_mat)
receptor_times <- receptor_mat[,1] #time stamps are in seconds since 1970-01-01 00:00:00Z
receptor_files <- receptor_mat[,2]
saveRDS(receptor_mat, paste0(out_path, "receptors.rds"))


if(aggregate_receptors){
  receptors_aggr <- aggregate_receptors(receptor_files, receptor_times)
  n_obs <- nrow(receptors_aggr)
  saveRDS(receptors_aggr, paste0(out_path, "receptors_aggr.rds"))
}

# 2. make sprior (prior emissions vector)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("loading prior emissions")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
sprior <- make_sprior(prior_file, lonlat_domain)
saveRDS(sprior, paste0(out_path, "sprior.rds"))


# 3. make sigma (prior uncertainty vector)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("loading prior uncertainty")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
sigma <- make_sigma(sigma_file, lonlat_domain)
saveRDS(sigma, paste0(out_path, "sigma.rds"))


if (include_outer) {
    # 4. make outer (outer-domain emissions vector)
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    print("loading outer-domain emissions")
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    outer_emiss <- make_outer(outer_file, lonlat_outer)
    saveRDS(outer_emiss, paste0(out_path, "outer.rds"))

}

if (include_bio) {
    # 5. make sbio - biogenic flux vector
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    print("loading bio fluxes")
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    bio_flux <- make_sbio(bio_file, lonlat_outer)
    saveRDS(bio_flux, paste0(out_path, "sbio.rds"))

}

if (clear_H) {
    # 6. make H (footprint matrices)
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    print("running Hsplit")
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    Hsplit(receptor_files, receptor_times, lonlat_domain, lonlat_outer)
}

if (include_outer | include_bio) {
    # 7. derive biogenic and outer-domain additions to bkgd
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    print("generating background enhancements")
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    Hs_bkgd <- make_Hs_bkgd(lonlat_domain, lonlat_outer, bio_flux, outer_emiss, n_obs)

    if(include_outer)
        saveRDS(Hs_bkgd[["Hs_outer"]], paste0(out_path, "Hs_outer.rds"))

    if(include_bio)
        saveRDS(Hs_bkgd[["Hs_bio"]], paste0(out_path, "Hs_bio.rds"))
}

# 8. make spatial covariance matrix
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("generating spatial covariance matrix, E")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
E <- make_sp_cov(lonlat_domain)
saveRDS(E, paste0(out_path, "sp_cov.rds"))


# 9. make temporal covariance matrix
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("generating temporal covariance matrix, D")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
D <- make_tmp_cov(lonlat_domain)
saveRDS(D, paste0(out_path, "tmp_cov.rds"))

# 10. make HQ
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("generating HQ files")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
HQHt <- make_HQ(sigma, D, E, lonlat_domain, n_obs)
saveRDS(HQHt, paste0(out_path, "HQHt.rds"))

# 11. make bkgd
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("generating background")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
bkgd <- make_bg(receptor_files, receptor_times, bkgd_mat[,2], bkgd_mat[,1],
                Hs_outer, Hs_bio, receptors_aggr)
saveRDS(bkgd, paste0(out_path, "bg.rds"))


# 12. make z (anthropogenic enhancement values)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("generating enhancements")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
z <- make_z(receptor_files, receptor_times, obs_mat, bkgd, receptors_aggr)
saveRDS(z, paste0(out_path, "z.rds"))

# 13. make R (model data mismatch)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("creating model-data mismatch matrix, R")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
R <- make_R(receptor_files, receptor_times, lonlat_domain, receptors_aggr)
saveRDS(R, paste0(out_path, "R.rds"))

# 14. make z - Hsp
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("getting model-obs enhancement differences")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
zHsp_list <- make_zHsp(z, sprior)
zHsp <- zHsp_list[["zHsp"]]
Hsprior <- zHsp_list[["Hsprior"]]
saveRDS(zHsp, paste0(out_path, "zHsp.rds"))
saveRDS(Hsprior, paste0(out_path, "Hsprior.rds"))


# 15. make s_hat (optimized emissions vector)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("calculating additive corrections")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
s_hat <- make_post(sprior, R, zHsp, HQHt)
saveRDS(s_hat, paste0(out_path, "s_hat.rds"))

# 16. make_Vshat (posterior uncertainty - technically makes Vshat-bar, which is
# the grid-scale aggregated uncertainty)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("calculating posterior uncertainty")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
Vshat_list <- make_Vshat(lonlat_domain, sigma, D, E, R, HQHt)
Qsum <- Vshat_list[["Qsum"]]
Vshat <- Vshat_list[["Vshat_bar"]]
percent_uncert_reduction <- Vshat_list[["perc_unc_red"]]
saveRDS(Qsum, paste0(out_path, "Qsum.rds"))
saveRDS(Vshat, paste0(out_path, "Vshat_bar.rds"))
saveRDS(percent_uncert_reduction, paste0(out_path, "perc_unc_red.rds"))


# 17. convert posterior emissions/uncertainty into netcdf format for interpreting results
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("saving results to netcdf")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
post_to_ncdf(lonlat_domain, s_hat, sigma, D, E, R, HQHt)


if (compute_chi_sq) {
    # 18. calculate Chi-squared
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    print("calculating reduced Chi squared")
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    reduced_chi_sq <- calc_chi_sq(sprior, s_hat, sigma, D, E, R, z)
    saveRDS(reduced_chi_sq, paste0(out_path, "chi_sq.rds"))
}


# ~~~~~~~~~~~~~~~~ get elapsed time data ~~~~~~~~~~~~~~~~#
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
ptm2_all <- proc.time()
elapsed_seconds <- as.numeric(ptm2_all["elapsed"] - ptm1_all["elapsed"])
e_mins <- round(elapsed_seconds/60)
e_secs <- round(elapsed_seconds%%60, digits = 1)
print(paste0("elapsed time: ", e_mins, " minutes, ", e_secs, " seconds"))
