# create H matrices and save to files in dense-data format.
# Each row of H contains footprint values for its corresponding observation in
# the obs vector. Columns correspond to footprint values of corresponding gridcells
# at all timesteps.
# Domain-sized files are stored in H/ directory outer-domain footprint files
# are stored in H_outer/ directory footprints and are used to calculate bkgd contributions

# author: Lewis Kunik

## prerequisite scripts:
## make_receptors.r
##
## output files:
## H/H*.rds - array of files, 1 per inversion timestep, containing condensed set
##  of footprint values [mixing ratio/flux unit]
##  (H data files contain 3 columns: obs_index (row in H matrix), cell_index (column
##  in H matrix), foot_value) (arranged this way to omit cells with foot value = 0,
##  saves a lot of space because H matrix is big)
##
## H_outer/H*.rds - array of files, 1 per inversion timestep, same format as
##  H/H*.rds but for footprints including the outer-domain. Larger files.
##  footprint values still expressed as [mixing ratio/flux unit]

# load package dependencies
library(ncdf4)
library(data.table)
library(lubridate)
library(raster)

# run dependent scripts
source("config.r")

# Start the clock!
ptm1 <- proc.time()

#if a lon/lat file is supplied, make H_outer files
make_H_outer <- !is.na(lonlat_outer_file)

# ~~~~~~~~~~~~~~~~~~~~~~~ load required files ~~~~~~~~~~~~~~~~~~~~~~~#

# load in receptor info
recep_file <- paste0(out_path, "receptors.rds")
recep_mat <- readRDS(recep_file)
recep_times <- as.numeric(recep_mat[, 1])
class(recep_times) <- c("POSIXt", "POSIXct")
attributes(recep_times)$tzone <- "UTC"
receptors <- recep_mat[, 2]
nobs <- length(receptors)

# site info
isites <- lapply(sites, FUN = function(x) grep(x, receptors))

# load in lonlat_domain file specified in config.r to make mask for domain cells
lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)

# load in sample foot file to get the footprint domain details
nc_f <- nc_open(receptors[1])
nc_lat <- round(ncvar_get(nc_f, "lat"), round_digs)  #naming based on STILT-R conventions
nc_lon <- round(ncvar_get(nc_f, "lon"), round_digs)  #naming based on STILT-R conventions
nlat <- length(nc_lat)
nlon <- length(nc_lon)
nc_close(nc_f)
lonlat_foot <- expand.grid(nc_lon, nc_lat)


# create index vector to act as a mask for the outer H domain (cells which
# overlap with bio fluxes)
iDomain <- apply(lonlat_domain, FUN = function(x) which((lonlat_foot[, 1] == x[1]) &
                                                          (lonlat_foot[, 2] == x[2])), MARGIN = 1)

#if outer domain is toggled (for far-field fluxes) define here
if (make_H_outer) {
  # load in lonlat_outer file to make mask for outer cells
  lonlat_outer <- readRDS(lonlat_outer_file)
  ncells_outer <- nrow(lonlat_outer)

  #apply rectangular index filter to the footprint domain
  ifoot_complete <- array(1:nrow(lonlat_foot), dim = c(nlon, nlat))
  in_lon <- unique(lonlat_outer[, 1])
  in_lat <- unique(lonlat_outer[, 2])
  imin_lat <- which(nc_lat == min(in_lat))
  imax_lat <- which(nc_lat == max(in_lat))
  imin_lon <- which(nc_lon == min(in_lon))
  imax_lon <- which(nc_lon == max(in_lon))

  ifoot_outer <- ifoot_complete[imin_lon:imax_lon, imin_lat:imax_lat]
  iOuter <- sort(as.vector(ifoot_outer))  #specify the indices for the model outer domain
} #end if make_H_outer


# ~~~~~~~~~~~~~~~~~~~~~~~ create time bins ~~~~~~~~~~~~~~~~~~~~~~~#

# specify start/end from config.r variables
y1 <- flux_year_start
y2 <- flux_year_end

m1 <- flux_month_start
m2 <- flux_month_end

d1 <- flux_day_start
d2 <- flux_day_end

h1 <- flux_hour_start
h2 <- flux_hour_end

mn1 <- flux_min_start
mn2 <- flux_min_end

# list all desired obs times
time_bins <- seq(from = ISOdatetime(y1, m1, d1, h1, mn1, 0, tz = "UTC"),
                 to = ISOdatetime(y2,m2, d2, h2, mn2, 0, tz = "UTC"), by = flux_t_res) # flux_t_res defined in config.r


# create H files for each timestep - we will use TXT files and later append to
# them with each footprint file. File I/O is faster this way than constantly
# saving/rewriting RDS files. At the end we will convert TXT to RDS for more
# compact storage/quicker loading in R for later

# note: ntimes defined in config.r
for (ii in 1:ntimes) {

  Hii <- array(0, dim = c(0, 3))
  colnames(Hii) <- c("obs_index", "cell_index", "foot_value")

  # make empty rds files for each timestep
  tmpfile <- paste0("H/H", formatC(ii, width = 3, flag = "0"), ".txt")
  fileWrite <- file(tmpfile, "wt")
  write.table(Hii, fileWrite, row.names = F)
  close(fileWrite)

  if (make_H_outer) {

    # do the same for outer files
    tmpfile_outer <- paste0("H_outer/H", formatC(ii, width = 3, flag = "0"), ".txt")
    fileWrite <- file(tmpfile_outer, "wt")
    write.table(Hii, fileWrite, row.names = F)
    close(fileWrite)

  } #end if make_H_outer

}


# ~~~~~~~~~~~~~~~~~~~~~~~ assemble H matrix ~~~~~~~~~~~~~~~~~~~~~~~#
for (ii in 1:nobs) {

    print(paste("loading footprint file:", receptors[ii]))

    # load in footprint file and obtain the vars we need
    nc_f <- nc_open(receptors[ii])
    nc_lat <- round(ncvar_get(nc_f, "lat"), round_digs)
    nc_lon <- round(ncvar_get(nc_f, "lon"), round_digs)
    nc_time <- ncvar_get(nc_f, "time")
    class(nc_time) <- c("POSIXt", "POSIXct")
    attributes(nc_time)$tzone <- "UTC"
    nc_foot <- ncvar_get(nc_f, "foot")  #naming based on STILT-R convention
    nc_close(nc_f)  #close the netcdf file

    nlat <- length(nc_lat)
    nlon <- length(nc_lon)

    # note: nc_foot has dims lat x lon x time
    if(foot_dim_lonxlat){
      default_xy_dims <- c(nlon, nlat)
    } else {
      default_xy_dims <- c(nlat, nlon)
    }

    # break up footprint times into corresponding time bins to appropriately assign columns
    times_cut <- as.POSIXct(cut(nc_time, breaks = time_bins), tz = "UTC")
    if (all(is.na(times_cut)))
        next

    # number of time bins this footprint's time dimension falls into
    nbins <- length(unique(times_cut))

    # nc_foot should have >1 timestep for all files, but this is a safety check
    if (length(nc_time) > 1) {

        # for each bin, sum the foot grids of all timesteps that fall into that bin
        nc_foot_sum <- array(0, dim = c(default_xy_dims, nbins))
        for (jj in 1:nbins) {
            ibin <- which(times_cut == unique(times_cut)[jj])  #get indices of timesteps corresponding to this bin
            if(length(ibin) > 1){
              nc_foot_sum[, , jj] <- apply(nc_foot[ , , ibin], FUN = sum, MARGIN = c(1,2))  #sum over bin's timesteps
            } else{
              nc_foot_sum[, , jj] <- nc_foot[ , , ibin] #get the lat/lon components even if 1 timestep
            }
        }

    } else {
        # if footprint has 1 time step, then dim(nc_foot) = lon x lat, and must be
        # converted to lon x lat x 1
        nc_foot_sum <- array(0, dim = c(default_xy_dims, 1))
        nc_foot_sum[, , 1] <- nc_foot
    }

  # create 2-D arrays with dims = c(ncells, ntimebins).
  if (foot_dim_lonxlat) {
    nc_foot_vec <- apply(nc_foot_sum, FUN = function(x) as.vector(x)[iDomain],
                         MARGIN = 3)  #iDomain denotes inner-domain cells

    if (make_H_outer) {

      nc_foot_vec_outer <- apply(nc_foot_sum, FUN = function(x) as.vector(x)[iOuter],
                                 MARGIN = 3)  #iOuter denotes outer-domain cells
    } #end if make_H_outer
  } else {
    nc_foot_vec <- apply(nc_foot_sum, FUN = function(x) as.vector(t(x))[iDomain],
                         MARGIN = 3)  #iDomain denotes inner-domain cells

    if (make_H_outer) {

      nc_foot_vec_outer <- apply(nc_foot_sum, FUN = function(x) as.vector(t(x))[iOuter],
                                 MARGIN = 3)  #iOuter denotes outer-domain cells
    } #end if make_H_outer
  }

  # append to Hfiles in the H directory
  for (jj in 1:nbins) {

    # itime denotes the index of inversion timestep i.e. which H file to append to
    itime <- itime <- which(time_bins == unique(times_cut)[jj])

    # INNER DOMAIN: format the data to include only nonzero foot values (dense format)
    inonzero <- which(nc_foot_vec[, jj] != 0)  #gets the cell indices of non-zero foot vals
    Hvec_nonzero <- nc_foot_vec[inonzero, jj]
    Hsave <- cbind(rep(ii, length(inonzero)), inonzero, Hvec_nonzero)  #ii is obs index

    # append this obs' values to the file
    fwrite(as.data.frame(Hsave), paste0("H/H", formatC(itime, width = 3, flag = "0"), ".txt"),
           row.names = F, col.names = F, append = T, sep = ' ')

    if (make_H_outer) {

      # OUTER DOMAIN: format the data to include only nonzero foot values (dense format)
      inonzero_outer <- which(nc_foot_vec_outer[, jj] != 0)  #gets the cell indices of non-zero foot vals
      Hvec_nonzero_outer <- nc_foot_vec_outer[inonzero_outer, jj]
      Hsave_outer <- cbind(rep(ii, length(inonzero_outer)), inonzero_outer, Hvec_nonzero_outer)  #ii is obs index

      # append this obs' values to the file
      fwrite(as.data.frame(Hsave_outer), paste0("H_outer/H", formatC(itime, width = 3, flag = "0"),
                                 ".txt"), row.names = F, col.names = F, append = T, sep = ' ')
    } #end if make_H_outer

  }
}


print("replacing footprint .txt files with .rds files")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if aggregating obs: now combine rows (receptors) into single daily averages
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load in time vars to specify the span of time over which we want subsetted obs
y1 <- obs_year_start
y2 <- obs_year_end

m1 <- obs_month_start
m2 <- obs_month_end

d1 <- obs_day_start
d2 <- obs_day_end

h1 <- obs_hour_start
h2 <- obs_hour_end

mn1 <- obs_min_start
mn2 <- obs_min_end

# prepare to aggregate footprint files based on day
time_bins_daily <- seq(from = ISOdatetime(y1, m1, d1, h1, mn1, 0, tz = "UTC"), to = ISOdatetime(y2,
                       m2, d2, h2, mn2, 0, tz = "UTC") + 3600, by = 24 * 3600)


times_cut_day <- as.POSIXct(cut(recep_times, breaks = time_bins_daily), tz = "UTC")

#read in receptor aggregation file:
receps_aggr <- readRDS(paste0(out_path, "receptors_aggr.rds"))
nobs_aggr <- nrow(receps_aggr)

# run through each timestep, load H txt file, combine times into daily
# footprints, and re-save as RDS

# note: ntimes defined in config.r
for (ii in 1:ntimes) {

  print(paste("re-constructing compact H matrix for timestep", ii))

  # inner H files
  Hi <- fread(paste0("H/H", formatC(ii, width = 3, flag = "0"), ".txt"))
  Hsave <- Hi

  if (make_H_outer) {
    # outer H files
    Hi_outer <- fread(paste0("H_outer/H", formatC(ii, width = 3, flag = "0"), ".txt"))
    Hsave_outer <- Hi_outer
  } #end if make_H_outer

  # if aggregation is toggled, fill in footprint values and aggregate daily by
  # subsetted times
  if (aggregate_obs) {

    #this will hold the old H matrix
    Hi_full <- array(0, dim = c(nobs, ncells))

    #this will hold the new H matrix
    Hi_new <- array(0, dim = c(nobs_aggr, ncells))

    iobs <- as.matrix(Hi[, 1]) #row indices
    icell <- as.matrix(Hi[, 2]) #column indices
    Hi_full[cbind(iobs, icell)] <- as.matrix(Hi[, 3]) #fill with the non-zero values

    #this will hold the new dense-format daily-aggregated footprint
    Hsave <- array(0, dim = c(0, 3))


    if (make_H_outer) {

      #this will hold the old H matrix
      Hi_full_outer <- array(0, dim = c(nobs, ncells_outer))

      #this will hold the new H matrix
      Hi_new_outer <- array(0, dim = c(nobs_aggr, ncells_outer))
      #Hi_new_outer$Hrows <- rep(0, nobs_aggr)

      iobs_outer <- as.matrix(Hi_outer[, 1]) #row indices
      icell_outer <- as.matrix(Hi_outer[, 2]) #column indices
      Hi_full_outer[cbind(iobs_outer, icell_outer)] <- as.matrix(Hi_outer[, 3]) #fill with non-zero vals

      #this will hold the new dense-format daily-aggregated footprint
      Hsave_outer <- array(0, dim = c(0, 3))
    } #end if make_H_outer

    row_count <- 1

    # note: ndays and nsites defined in config.r
    for (jj in 1:ndays) {
      for (kk in 1:nsites) {

        isite <- isites[[kk]]  #find which receps correspond to this site
        iday <- which(times_cut_day == time_bins_daily[jj])  #find which recep times correspond to this day
        idaysite <- isite[which(isite %in% iday)]  #get the overlap of these two - THIS site on THIS day

        # if this there are no obs for this site on this day, skip
        if (length(idaysite) < min_agg_obs)
          next

        # now take the average footprint for this day and set this as the row for this
        # day/site in the new H matrix
        if (length(idaysite) > 1) {
          Hi_new[row_count,] <- colMeans(Hi_full[idaysite,])
        } else {
          Hi_new[row_count,] <- Hi_full[idaysite, ]
        }

        if (make_H_outer) {

          if (length(idaysite) > 1) {
            Hi_new_outer[row_count,] <- colMeans(Hi_full_outer[idaysite, ])
          } else {
            Hi_new_outer[row_count,] <- Hi_full_outer[idaysite, ]
          }

        } #end if make_H_outer

        row_count <- row_count + 1

      }
    }

    for (jj in 1:nrow(Hi_new)) {

      inonzero <- which(Hi_new[jj, ] != 0)  #gets the cell indices of non-zero foot vals
      Hvec_nonzero <- Hi_new[jj, inonzero]
      Hsave <- rbind(Hsave, cbind(rep(jj, length(inonzero)), inonzero, Hvec_nonzero))  #jj is obs index

      if (make_H_outer) {

        inonzero_outer <- which(Hi_new_outer[jj, ] != 0)  #gets the cell indices of non-zero foot vals
        Hvec_nonzero_outer <- Hi_new_outer[jj, inonzero_outer]
        Hsave_outer <- rbind(Hsave_outer, cbind(rep(jj, length(inonzero_outer)),
                                                inonzero_outer, Hvec_nonzero_outer))  #jj is obs index

      } #end if make_H_outer

    }
  }

  # save files whether or not aggregation has occurred
  saveRDS(Hsave, paste0("H/H", formatC(ii, width = 3, flag = "0"), ".rds"))
  system(paste0("rm H/H", formatC(ii, width = 3, flag = "0"), ".txt"))  #remove original text file which is large

  if (make_H_outer) {

    saveRDS(Hsave_outer, paste0("H_outer/H", formatC(ii, width = 3, flag = "0"), ".rds"))
    system(paste0("rm H_outer/H", formatC(ii, width = 3, flag = "0"), ".txt"))  #remove original text file which is large

  } #end if make_H_outer
}


# get elapsed time data
ptm2 <- proc.time()
elapsed_seconds <- as.numeric(ptm2["elapsed"] - ptm1["elapsed"])
e_mins <- round(elapsed_seconds/60)
e_secs <- round(elapsed_seconds%%60, digits = 1)
print(paste0("Hsplit elapsed time: ", e_mins, " minutes, ", e_secs, " seconds"))
