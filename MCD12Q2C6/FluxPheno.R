#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Josh Gray May 2020
# Smoother/adaptor so we can run MCD12Q2C6 algo on FLUXNET2015 GPP time series
# creates .Rdata files containing all Q2-style flux metrics from GPP
# Some FLUXNET2015 variables:
# GPP_DT_VUT_REF: GPP from datyime partioning method, ref selected from GPP versions w/ MEF
# GPP_DT_VUT_MEAN: GPP from datyime partioning method, avg from GPP versions "each from corresponding NEE_VUT_XX version"
# GPP_DT_VUT_SE: s.e. for GPP "calculated as (SD(GPP_DT_VUT_XX) / SQRT(40))"
# there are _CUT_ versions of all of the above, with the same interpretation, but constant ustar
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# packages
library(data.table)
library(viridis)
library(scales)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# plots a background spanning a calendar year
PlotYearBackground <- function(year, col=rgb(0.92, 0.92, 0.92)){
    poly_x = c(rep(as.Date(paste(year, "-1-1", sep="")), 2), rep(as.Date(paste(year, "-12-31", sep="")), 2))
    poly_y <- c(par()$usr[c(3:4, 4:3)])
    polygon(poly_x, poly_y, col=col, border=NA)    
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# helper function: how many days in year
DaysInYear <- function(year) as.integer(diff(as.Date(paste(year, c("-1-1", "-12-31"), sep="")))) + 1

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Spline smooths with optional weighting (weight_var, weight_func, weight_func_args) over
# the year_of_interest and pad_days before/after. If fake_head_tail is TRUE, and year_of_interest
# is the first or last year, then pad_days data is permuted from year_of_interest
# Returns a list with elements: date, x, and spl_weights. The object x mimics what
# Get3YearDataChunk returns: c(gpp, resids, qc) QC=0 observed data, Qc=1 permuted,
# QC=missing. This allows for using MCD12Q2C6 AnnualPhenology() function
GetSmoothFluxYear <- function(flux_dataset, year_of_interest, gpp_metric, weight_var=NA, time_var="TIMESTAMP", weight_func = function(x,c) 1/(x+c$c), weight_func_args = list(c=0.1), const_weight=1, fake_head_tail=FALSE, pad_days=190, max_span=45, spline_spar=0.5, make_plot=FALSE){
    #---------------------
    # Read file and check for gpp_metric, weight_var, and time_var. Returns NA for failure
    # Subset to date, gpp, and weight columns, apply weight func to weight_var
    fluxDT <- fread(flux_dataset, na.strings=c("-9999", "9999"))
    sdcols <- as.character(na.omit(c(time_var, gpp_metric, weight_var)))
    if(!all(sdcols %in% names(fluxDT))) return(NA) # missing one or more variables
    fluxDT <- fluxDT[, .SD, .SDcols=sdcols]
    if(ncol(fluxDT) != 3) fluxDT[, w := const_weight] # no weight_var: set const weights
    setnames(fluxDT, c("date", "gpp", "w"))
    fluxDT[, spl_weight := weight_func(w, weight_func_args)]
    fluxDT[, qc := 0] # all data retrieved get QA=0
    
    #---------------------
    # Create YOI DT; return NA if YOI not in data. Fill w/ NA to get daily YOI DT
    fluxDT[, date := as.Date(as.character(date), format="%Y%m%d")]
    yoi_min_max <- as.Date(paste(year_of_interest, c("-1-1", "-12-31"), sep=""))
    fluxDT_yoi <- fluxDT[between(date, yoi_min_max[1], yoi_min_max[2])]
    if(nrow(fluxDT_yoi) == 0) return(NA) # yoi not in data range
    tmpDT_yoi <- data.table(date=seq.Date(yoi_min_max[1], yoi_min_max[2], by="day"))
    fluxDT_yoi <- merge(fluxDT_yoi, tmpDT_yoi, all=TRUE) # add NA records for missing days

    #---------------------
    # Get pad_days of data from year before and after and append to fluxDT_yoi
    # Will be faked from yoi if(fake_head_tail) AND yoi is first/last year in data
    head_min_max <- as.Date(paste(c(year_of_interest - 1, year_of_interest), c(DaysInYear(year_of_interest - 1) - pad_days, 1), sep=""), format="%Y%j")
    tail_min_max <- as.Date(paste(c(year_of_interest, year_of_interest + 1), c(DaysInYear(year_of_interest), pad_days + 1), sep=""), format="%Y%j")
    
    # do the head
    if(fake_head_tail & (year_of_interest == fluxDT[, as.integer(strftime(min(date, na.rm=T), format="%Y"))])){
        # if first year of data, and we're faking, grab pad_days from end of year_of_interest
        fluxDT_head <- fluxDT_yoi[between(date, yoi_min_max[2] - pad_days, tail_min_max[1] + 1, incbounds = FALSE)]
        fluxDT_head[, date := date - DaysInYear(year_of_interest)] # fake the date
        fluxDT_head[, qc := 1] # indicate this is fill data w/ QC=1
    }else{
        # if not first year of data, or we're not faking, get available real data (may be empty)
        fluxDT_head <- fluxDT[between(date, head_min_max[1], head_min_max[2], incbounds = FALSE)]
    }
    # create daily series, filling w/ NA as necessary
    tmpDT_head <- data.table(date=seq.Date(head_min_max[1] + 1, head_min_max[2] - 1, by="day"))
    fluxDT_head <- merge(fluxDT_head, tmpDT_head, all=TRUE)

    # do the tail
    if(fake_head_tail & (year_of_interest == fluxDT[, as.integer(strftime(max(date, na.rm=T), format="%Y"))])){
        # if last year of data, and we're faking, grab pad_days from start of year_of_interest
        fluxDT_tail <- fluxDT_yoi[between(date, head_min_max[2] - 1, yoi_min_max[1] + pad_days, incbounds = FALSE)]
        fluxDT_tail[, date := date + DaysInYear(year_of_interest)] # fake the date
        fluxDT_tail[, qc := 1] # indicate this is fill data w/ QC=1
    }else{
        # if not last year of data, or we're not faking, get available real data (may be empty)
        fluxDT_tail <- fluxDT[between(date, tail_min_max[1], tail_min_max[2], incbounds = FALSE)]
    }
    # create daily series, filling w/ NA as necessary
    tmpDT_tail <- data.table(date=seq.Date(tail_min_max[1] + 1, tail_min_max[2] - 1, by="day"))
    fluxDT_tail <- merge(fluxDT_tail, tmpDT_tail, all=TRUE)

    # append head/tail to yoi and redefine fluxDT as final data
    fluxDT <- rbind(fluxDT_yoi, fluxDT_head, fluxDT_tail)
    fluxDT[is.na(spl_weight), spl_weight := 0] # added days will have NA weights, which we set to 0
    fluxDT[is.na(qc) | is.na(gpp), qc := 2] # QC=2 for missing days
    setkey(fluxDT, date)
    
    #---------------------
    # check for exceedance of maximum missing data span
    na_rle <- fluxDT[, rle(is.na(gpp))] # convert to run length encoding
    if(any(na_rle$lengths[na_rle$values] > max_span)) return(NA)
    
    #---------------------
    # fit a weighted smoothing spline, and then predict over the range
    flux_spl <- smooth.spline(fluxDT[!is.na(gpp), gpp] ~ fluxDT[!is.na(gpp), as.integer(date)], w=fluxDT[!is.na(gpp), spl_weight], spar=spline_spar)
    fluxDT[, gpp_smooth := predict(flux_spl, as.integer(date))$y]
    fluxDT[, gpp_resids := gpp_smooth - gpp]

    #---------------------
    # make a plot of the splined and raw data
    if(make_plot){
        pt_col <- "grey60"
        line_col <- "grey30"
        cexs <- rescale(fluxDT[, spl_weight], c(0.2, 1.5))
        pchs <- c(16, 1, NA)[fluxDT[, qc + 1]]
        plot(fluxDT[, .(date, gpp)], type="n", xlab="", ylab=gpp_metric)
        PlotYearBackground(year_of_interest)
        points(fluxDT[, .(date, gpp)], pch=pchs, cex=cexs, col=pt_col, xlab="", ylab=gpp_metric)
        points(fluxDT[, .(date, gpp_smooth)], type="l", col=line_col, lwd=2)
        legend("topleft", legend=c("raw", "permuted", "spline"), pch=c(16, 1, NA), lwd=c(NA, NA, 2), col=c(rep(pt_col, 2), line_col), bg="white")
        box()
    }

    #---------------------
    # return a list of dates ("date"), spline weights ("spl_weights"), and an object ("x")
    # that mimics Get3YearDataChunk() return type: c(gpp_smooth, gpp_resid, qc)
    return(
        list(date=fluxDT[, date], x=c(fluxDT[, gpp_smooth], fluxDT[, gpp_resids], fluxDT[, qc]), spl_weights=fluxDT[, spl_weight]))
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Wraps up the soup-to-nuts Q2 approach for a given FLUXNET site-year
# extra parameters are passed to GetSmoothFluxYear()
# the plot argument determine whether GetSmoothFluxYear plots it's 24 month data and spline
# addPlot will add the raw data and the smoothed line to an assumed existing plot
# addPlot can be: single T/F and either nothing or everything will plot; or a vector of T/F
# for each of: points, spline, phenometrics
DoFluxPheno <- function(flux_file, year_of_interest, pheno_pars, addPlot=F, gpp_metric="GPP_DT_CUT_REF", weight_var="GPP_DT_CUT_SE", fake_head_tail=FALSE, pt_col=rgb(0.3, 0.3, 0.3), line_col="black", metric_colors = magma(8)[1:7]){
    site_to_do <- gsub("FLX_(.*)_FLUXNET2015_.*.csv$", "\\1", basename(flux_file)) # determine site name
    pheno_period_start <- as.Date(paste(year_of_interest, "-1-1", sep=""))
    pheno_period_end <- as.Date(paste(year_of_interest, "-12-31", sep=""))

    # read the file, get a smooth data chunk, and process the phenology
    flux_x <- GetSmoothFluxYear(flux_file, year_of_interest, gpp_metric=gpp_metric, weight_var=weight_var, spline_spar = pheno_pars$spline_spar, max_span = pheno_pars$max_span, fake_head_tail = fake_head_tail)
    if(!all(is.na(flux_x))){
        flux_metrics <- AnnualPhenologyC6(flux_x$x, dates=flux_x$date, pheno_pars=pheno_pars, pheno_period_start=pheno_period_start, pheno_period_end=pheno_period_end)
        flux_metrics[flux_metrics == 32767] <- NA # b/c Q2 default is to use 32767 for NA values
    }else{
        # we didn't get any phenology for this one, so we store all NA for these
        flux_metrics <- rep(NA, 25)
    }

    # add plot(s) if requested
    if(length(addPlot) == 1 & any(addPlot == TRUE)){
        # assuming there's no existing plot, so we setup the plot
        day_seq <- seq.Date(pheno_period_start, pheno_period_end, by="day")
        if(!all(is.na(flux_x))){
            ylim <- range(flux_x$x[1:(length(flux_x$x) / 3)])
        }else{
            ylim <- c(0, 28) # if no points/line then arbitrary
        }
        plot(day_seq, rep(NA, length(day_seq)), type="n", ylim=ylim, xlab="", ylab=gpp_metric, yaxt="n")
        axis(side=2)
        addPlot <- rep(TRUE, 3) # plot points, spline, and phenometrics
    }

    # add points/line/metrics if requested
    if(any(addPlot)){
        if(!all(is.na(flux_x))){ # check if smoothing was successful; otherwise, nothing to plot
            # add the flux metrics to the graph if available and requested
            if(addPlot[3]){
                if(!all(is.na(flux_metrics))){
                    abline(v=as.Date(flux_metrics[c(5:11, 17:23)], origin="1970-1-1"), col=metric_colors, lwd=1.2, lty=rep(c(1, 2), each=7)) # add the lines to the plot
                }
            }
            # extract components of x if we have to plot points/line
            if(any(addPlot[1:2])){    
                data_length <- length(flux_x$x) / 3
                gpp <- flux_x$x[1:data_length]
                resids <- flux_x$x[(data_length + 1):(2 * data_length)]
                qc <- flux_x$x[(2 * data_length + 1):(3 * data_length)]
            }
            if(addPlot[1]){
                # we add only points in the yoi and those that were permuted (qc=1)
                spl_weights_breaks <- seq(0, 10, by=2)
                cexs <- seq(0.25, 1.25, len=5)
                plot_cexs <- cexs[findInterval(flux_x$spl_weights, spl_weights_breaks, all.inside=T)]
                yoi_inds <- as.integer(strftime(flux_x$date, format="%Y")) == year_of_interest
                yoi_inds[qc == 1] <- TRUE
                pchs <- c(16, 1)[qc[yoi_inds] + 1] # plot qc=0 as filled, qc=1 as open
                points(flux_x$date[yoi_inds], gpp[yoi_inds] - resids[yoi_inds], pch=pchs, col=pt_col, cex=plot_cexs)
            }
            if(addPlot[2]){
                # add the spline fit if requested
                points(flux_x$date, gpp, type="l", lwd=2, col=line_col)
            }
        } # end check for good smoothing result
    } # end plotting

    # pack results into a list for return value
    flux_results <- c(
        list(site_to_do, year_of_interest, gpp_metric),
        as.list(flux_metrics[2] / 10), # Cycle 1 GPParea (is scaled by 10)
        as.list(flux_metrics[3:4]), # Cycle 1 GPPamp, GPPmin
        as.list(as.Date(flux_metrics[5:11], origin="1970-1-1")), # cycle1 DOY metrics
        as.list(flux_metrics[12]), # Cycle 1 overall QA
        as.list(flux_metrics[13]), # Cycle 1 detailed QA
        as.list(flux_metrics[14] / 10), # Cycle 2 GPParea (is scaled by 10)
        as.list(flux_metrics[15:16]), # Cycle 2 GPPamp, GPPmin
        as.list(as.Date(flux_metrics[17:23], origin="1970-1-1")), # cycle1 DOY metrics
        as.list(flux_metrics[24]), # Cycle 2 overall QA
        as.list(flux_metrics[25]) # Cycle 2 detailed QA
    )

    flux_results_names <- c(
        "site_to_do", "year", "gpp_var",
        "flux_gpparea_cycle1",
        "flux_gppamp_cycle1",
        "flux_gppmin_cycle1",
        "flux_gup_cycle1", "flux_midgup_cycle1", "flux_maturity_cycle1", "flux_peak_cycle1", "flux_senescence_cycle1", "flux_midgreendown_cycle1", "flux_dormancy_cycle1",
        "flux_qa_cycle1",
        "flux_qadetailed_cycle1",
        "flux_gpparea_cycle2",
        "flux_gppamp_cycle2",
        "flux_gppmin_cycle2",
        "flux_gup_cycle2", "flux_midgup_cycle2", "flux_maturity_cycle2", "flux_peak_cycle2", "flux_senescence_cycle2", "flux_midgreendown_cycle2", "flux_dormancy_cycle2",
        "flux_qa_cycle2",
        "flux_qadetailed_cycle3"
    )

    names(flux_results) <- flux_results_names

    # return(flux_metrics)
    return(flux_results)
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Run the code over flux sites and produce plots
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#---------------------
# prelims: source the Q2 algorithm functions; get and modify the default parameters
# set some constants and paths
source("/rsstu/users/j/jmgray2/SEAL/hank_code/MCD12Q2C6_AnnualPhenologyFunctions.R")
pheno_pars <- DefaultPhenoParameters()
pheno_pars$min_seg_amplitude <- 4 # just a first guess, need to test
pheno_pars$nbar_scale_factor <- 1 # b/c we're dealing w/ unscaled GPP values
pheno_pars$rel_amp_frac <- 0.35 # back to Q2 default
pheno_pars$nbar_NA_value <- NA # not sure how this is used...
pheno_pars$out_float_scale <- 1 # b/c we don't want to scale the output values
pheno_pars$nbar_NA_value <- NA
pheno_pars$spline_spar <- 0.6 # default is 0.5, bit more smoothing for GPP b/c greater variability
pheno_pars$max_span <- 45 # max number of missing GPP days we allow (unique to FLUXNET implementation)

# gpp_metric <- "GPP_DT_CUT_REF"
# weight_var <- "GPP_DT_CUT_SE"
gpp_metric <- "GPP_DT_VUT_REF"
weight_var <- "GPP_DT_VUT_SE"
time_var <- "TIMESTAMP"
fake_head_tail <- TRUE # if first/last year in series then we permute data

flux_data_dir <- "/rsstu/users/j/jmgray2/SEAL/flux_analysis/data_fluxnet2015/tier1"
out_dir <- "/rsstu/users/j/jmgray2/SEAL/flux_analysis"
out_file <- file.path(out_dir, paste("flux_q2algo_", gpp_metric, ".Rdata", sep=""))
out_plot <- file.path(out_dir, paste("flux_q2algo_", gpp_metric, ".pdf", sep=""))

modis_years <- 2001:2015 # the years in which we'll try and grab phenology
modis_start <- as.Date("2001-1-1")
modis_end <- as.Date("2015-12-31")
gpp_plot_range <- c(0, 30) # ylim for all plots; fixed
pt_col <- rgb(0.4, 0.4, 0.4, alpha=0.5)
# pt_col <- rgb(0.4, 0.4, 0.4)
line_col <- rgb(0.15, 0.15, 0.15, alpha=0.7)
# line_col <- rgb(0.15, 0.15, 0.15)
metric_colors <- c("#000004FF", "#231151FF", "#5F187FFF", "#982D80FF", "#D3436EFF", "#F8765CFF", "#FEBA80FF") # from magma()
plots_per_page <- 4

#---------------------
# collect all flux files and get site names
all_flux_files <- dir(flux_data_dir, pattern=paste(".*_FULLSET_DD_.*.csv$", sep=""), full=T)
all_sites <- gsub("FLX_(.*)_FLUXNET2015_.*.csv$", "\\1", basename(all_flux_files))

# setup the plot
pdf(file=out_plot, width=15, height=10)
par(mar=c(1, 2, 0.5 ,0.5), oma=rep(1, 4), tcl=0.25)

full_results_list <- list()
i <- 1 # site-yoi counter
site_counter <- 1
for(site_to_do in all_sites){
    #-----------------------
    # set up the plot for this site
    # check if we need a new page, and get one w/ layout()
    if(site_counter %% plots_per_page == 1) layout(matrix(1:plots_per_page, ncol=1))

    # set up the blank plot and and annotate w/ background, legend, sitename, etc.
    plot(c(modis_start, modis_end), rep(NA, 2), type="n", ylim=gpp_plot_range, ylab="", xlab="", xaxt="n", yaxt="n")
    trash <- lapply(seq(2001, 2015, by=2), PlotYearBackground) # apply the background
    # axes/labels
    axis(side=1, padj=-1.5, at=seq.Date(modis_start, modis_end, by="year"), labels=strftime(seq.Date(modis_start, modis_end, by="year"), format="%Y"))
    axis(side=2, padj=1.5)
    axis(side=3, at=seq.Date(modis_start, modis_end, by="year"), labels=FALSE)
    axis(side=4, labels=FALSE)
    mtext(gpp_metric, side=2, line=1.5)
    text(par()$usr[2], par()$usr[4], labels=site_to_do, cex=2, adj=c(1.1, 1.3)) # sitename
    box()

    # find the FLUXNET2015 file for this site, and then do pheno extraction for every requested year
    flux_file <- dir(flux_data_dir, pattern=paste(".*", site_to_do, ".*_FULLSET_DD_.*.csv$", sep=""), full=T)
    for(year_of_interest in modis_years){
        print(paste("Doing", site_to_do, "for", year_of_interest)) # DEBUG

        # get the flux metrics, and plot the line and the flux_metrics
        # NOTE: if it weren't for plotting, the code would essentially be these next few lines...
        flux_metrics <- DoFluxPheno(flux_file=flux_file, year_of_interest=year_of_interest, pheno_pars=pheno_pars, gpp_metric=gpp_metric, weight_var=weight_var, fake_head_tail = fake_head_tail, addPlot=c(TRUE, TRUE, TRUE), line_col=line_col, pt_col=pt_col)
        full_results_list[[i]] <- flux_results
        i <- i + 1
    }

    # is it time to make the legend?
    if(site_counter %% plots_per_page == 1){
        legend("topleft", legend=c("splinefit", "obs GPP", "perm GPP", "1st Cycle", "2nd Cycle", "Gup", "MidGup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy"), pch=c(NA, 16, 1, rep(NA, 9)), lwd=c(2, NA, NA, rep(1.5, 9)), lty=c(1, NA, NA, 1, 2, rep(1, 7)), col=c(line_col, pt_col, pt_col, metric_colors[1], metric_colors[1], metric_colors), ncol=2, bg="white")
    }
    site_counter <- site_counter + 1
}
dev.off() # stop plotting

# aggregate all the results and save the file
fullDT <- data.table(do.call(rbind, full_results_list))
save(fullDT, file=out_file)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Analysis: how many site years?
load(file.path(out_dir, "flux_q2algo_GPP_DT_VUT_REF.Rdata"))
fullDT[, has_cycle1 := !all(is.na(.SD)), .SDcols=names(fullDT)[4:15], by=1:nrow(fullDT)]
fullDT[, has_cycle2 := !all(is.na(.SD)), .SDcols=names(fullDT)[16:27], by=1:nrow(fullDT)]
fullDT[, sum(has_cycle1)] # 866
fullDT[, sum(has_cycle2)] # 16
fullDT[, sum(has_cycle1 + has_cycle2)] # 882

load(file.path(out_dir, "flux_q2algo_GPP_DT_CUT_REF.Rdata"))
fullDT[, has_cycle1 := !all(is.na(.SD)), .SDcols=names(fullDT)[4:15], by=1:nrow(fullDT)]
fullDT[, has_cycle2 := !all(is.na(.SD)), .SDcols=names(fullDT)[16:27], by=1:nrow(fullDT)]
fullDT[, sum(has_cycle1)] # 841
fullDT[, sum(has_cycle2)] # 15
fullDT[, sum(has_cycle1 + has_cycle2)] # 856


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# over site years, no plotting
source("/rsstu/users/j/jmgray2/SEAL/hank_code/MCD12Q2C6_AnnualPhenologyFunctions.R")

pheno_pars <- DefaultPhenoParameters()
pheno_pars$min_seg_amplitude <- 4 # just a first guess, need to test
pheno_pars$nbar_scale_factor <- 1 # b/c we're dealing w/ unscaled GPP values
pheno_pars$rel_amp_frac <- 0.35 # back to Q2 default
pheno_pars$nbar_NA_value <- NA # not sure how this is used...
pheno_pars$out_float_scale <- 1 # b/c we don't want to scale the output values
pheno_pars$nbar_NA_value <- NA
pheno_pars$spline_spar <- 0.6 # default is 0.5, bit more smoothing for GPP b/c greater variability
pheno_pars$max_span <- 45 # max number of missing GPP days we allow (unique to FLUXNET implementation)

# gpp_metric <- "GPP_DT_CUT_REF"
# weight_var <- "GPP_DT_CUT_SE"
gpp_metric <- "GPP_DT_VUT_REF"
weight_var <- "GPP_DT_VUT_SE"
time_var <- "TIMESTAMP"
fake_head_tail <- TRUE # if first/last year in series then we permute data

flux_data_dir <- "/rsstu/users/j/jmgray2/SEAL/flux_analysis/data_fluxnet2015/tier1"
out_dir <- "/rsstu/users/j/jmgray2/SEAL/flux_analysis"
out_file <- file.path(out_dir, paste("flux_q2algo_", gpp_metric, "TRASHME_TESTING", ".Rdata", sep=""))

#---------------------
# collect all flux files and get site names
all_flux_files <- dir(flux_data_dir, pattern=paste(".*_FULLSET_DD_.*.csv$", sep=""), full=T)
all_sites <- gsub("FLX_(.*)_FLUXNET2015_.*.csv$", "\\1", basename(all_flux_files))
all_years <- 2001:2014
all_site_years <- expand.grid(flux_file=all_flux_files, year=all_years)
all_site_years$flux_file <- as.character(all_site_years$flux_file)

# just went with a big for loop for the sake of coding time vs. run time, it only takes a few minutes per metric
# however, using an apply() structure with a parallel environment would be much faster
full_results_list <- vector("list", nrow(all_site_years))
for(i in 1:nrow(all_site_years)){
    print(paste("Doing", basename(all_site_years[i, 1]), "for", all_site_years[i, 2]))
    full_results_list[[i]] <- DoFluxPheno(flux_file=all_site_years[i, 1], year_of_interest=all_site_years[i, 2], pheno_pars=pheno_pars, gpp_metric=gpp_metric, weight_var=weight_var, fake_head_tail = fake_head_tail, addPlot=FALSE, line_col=1, pt_col=1)
}
fullDT <- data.table(do.call(rbind, full_results_list))
save(fullDT, file=out_file)

their_mets <- c(127, 158, 189, 255, 267, 276, 286) - (365/2)
our_mets <- as.integer(as.Date(c("2008-10-22", "2008-12-07", "2009-01-18", "2009-03-03", "2009-03-19", "2009-04-06", "2009-04-25")) - as.Date("2009-1-1"))
plot(our_mets, their_mets, xlab="SEAL fluxmetrics", ylab="Asko fluxmetrics")
abline(a=0, b=1, lty=2, col=2)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Linqing Yang & Asko Noormets data
library(data.table)

# load our data: fullDT
fluxmetrics_data_dir <- "/rsstu/users/j/jmgray2/SEAL/flux_analysis"
load(file.path(fluxmetrics_data_dir, "flux_q2algo_GPP_DT_VUT_REF.Rdata"))
# WTF is going on here?
fullDT <- fullDT[, lapply(.SD, unlist), .SDcols=names(fullDT)]
setnames(fullDT, old="site_to_do", new="site")
setkey(fullDT, site, year, gpp_var)

# I can't quickly figure out how to do this w/ lapply and .SDcols, so...onward w/ brute force
fullDT[, DOY_flux_gup_cycle1 := as.integer(strftime(as.Date(flux_gup_cycle1, origin="1970-1-1"), format="%j"))]
fullDT[, DOY_flux_midgup_cycle1 := as.integer(strftime(as.Date(flux_midgup_cycle1, origin="1970-1-1"), format="%j"))]
fullDT[, DOY_flux_maturity_cycle1 := as.integer(strftime(as.Date(flux_maturity_cycle1, origin="1970-1-1"), format="%j"))]
fullDT[, DOY_flux_peak_cycle1 := as.integer(strftime(as.Date(flux_peak_cycle1, origin="1970-1-1"), format="%j"))]
fullDT[, DOY_flux_senescence_cycle1 := as.integer(strftime(as.Date(flux_senescence_cycle1, origin="1970-1-1"), format="%j"))]
fullDT[, DOY_flux_midgreendown_cycle1 := as.integer(strftime(as.Date(flux_midgreendown_cycle1, origin="1970-1-1"), format="%j"))]
fullDT[, DOY_flux_dormancy_cycle1 := as.integer(strftime(as.Date(flux_dormancy_cycle1, origin="1970-1-1"), format="%j"))]


# read in all of the Yang & Noormets data
yang_seasonality_data_dir <- "/rsstu/users/j/jmgray2/SEAL/flux_analysis/Flux_Seasonality_Metrics Database_1992-2014/Flux Seasonality Metrics Database_1992-2014/GPP_DT_VUT_REF"
yang_files <- dir(yang_seasonality_data_dir, pattern=".*metrics.*.txt$", full=T)
yangDT <- do.call(rbind, lapply(yang_files, fread))
yang_names <- names(yangDT)
yang_names[c(1, 5, 6)] <- c("site", "gpp_var", "year")
setnames(yangDT, yang_names)
setkey(yangDT, site, year, gpp_var)

# I can't quickly figure out how to do this w/ lapply and .SDcols, so...onward w/ brute force
yangDT[, DateSFD := as.Date(paste(year, DOY_SFD, sep="-"), format="%Y-%j")]
yangDT[, DateMFD := as.Date(paste(year, DOY_MFD, sep="-"), format="%Y-%j")]
yangDT[, DateEFD := as.Date(paste(year, DOY_EFD, sep="-"), format="%Y-%j")]
yangDT[, DateFmax := as.Date(paste(year, DOY_Fmax, sep="-"), format="%Y-%j")]
yangDT[, DateSFR := as.Date(paste(year, DOY_SFR, sep="-"), format="%Y-%j")]
yangDT[, DateMFR := as.Date(paste(year, DOY_MFR, sep="-"), format="%Y-%j")]
yangDT[, DateEFR := as.Date(paste(year, DOY_EFR, sep="-"), format="%Y-%j")]

# merge our data
bothDT <- merge(yangDT, fullDT)
# make a column for whether or not SEAL flux metrics had more than 1 cycle
bothDT[, multicycle := !is.na(flux_gup_cycle2)]
# make a column for northern (1) vs southern (2) hemispheres; for plotting
bothDT[, SH := LAT < 0]

# calculate corrected version sof Yang SH points
yangDT[, DOY_SFD := 
bothDT[, DOY_MFD_cor := DOY_MFD]
bothDT[LAT < 0, DOY_MFD_cor := DOY_MFD_cor + 182.5]
yangDT[, DOY_EFD := 
yangDT[, DOY_Fmax :=
yangDT[, DOY_SFR := 
yangDT[, DOY_MFR := 
yangDT[, DOY_EFR := 

plot(bothDT[, .(DOY_flux_gup_cycle1, DOY_SFD)])
points(bothDT[(SH), .(DOY_flux_gup_cycle1, DOY_SFD)], col=2)
points(bothDT[(multicycle), .(DOY_flux_gup_cycle1, DOY_SFD)], col=4)

plot(bothDT[, .(DOY_flux_midgup_cycle1, DOY_MFD)], xlab="SEAL FLUXNET MidGup Cycle 1", ylab="Yang & Noormets DOY_SFD")
points(bothDT[(SH), .(DOY_flux_midgup_cycle1, DOY_MFD)], col=2)
points(bothDT[(multicycle), .(DOY_flux_midgup_cycle1, DOY_MFD)], col=4)
legend("topright", legend=c("Southern Hemisphere", "Multiple cycles in Year"), pch=1, col=c(2, 4))

plot(bothDT[, .(DOY_flux_midgup_cycle1, DOY_MFD_cor)], xlab="SEAL FLUXNET MidGup Cycle 1", ylab="Yang & Noormets DOY_SFD, corrected SH")
points(bothDT[(SH), .(DOY_flux_midgup_cycle1, DOY_MFD_cor)], col=2)
points(bothDT[(multicycle), .(DOY_flux_midgup_cycle1, DOY_MFD_cor)], col=4)
legend("topright", legend=c("Southern Hemisphere", "Multiple cycles in Year"), pch=1, col=c(2, 4))
abline(a=0,b=1)

PasteDate <- function(doy, year) as.Date(paste(year, doy, sep="-"), format="%Y-%j")

site <- "AU-DaP"
year_to_compare <- 2009

# load the Yang seasonality data
yangDT <- fread(dir(yang_seasonality_data_dir, pattern=paste(site, ".*metrics.*.txt$", sep=""), full=T))

# compare our metrics to Yang & Noormets
yang_metrics <- yangDT[Year == year_to_compare, lapply(.SD, PasteDate, year=year_to_compare), .SDcols=names(DT)[c(11, 8, 12, 20, 17, 9, 18)]]
seal_metrics <- fullDT[site_to_do=="AU-DaP" & year == 2009, .SD, .SDcol=names(fullDT)[7:13]]
plot(seal_metrics, yang_metrics)