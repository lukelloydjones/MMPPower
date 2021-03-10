# =================================================
# Increasing sites versus increasing sampling 
# density within site
# Author: Luke Lloyd-Jones
# Date started: 09/11/2020
# Date updated: 02/12/2020
# =================================================
args <- commandArgs(trailingOnly = TRUE)
# args[1]: NRM name - "burdekin", "russ_mull", "tully", "wtsndys"
# args[2]: analyte - "CHL_QAQC" "SECCHI_DEPTH" "PP_QAQC" "PN_SHIM_QAQC" "SS_QAQC" 
args <- c("tully", "PP_QAQC")
# Load some libraries
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
#library(MMP)
library(lubridate)
library(dplyr)
library(ggplot2)
library(mgcv)
library(car)
source("rscripts/power_boot_functions.R")
# File output path setup
# ----------------------
# args        <- c("russ_mull", "SECCHI_DEPTH", "post")
top.pth      <- "sites_vs_samples"
nrm.nm       <- as.character(args[1])
analyte      <- as.character(args[2])
res.nrm.pth  <- paste0(nrm.nm, "_phase_3_analysis")
figs.pth     <- "figures"
mdls.pth     <- "model_fits"
power.st.pth <- "power_sets"
ifelse(!dir.exists(file.path(top.pth)), dir.create(file.path(top.pth)), FALSE)
ifelse(!dir.exists(file.path(top.pth, res.nrm.pth)), 
        dir.create(file.path(top.pth, res.nrm.pth)), FALSE)
ifelse(!dir.exists(file.path(top.pth, res.nrm.pth, figs.pth)), 
        dir.create(file.path(top.pth, res.nrm.pth, figs.pth)), FALSE)
ifelse(!dir.exists(file.path(top.pth, res.nrm.pth, mdls.pth)), 
        dir.create(file.path(top.pth, res.nrm.pth, mdls.pth)), FALSE)
ifelse(!dir.exists(file.path(top.pth, res.nrm.pth, power.st.pth)), 
        dir.create(file.path(top.pth, res.nrm.pth, power.st.pth)), FALSE)
# ----------------------------------------------------------
# Special function to choose the minimum value to add to the 
# zeroes of the logged NOX distribution based on the qqplot 
# ----------------------------------------------------------
# Load the nrmekin data
nrm <- read.csv(paste0("data/raw/MMPDataset-sans-loggers-allfields_", nrm.nm,"_subset.csv"))
# How many analytes
analytes <- nrm[, grep("_QAQC", colnames(nrm))]
# 17 analytes and 1133 measurements 
analyte.nms <- data.frame(t(matrix(c("DIP", "Dissolved inorganic phosphorus" 
                                     , " DOC ",  "Dissolved organic carbon"
                                     , " NH4 ",  "Ammonium"
                                     , " NO2 ",  "Nitrite"
                                     , " NO3 ",  "Nitrate"
                                     , " PN_SHIM ",  "Particulate nitrogen-AIMS-Shimadzu "
                                     , " POC ",      "Particulate organic carbon"
                                     , " PP ",       "Particulate phosphorus"
                                     , " SAL ",      "Salinity"
                                     , " SI ",       "Dissolved silica"
                                     , " SS ",       "Total suspended solids-gravimetrically"
                                     , " TDN_PER ",  "Total dissolved nitrogen-persulphate digestion"
                                     , " TDN_SHIM ", "Total dissolved nitrogen-Shimadzu"
                                     , " TDP_PER ",  "Total dissolved phosphorus-persulphate digestion"
                                     , " TEMP ",     "Temperature"
                                     , " CHL ",      "Chlorophyll a measured via filtration and fluorescence"
                                     , " PHAEO ",    "Phaeophytin a measured via filtration and fluorescence"
                                     , " SECCHI_DEPTH ",  "Secchi depth as a measurement of water clarity"), nrow = 2)))
# Collection dates
dates    <- as.Date(nrm$COLLECTION_START_DATE, "%d/%m/%y")
nrm$Date <- dates
# ==========================================
# Subset the data to the analyte of interest
# ==========================================
#"CHL_QAQC", "SECCHI_DEPTH", "PP_QAQC",
# Make the analyte specific path
ifelse(!dir.exists(file.path(top.pth, res.nrm.pth, figs.pth, analyte)),
        dir.create(file.path(top.pth, res.nrm.pth, figs.pth, analyte)), FALSE)
ifelse(!dir.exists(file.path(top.pth, res.nrm.pth, mdls.pth, analyte)),
        dir.create(file.path(top.pth, res.nrm.pth, mdls.pth, analyte)), FALSE)
ifelse(!dir.exists(file.path(top.pth, res.nrm.pth, power.st.pth, analyte)),
        dir.create(file.path(top.pth, res.nrm.pth, power.st.pth, analyte)), FALSE)
fig.pth <- file.path(top.pth, res.nrm.pth, figs.pth, analyte)
mdl.pth <- file.path(top.pth, res.nrm.pth, mdls.pth, analyte)
pwr.pth <- file.path(top.pth, res.nrm.pth, power.st.pth, analyte)
# -----------------------------
# Grab the data for the analyte
# -----------------------------
nrm.analyte <- nrm[, c("STATION_NAME", "DEPTH_CODE", "DUPLICATE",
                       "DEPTH", analyte, "LATITUDE", "LONGITUDE",
                       "ACOUSTIC_DEPTH", "SHORT_NAME", "PROJECT", "Date")]
# Rename the column to make downstream more general
colnames(nrm.analyte)[grep(analyte, colnames(nrm.analyte))] <- "ANALYTE"
# ------------------------------------------------------------
# What to do with duplicates? Let's just average them for now?
# ------------------------------------------------------------
nrm.analyte$Date <- as.Date(nrm.analyte$Date, "%Y-%m-%d")
nrm.analyte.avg  <- as.data.frame(nrm.analyte %>%
                                  group_by(STATION_NAME, PROJECT, DEPTH_CODE, Date) %>%
                                  mutate(ANALYTE_MN = mean(ANALYTE, na.rm = T)))
# Change NaN to NAs. WHy NaNs? No idea.
nrm.analyte.avg$ANALYTE_MN[is.nan(nrm.analyte.avg$ANALYTE_MN)] <- NA
head(nrm.analyte.avg)
dim(nrm.analyte.avg)
# -----------------------------------------
# Subset to just the unique average values.
# Could they be the same by chance I
# suppose the could Oh god.
# -----------------------------------------
pull.vars <- unique(nrm.analyte.avg[, c("STATION_NAME", "PROJECT", "Date", "DEPTH_CODE")])
for (i in seq(1, dim(pull.vars)[1]))
{
  print(paste0("Doing ", i))
  nrm.analyte.avg.i <- filter(nrm.analyte.avg,
                              STATION_NAME == pull.vars[i, 1] &
                                PROJECT      == pull.vars[i, 2] &
                                Date         == pull.vars[i, 3] &
                                DEPTH_CODE   == pull.vars[i, 4])
  if (dim(nrm.analyte.avg.i)[1] > 1)
  {
    nrm.analyte.avg.i.nd <- nrm.analyte.avg.i[-which(duplicated(nrm.analyte.avg.i$ANALYTE_MN)), ]
  } else {
    nrm.analyte.avg.i.nd <- nrm.analyte.avg.i
  }
  print(paste0("Dim is ", dim(nrm.analyte.avg.i.nd)[1]))
  if (i == 1)
  {
    nrm.analyte.avg.red <- nrm.analyte.avg.i.nd
  } else {
    nrm.analyte.avg.red <- rbind(nrm.analyte.avg.red, nrm.analyte.avg.i.nd)
  }
}
dim(nrm.analyte.avg.red)
# --------------------------------------------
# Average over the depths
#  - Something interesting to think about here
#    what depth averaging means
# --------------------------------------------
nrm.analyte.avg.red.mn.dpth  <- as.data.frame(nrm.analyte.avg.red %>%
                                              group_by(STATION_NAME, PROJECT, Date) %>%
                                              mutate(ANALYTE_MN_DPTH = mean(ANALYTE_MN, na.rm = T)))
# Not really sure about the weighted mean?? Out of two measurement why should the deeper one have more relavance.
pull.vars <- unique(nrm.analyte.avg.red[, c("STATION_NAME", "PROJECT", "Date")])
for (i in seq(1, dim(pull.vars)[1]))
{
  print(paste0("Doing ", i))
  nrm.analyte.avg.i <- filter(nrm.analyte.avg.red.mn.dpth,
                              STATION_NAME == pull.vars[i, 1] &
                                PROJECT      == pull.vars[i, 2] &
                                Date         == pull.vars[i, 3])
  if (dim(nrm.analyte.avg.i)[1] > 1)
  {
    nrm.analyte.avg.i.nd <- nrm.analyte.avg.i[-which(duplicated(nrm.analyte.avg.i$ANALYTE_MN_DPTH)), ]
  } else {
    nrm.analyte.avg.i.nd <- nrm.analyte.avg.i
  }
  if (i == 1)
  {
    nrm.analyte.avg.red.mn.dpth.2 <- nrm.analyte.avg.i.nd
  } else {
    nrm.analyte.avg.red.mn.dpth.2 <- rbind(nrm.analyte.avg.red.mn.dpth.2, nrm.analyte.avg.i.nd)
  }
}
# ------------------------------
# Process zeros if present
# ------------------------------
ana.vals <- nrm.analyte.avg.red.mn.dpth.2$ANALYTE_MN_DPTH
zeros    <- which(ana.vals == 0)
if (length(which(ana.vals == 0)) > 0)
{
  nrm.analyte.avg.red.mn.dpth.2$ANALYTE_MN_DPTH[zeros] <- NA
}
# ------------------------------
# Separate pre-2015 to post 2015
# ------------------------------
# Compute some common variables
nrm.analyte.avg.red.mn.dpth.2$PROJECT    <- as.factor(nrm.analyte.avg.red.mn.dpth.2$PROJECT)
nrm.analyte.avg.red.mn.dpth.2$SHORT_NAME <- as.factor(nrm.analyte.avg.red.mn.dpth.2$SHORT_NAME)
# Separate into pre and post
# --------------------------
nrm.analyte.pre.2015 <- filter(nrm.analyte.avg.red.mn.dpth.2,
                               Date < as.Date("2015-01-01"))
nrm.analyte.pst.2015 <- filter(nrm.analyte.avg.red.mn.dpth.2,
                               Date >= as.Date("2015-01-01"))
max(nrm.analyte.pst.2015$Date) - min(nrm.analyte.pst.2015$Date)
# --------------------------------------------------------------
# Down sample pst-2015 to
#  - same sites and site visits per year as pre-2015
#  - add in extra sites but same sampling frequency
#  - add in increased sampling frequency but same site as
#    pre-2015
#  Run the regional model
# --------------------------------------------------------------
# What was the sample frequency pre-2015
pre.sts    <- colSums(table(format(nrm.analyte.pre.2015$Date, "%Y-%m"), 
                                   nrm.analyte.pre.2015$SHORT_NAME))
pre.st.nms <- names(which(pre.sts > 0))
# Downsample to nrm1, nrm2, and nrm4 on the 2, 6, 10 months
if (nrm.nm == "burdekin")
{
  nrm.analyte.pst.2015.pre.sts <- filter(nrm.analyte.pst.2015,
                                         SHORT_NAME == "BUR1" |
                                         SHORT_NAME == "BUR2" |
                                         SHORT_NAME == "BUR4")
} else if (nrm.nm == "tully") {
  nrm.analyte.pst.2015.pre.sts <- filter(nrm.analyte.pst.2015,
                                         SHORT_NAME == "TUL3")
} else if (nrm.nm == "russ_mull") {
  nrm.analyte.pst.2015.pre.sts <- filter(nrm.analyte.pst.2015,
                                         SHORT_NAME == "RM1" |
                                         SHORT_NAME == "RM7" |
                                         SHORT_NAME == "RM8")
} else if (nrm.nm == "wtsndys") {
  # There were monitoring at 3 sites prior but ceased one. We are 
  # going to keep on new one post.
  nrm.analyte.pst.2015.pre.sts <- filter(nrm.analyte.pst.2015,
                                         SHORT_NAME == "WHI1" |
                                         SHORT_NAME == "WHI4" |
                                         SHORT_NAME == "WHI5")
}
# Look at the sampling density
table(format(nrm.analyte.pre.2015$Date, "%Y-%m"), nrm.analyte.pre.2015$SHORT_NAME)
table(format(nrm.analyte.pst.2015$Date, "%Y-%m"), nrm.analyte.pst.2015$SHORT_NAME)
# Take the months needed 
if (nrm.nm == "burdekin")
{
  yr.mths <- c("2015-03", "2015-05", "2015-06", "2015-09",
               "2016-02", "2016-06", "2016-09",
               "2017-02", "2017-05", "2017-09",
               "2018-03", "2018-05", "2018-06", "2018-09",
               "2019-03", "2019-05", "2019-06")
} else if (nrm.nm == "tully") {
  yr.mths <- c("2015-03", "2015-06", "2015-07", "2015-09",
               "2016-03", "2016-07", "2016-09",
               "2017-03", "2017-05", "2017-09",
               "2018-02", "2018-07", "2018-10",
               "2019-02", "2019-05")
} else if (nrm.nm == "russ_mull") {
  yr.mths <- c("2015-03", "2015-07", "2015-10",
               "2016-03", "2016-07", "2016-10",
               "2017-03", "2017-05", "2017-09",
               "2018-02", "2018-05", "2018-10",
               "2019-02", "2019-05")
} else if (nrm.nm == "wtsndys") {
  yr.mths <- c("2015-03", "2015-07", "2015-10",
               "2016-03", "2016-06", "2016-09",
               "2017-03", "2017-05", "2017-09",
               "2018-02", "2018-06", "2018-09",
               "2019-02", "2019-05")
}
nrm.analyte.pst.2015$YR_MNTH   <- format(nrm.analyte.pst.2015$Date, "%Y-%m")
nrm.analyte.pst.2015.sub.samps <- nrm.analyte.pst.2015[which(nrm.analyte.pst.2015$YR_MNTH %in% yr.mths), ]
table(format(nrm.analyte.pst.2015.sub.samps$Date, "%Y-%m"), nrm.analyte.pst.2015.sub.samps$SHORT_NAME)
table(nrm.analyte.pst.2015.sub.samps$Date, nrm.analyte.pst.2015.sub.samps$SHORT_NAME)
# Dates to remove as we get double ups
if (nrm.nm == "burdekin")
{
  dts.rm <- as.Date(c("2015-06-15", "2015-06-18", "2016-06-01", "2018-05-15", "2018-06-21", "2019-05-07"))
  nrm.analyte.pst.2015.sub.samps <- nrm.analyte.pst.2015.sub.samps[-which(nrm.analyte.pst.2015.sub.samps$Date %in% dts.rm), ]
} else if (nrm.nm == "tully") {
  # No duplicates for Tully
  nrm.analyte.pst.2015.sub.samps <- nrm.analyte.pst.2015.sub.samps
} else if (nrm.nm == "russ_mull") {
  dts.rm <- as.Date(c("2017-03-23"))
  nrm.analyte.pst.2015.sub.samps <- nrm.analyte.pst.2015.sub.samps[-which(nrm.analyte.pst.2015.sub.samps$Date %in% dts.rm), ]
} else if (nrm.nm == "wtsndys") {
  dts.rm <- as.Date(c("2015-03-07", "2019-02-01", "2019-02-02"))
  nrm.analyte.pst.2015.sub.samps <- nrm.analyte.pst.2015.sub.samps[-which(nrm.analyte.pst.2015.sub.samps$Date %in% dts.rm), ]
} 
table(format(nrm.analyte.pst.2015.sub.samps$Date, "%Y-%m"), nrm.analyte.pst.2015.sub.samps$SHORT_NAME)
# Subsample this one to just three site
if (nrm.nm == "burdekin")
{
  nrm.analyte.pst.2015.sub.samps.pre.sts <- filter(nrm.analyte.pst.2015.sub.samps,
                                                   SHORT_NAME == "BUR1" |
                                                   SHORT_NAME == "BUR2" |
                                                   SHORT_NAME == "BUR4")
} else if (nrm.nm == "tully") {
  nrm.analyte.pst.2015.sub.samps.pre.sts <- filter(nrm.analyte.pst.2015.sub.samps,
                                                   SHORT_NAME == "TUL3")
} else if (nrm.nm == "russ_mull") {
  nrm.analyte.pst.2015.sub.samps.pre.sts <- filter(nrm.analyte.pst.2015.sub.samps,
                                                   SHORT_NAME == "RM1" |
                                                   SHORT_NAME == "RM7" |
                                                   SHORT_NAME == "RM8")
} else if (nrm.nm == "wtsndys") {
  nrm.analyte.pst.2015.sub.samps.pre.sts <- filter(nrm.analyte.pst.2015.sub.samps,
                                                   SHORT_NAME == "WHI1" |
                                                   SHORT_NAME == "WHI4" |
                                                   SHORT_NAME == "WHI5")
}
table(format(nrm.analyte.pst.2015.sub.samps.pre.sts$Date, "%Y-%m"), nrm.analyte.pst.2015.sub.samps.pre.sts$SHORT_NAME)
# The data are
# NULL -
dim(nrm.analyte.pst.2015.sub.samps.pre.sts)
# NULL + increased sites
dim(nrm.analyte.pst.2015.sub.samps)
# NULL + increased sampling
dim(nrm.analyte.pst.2015.pre.sts)
# All
dim(nrm.analyte.pst.2015)
# --------------------------------------------------------------------
# Visualise each one to make sure it what I think
# --------------------------------------------------------------------
p.null <-  ggplot(nrm.analyte.pst.2015.sub.samps.pre.sts,
                  aes_string(x = "Date",
                             y = "ANALYTE_MN_DPTH")) +
                  geom_smooth(method="lm", se = F) +
                  facet_wrap(.~SHORT_NAME, ncol = 2) +
                  geom_point() + ylab(analyte) +
                  xlab("Year")
png(filename =  paste0(fig.pth, "/", analyte,"_site_vs_samples_null_set.png"),
    bg = "white", width =  10, height = 6, units = 'in', res = 300)
print(p.null)
dev.off()
p.add.sites <-  ggplot(nrm.analyte.pst.2015.sub.samps,
                       aes_string(x = "Date",
                                  y = "ANALYTE_MN_DPTH")) +
                       geom_smooth(method="lm", se = F) +
                       facet_wrap(.~SHORT_NAME, ncol = 2) +
                       geom_point() + ylab(analyte) +
                       xlab("Year")
png(filename = paste0(fig.pth, "/", analyte,"_site_vs_samples_null_plus_sites.png"),
    bg = "white", width =  10, height = 6, units = 'in', res = 300)
print(p.add.sites)
dev.off()
p.add.samples <-  ggplot(nrm.analyte.pst.2015.pre.sts,
                         aes_string(x = "Date",
                                    y = "ANALYTE_MN_DPTH")) +
                         geom_smooth(method="lm", se = F) +
                         facet_wrap(.~SHORT_NAME, ncol = 2) +
                         geom_point() + ylab(analyte) +
                         xlab("Year")
png(filename = paste0(fig.pth, "/", analyte, "_site_vs_samples_null_plus_samples.png"),
    bg = "white", width =  10, height = 6, units = 'in', res = 300)
print(p.add.samples)
dev.off()

# --------------------------------------------------------------------
# For each of the sets compute the minimum date and the seasonal terms
# relative to a Jan 1st date from the first year
# Pst all
nrm.analyte.pst.2015$Date_2  <- as.numeric(nrm.analyte.pst.2015$Date - min(nrm.analyte.pst.2015$Date))
nrm.analyte.pst.2015$Date_3  <- as.numeric(nrm.analyte.pst.2015$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pst.2015$Date)),"-01-01")))
nrm.analyte.pst.2015$SEAS_C1 <- cos(2 * pi * nrm.analyte.pst.2015$Date_3 / 365.25)
nrm.analyte.pst.2015$SEAS_C2 <- sin(2 * pi * nrm.analyte.pst.2015$Date_3 / 365.25)
# Null
nrm.analyte.pst.2015.sub.samps.pre.sts$Date_2  <- as.numeric(nrm.analyte.pst.2015.sub.samps.pre.sts$Date - min(nrm.analyte.pst.2015.sub.samps.pre.sts$Date))
nrm.analyte.pst.2015.sub.samps.pre.sts$Date_3  <- as.numeric(nrm.analyte.pst.2015.sub.samps.pre.sts$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pst.2015.sub.samps.pre.sts$Date)),"-01-01")))
nrm.analyte.pst.2015.sub.samps.pre.sts$SEAS_C1 <- cos(2 * pi * nrm.analyte.pst.2015.sub.samps.pre.sts$Date_3 / 365.25)
nrm.analyte.pst.2015.sub.samps.pre.sts$SEAS_C2 <- sin(2 * pi * nrm.analyte.pst.2015.sub.samps.pre.sts$Date_3 / 365.25)
# Null + sites
nrm.analyte.pst.2015.sub.samps$Date_2  <- as.numeric(nrm.analyte.pst.2015.sub.samps$Date - min(nrm.analyte.pst.2015.sub.samps$Date))
nrm.analyte.pst.2015.sub.samps$Date_3  <- as.numeric(nrm.analyte.pst.2015.sub.samps$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pst.2015.sub.samps$Date)),"-01-01")))
nrm.analyte.pst.2015.sub.samps$SEAS_C1 <- cos(2 * pi * nrm.analyte.pst.2015.sub.samps$Date_3 / 365.25)
nrm.analyte.pst.2015.sub.samps$SEAS_C2 <- sin(2 * pi * nrm.analyte.pst.2015.sub.samps$Date_3 / 365.25)
# Null + samples
nrm.analyte.pst.2015.pre.sts$Date_2  <- as.numeric(nrm.analyte.pst.2015.pre.sts$Date - min(nrm.analyte.pst.2015.pre.sts$Date))
nrm.analyte.pst.2015.pre.sts$Date_3  <- as.numeric(nrm.analyte.pst.2015.pre.sts$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pst.2015.pre.sts$Date)),"-01-01")))
nrm.analyte.pst.2015.pre.sts$SEAS_C1 <- cos(2 * pi * nrm.analyte.pst.2015.pre.sts$Date_3 / 365.25)
nrm.analyte.pst.2015.pre.sts$SEAS_C2 <- sin(2 * pi * nrm.analyte.pst.2015.pre.sts$Date_3 / 365.25)
# ===========================================
# Region scale
#  - Run the seasonal linear models
# ===========================================
# --------
# NULL
# --------
if ((length(which(nrm.analyte.pst.2015.sub.samps.pre.sts$PROJECT == "MMP-JCU")) > 0) & nrm.nm != "tully" & analyte != "PN_SHIM_QAQC")
{
  fit.nrm.null               <- lm(log(ANALYTE_MN_DPTH) ~ PROJECT + SHORT_NAME  + Date_2 + SEAS_C1 + SEAS_C2,
                                   data = nrm.analyte.pst.2015.sub.samps.pre.sts)
} else if ((length(which(nrm.analyte.pst.2015.sub.samps.pre.sts$PROJECT == "MMP-JCU")) > 0) & nrm.nm == "tully" & analyte != "PN_SHIM_QAQC" & analyte != "PP_QAQC") {
  # Only one site for tully
  fit.nrm.null               <- lm(log(ANALYTE_MN_DPTH) ~ PROJECT + Date_2 + SEAS_C1 + SEAS_C2,
                                   data = nrm.analyte.pst.2015.sub.samps.pre.sts)
} else {
  fit.nrm.null               <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2,
                                   data = nrm.analyte.pst.2015.sub.samps.pre.sts)
}
summary(fit.nrm.null)
png(filename = paste0(fig.pth, "/regional_", analyte, "_null_qqplot_log_lm.png"), 
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
qqPlot(fit.nrm.null$residuals)
dev.off()
# Null + sites
# ------------
if ((length(which(nrm.analyte.pst.2015.sub.samps$PROJECT == "MMP-JCU")) > 0) & analyte != "PN_SHIM_QAQC")
{
  fit.nrm.null.sts               <- lm(log(ANALYTE_MN_DPTH) ~ PROJECT + SHORT_NAME  + Date_2 + SEAS_C1 + SEAS_C2,
                                       data = nrm.analyte.pst.2015.sub.samps)
} else {
  fit.nrm.null.sts               <- lm(log(ANALYTE_MN_DPTH) ~  SHORT_NAME  + Date_2 + SEAS_C1 + SEAS_C2,
                                       data = nrm.analyte.pst.2015.sub.samps)
}
summary(fit.nrm.null.sts)
png(filename = paste0(fig.pth, "/regional_", analyte, "_null_pls_sites_qqplot_log_lm.png"), 
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
qqPlot(fit.nrm.null.sts$residuals)
dev.off()
# Null + samples
# ------------
# nrm.analyte.pst.2015.pre.sts <- nrm.analyte.pst.2015.pre.sts[-83, ]
if ((length(which(nrm.analyte.pst.2015.pre.sts$PROJECT == "MMP-JCU")) > 0) & nrm.nm != "tully" & analyte != "PN_SHIM_QAQC")
{
  fit.nrm.null.smps                <- lm(log(ANALYTE_MN_DPTH) ~ PROJECT + SHORT_NAME  + Date_2 + SEAS_C1 + SEAS_C2,
                                         data = nrm.analyte.pst.2015.pre.sts)
} else if ((length(which(nrm.analyte.pst.2015.pre.sts$PROJECT == "MMP-JCU")) > 0) & nrm.nm == "tully" & analyte != "PN_SHIM_QAQC") {
  # Only one site for tully
  fit.nrm.null.smps                <- lm(log(ANALYTE_MN_DPTH) ~ PROJECT + Date_2 + SEAS_C1 + SEAS_C2,
                                        data = nrm.analyte.pst.2015.pre.sts)
} else if (nrm.nm == "tully" & analyte == "PN_SHIM_QAQC") {
  # Only one site for tully
  fit.nrm.null.smps                <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2,
                                         data = nrm.analyte.pst.2015.pre.sts)
} else {
  fit.nrm.null.smps                <- lm(log(ANALYTE_MN_DPTH) ~ SHORT_NAME  + Date_2 + SEAS_C1 + SEAS_C2,
                                         data = nrm.analyte.pst.2015.pre.sts)
}
summary(fit.nrm.null.smps)
png(filename = paste0(fig.pth, "/regional_", analyte, "_null_pls_samples_qqplot_log_lm.png"), 
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
qqPlot(fit.nrm.null.smps$residuals)
dev.off()
# All
# ----------
if (nrm.nm != "wtsndys")
{
  fit.nrm.pst                <- lm(log(ANALYTE_MN_DPTH) ~ PROJECT + SHORT_NAME  + Date_2 + SEAS_C1 + SEAS_C2,
                                   data = nrm.analyte.pst.2015)
} else {
  nrm.analyte.pst.2015       <- nrm.analyte.pst.2015[-which(nrm.analyte.pst.2015$SHORT_NAME == "WHI0"), ]
  fit.nrm.pst                <- lm(log(ANALYTE_MN_DPTH) ~ SHORT_NAME  + Date_2 + SEAS_C1 + SEAS_C2,
                                   data = nrm.analyte.pst.2015)
}
summary(fit.nrm.pst)
png(filename = paste0(fig.pth, "/regional_", analyte, "_all_samples_qqplot_log_lm.png"), 
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
qqPlot(fit.nrm.pst$residuals)
dev.off()
# ------------------------------------------
# Write out models and informative model output
# ------------------------------------------
save(fit.nrm.null,      file = paste0(mdl.pth, "/", analyte, "_site_vs_samples_null.Rdata"))
save(fit.nrm.null.sts,  file = paste0(mdl.pth, "/", analyte, "_site_vs_samples_null_plus_sites.Rdata"))
save(fit.nrm.null.smps, file = paste0(mdl.pth, "/", analyte, "_site_vs_samples_nul_plus_samples.Rdata"))
save(fit.nrm.pst,       file = paste0(mdl.pth, "/", analyte, "_site_vs_samples_all.Rdata"))
# -------------------------
# Let's do the power
# -------------------------
deltas        <- seq(-0.2, 0, 0.02)
nSim          <- 200
nrm.null      <- sapply(deltas, function(x) powerBoot2Regional(fit.nrm.null,      nSim, x, "Date_2", 1))
nrm.null.sts  <- sapply(deltas, function(x) powerBoot2Regional(fit.nrm.null.sts,  nSim, x, "Date_2", 1))
nrm.null.smps <- sapply(deltas, function(x) powerBoot2Regional(fit.nrm.null.smps, nSim, x, "Date_2", 1))
nrm.pst       <- sapply(deltas, function(x) powerBoot2Regional(fit.nrm.pst,       nSim, x, "Date_2", 1))
# -------------------------------------
# Write out the power sets
# -------------------------------------
colnames(nrm.null)      <- deltas
write.csv(nrm.null, paste0(pwr.pth, "/",      analyte, "_site_vs_samples_null.csv"))
colnames(nrm.null.sts)  <- deltas
write.csv(nrm.null.sts, paste0(pwr.pth, "/",  analyte, "_site_vs_samples_null_plus_sites.csv"))
colnames(nrm.null.smps) <- deltas
write.csv(nrm.null.smps, paste0(pwr.pth, "/", analyte, "_site_vs_samples_nul_plus_samples.csv"))
colnames(nrm.pst)       <- deltas
write.csv(nrm.pst, paste0(pwr.pth, "/",       analyte, "_site_vs_samples_all.csv"))
# ---------------------
# Plot the power curves
# ---------------------
nrm.site.samp  <- rbind(data.frame(deltas,  scenario = "Null",    power = colMeans(nrm.null < 0.05)),
                        data.frame(deltas,  scenario = "Sites",   power = colMeans(nrm.null.sts < 0.05)),
                        data.frame(deltas,  scenario = "Samples", power = colMeans(nrm.null.smps < 0.05)),
                        data.frame(deltas,  scenario = "All",     power = colMeans(nrm.pst < 0.05)))
p.nrm.res.seas <-  ggplot(nrm.site.samp,
                          aes_string(x = "deltas",
                                     y = "power",
                                     group = "scenario",
                                     col = "scenario")) +
                          geom_line(lwd = 2) +
                          ylab("Power") +
                          xlab("Faction year-on-year change in linear trend")
png(filename =  paste0(fig.pth, "/", analyte, "_sites_vs_samples.png"),
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
print(p.nrm.res.seas)
dev.off()
# ------------------------------------------
# Visualise difference to null
# ------------------------------------------
# null.all <- data.frame(deltas   = deltas,
#                        scenario = "All (n=227) vs NULL (n=42)",
#                        power_diff = filter(nrm.site.samp, scenario == "All")$power/filter(nrm.site.samp, scenario == "Null")$power)
# null.samples <- data.frame(deltas   = deltas,
#                            scenario = "Samples (n=108) vs NULL",
#                            power_diff = filter(nrm.site.samp, scenario == "Samples")$power/filter(nrm.site.samp, scenario == "Null")$power)
# null.sites <- data.frame(deltas   = deltas,
#                          scenario = "Sites (n=83) vs NULL",
#                          power_diff = filter(nrm.site.samp, scenario == "Sites")$power/filter(nrm.site.samp, scenario == "Null")$power)
# pwr.ratio <- rbind(null.all, null.samples, null.sites)
# p.nrm.pwr.ratio <-  ggplot(pwr.ratio,
#                            aes_string(x = "deltas",
#                                       y = "power_diff",
#                                       group = "scenario",
#                                       col = "scenario")) +
#   geom_line(lwd = 2) +
#   ylab("Power ratio to null") +
#   xlab("Faction year-on-year change in linear trend")
# # Write out the plot
# # ------------------
# png(filename =  paste0(fig.pth, "/", analyte, "_power_ratios.png"),
#     bg = "white", width =  10, height = 5, units = 'in', res = 300)
# print(p.nrm.pwr.ratio)
# dev.off()

