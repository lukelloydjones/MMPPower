# =================================================
# Phase 3 investigation just for
#   - TSS, Secchi, Chl A, PN, PP, NOX
# Partition into pre, post, s1, s2, s3 to make it
# more safe and faster to run on cluster
# Author: Luke Lloyd-Jones
# Date started: 06/11/2020
# Date updated: 20/01/2021
# =================================================
args <- commandArgs(trailingOnly = TRUE)
# Example trailing arguments
args <- c("russ_mull", "PP_QAQC", "pre")
# Arguments are 1. nrm region # burdekin, russ_mull, tully, wtsndys
#               2. analyte - "CHL_QAQC" "SECCHI_DEPTH" "PP_QAQC" "PN_SHIM_QAQC" "SS_QAQC" "NOX_QAQC"
#               3. pre, post, s1, s2, s3
#               4. Maybe site but we could set this up as a separate thing later on
# Load some libraries
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
#library(MMP)
library(lubridate)
library(dplyr)
library(ggplot2)
library(car)
library(car)
library(mgcv)
library(ggwebthemes)
source("rscripts/power_boot_functions.R")
# A function to do the quantile residuals
# ----------------------
# File output path setup
# ----------------------
# Specify the region burdekin, fitzroy (not-analysed), russ_mull, tully, wtsndys
nrm.nm   <- as.character(args[1])
scenario <- as.character(args[3])
# CHOOSE A SITE TO FOCUS ON
site.sel     <- "RM8"
top.pth      <- "results"
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
nrm <- read.csv(paste0("data/raw/MMPDataset-sans-loggers-allfields_", nrm.nm, "_subset.csv"))
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
                                     , " PP ",       "Particulate phosphorus (micrograms per litre)"
                                     , " SAL ",      "Salinity"
                                     , " SI ",       "Dissolved silica"
                                     , " SS ",       "Total suspended solids-gravimetrically"
                                     , " TDN_PER ",  "Total dissolved nitrogen-persulphate digestion"
                                     , " TDN_SHIM ", "Total dissolved nitrogen-Shimadzu"
                                     , " TDP_PER ",  "Total dissolved phosphorus-persulphate digestion"
                                     , " TEMP ",     "Temperature"
                                     , " CHL ",      "Chlorophyll a (micrograms per litre)"
                                     , " PHAEO ",    "Phaeophytin a measured via filtration and fluorescence"
                                     , " SECCHI_DEPTH ",  "Secchi depth (metres)"), nrow = 2)))
# Collection dates
dates    <- as.Date(nrm$COLLECTION_START_DATE, "%d/%m/%y")
nrm$Date <- dates
# ==========================================
# Subset the data to the analyte of interest
# ==========================================
# Do a site or not ?
do.site <- 0
# Analyse priority constituents #"CHL_QAQC", "SECCHI_DEPTH"c("PP_QAQC", "PN_SHIM_QAQC", "SS_QAQC")
analyte <- as.character(args[2])
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
# ---------------------------------
# Subset to the analyte of interest
# ---------------------------------
if (analyte == "NOX_QAQC")
{
  # Subset for NOX and just keep BDLs as 0
  nrm.analyte <- nrm[, c("STATION_NAME", "DEPTH_CODE", "DUPLICATE",
                         "DEPTH", "NO2_QAQC", "NO3_QAQC", "LATITUDE", "LONGITUDE",
                         "ACOUSTIC_DEPTH", "SHORT_NAME", "PROJECT", "Date")]
  # Replace below detection limit with 0 then do something with the log normal
  # qqPlot for each analysis separately
  nrm.analyte$NOX_QAQC  <- nrm.analyte$NO2_QAQC + nrm.analyte$NO3_QAQC
  # Remove the JCU project
  if (nrm.nm != "wtsndys")
  {
    nrm.analyte <- nrm.analyte[-which(nrm.analyte$PROJECT == "MMP-JCU"), ]
  }
} else {
  nrm.analyte <- nrm[, c("STATION_NAME", "DEPTH_CODE", "DUPLICATE",
                         "DEPTH", analyte, "LATITUDE", "LONGITUDE", 
                         "ACOUSTIC_DEPTH", "SHORT_NAME", "PROJECT", "Date")]
}
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
# --------------------------
# OK this is the reduced set
# --------------------------
# Table depth by station name
summary(nrm.analyte.avg.red$DEPTH)
summary(nrm.analyte.avg.red$ACOUSTIC_DEPTH)
# -----------------------------------------------------------------
# Make a directory specific to the analyte to store all the figures
# -----------------------------------------------------------------
# Plot and save a depth versus measurement plot
if (scenario == "pre")
{
  plots <- 1
} else {
  plots <- 0
}
if (plots == 1)
{
  lm.dp <- lm(ANALYTE_MN ~ DEPTH, data = nrm.analyte.avg.red)
  main.val <- paste0("Depth coefs -- ", paste0("Effect = ", round(summary(lm.dp)$coefficients["DEPTH", 1], 3),
                                               ", P-val = ",  round(summary(lm.dp)$coefficients["DEPTH", 4], 3)))
  dpth.p <- ggplot(nrm.analyte.avg.red, aes(x = DEPTH, y = ANALYTE_MN)) + 
    geom_point() + ylab(analyte) + xlab("Depth (m)")  + ggtitle(main.val) +
    geom_smooth(method = lm, se = FALSE)
  png(filename = paste0(fig.pth, "/", analyte ,"_depth_versus_analyte.png"), 
      bg = "white", width =  8, height = 6, units = 'in', res = 300)
  print(dpth.p)
  dev.off()
}
# --------------------------------------------
# Average over the depths
#  - Something interesting to think about here
#    what depth averaging means 
# --------------------------------------------
nrm.analyte.avg.red.mn.dpth  <- as.data.frame(nrm.analyte.avg.red %>% 
                                                group_by(STATION_NAME, PROJECT, Date) %>% 
                                                mutate(ANALYTE_MN_DPTH = mean(ANALYTE_MN, na.rm = T)))
# Not really sure about the weighted mean?? Out of two measurement why should the deeper one have more relevance. 
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
# Remove DEPTH_CODE DEPTH TDN_MN
nrm.analyte.avg.red.mn.dpth.2 <- subset(nrm.analyte.avg.red.mn.dpth.2, 
                                        select = -c(DEPTH_CODE, DEPTH, ANALYTE_MN))
# --------------------------------------------------
# For NOX what does the distribution now look like ?
# --------------------------------------------------
# -------------------------------------------
# Produce a similar time series plots by site
#  This script gets run a lot so don't plot each time
# -------------------------------------------
if (plots == 1)
{
  p1 <-  ggplot(nrm.analyte.avg.red.mn.dpth.2, aes_string(x = "Date", 
                                                          y = "ANALYTE_MN_DPTH", 
                                                          group = "PROJECT", 
                                                          col   = "PROJECT")) +
    facet_wrap(.~SHORT_NAME, ncol = 2, scales = "free_y") + theme(legend.position = "top") + scale_colour_brewer(palette = "Set1") +
    theme(legend.position = "top") +
    geom_point() + theme_web_bw() +
    ylab(analyte.nms[grep(gsub("_QAQC", "", analyte), 
                          analyte.nms[, 1]), 2]) +
    xlab("Year") 
  png(filename = paste0(fig.pth, "/", analyte, "_averaged_rep_depth.png"), 
      bg = "white", width =  14, height = 8, units = 'in', res = 300)
  print(p1)
  dev.off()
}
# ------------------------------
# Process zeros if present
# ------------------------------
ana.vals <- nrm.analyte.avg.red.mn.dpth.2$ANALYTE_MN_DPTH
zeros <- which(ana.vals == 0)
if (length(which(ana.vals == 0)) > 0)
{
  nrm.analyte.avg.red.mn.dpth.2$ANALYTE_MN_DPTH[zeros] <- min(ana.vals[-zeros], na.rm = T) / 2
}
# ------------------------------
# Separate pre-2015 to post 2015
# ------------------------------
# Compute some common variables
nrm.analyte.avg.red.mn.dpth.2$PROJECT    <- as.factor(nrm.analyte.avg.red.mn.dpth.2$PROJECT)
nrm.analyte.avg.red.mn.dpth.2$SHORT_NAME <- as.factor(nrm.analyte.avg.red.mn.dpth.2$SHORT_NAME)
# Separate into pre and post
# --------------------------
nrm.analyte.pre.2015 <- filter(nrm.analyte.avg.red.mn.dpth.2, Date < as.Date("2015-01-01"))
nrm.analyte.pst.2015 <- filter(nrm.analyte.avg.red.mn.dpth.2, Date >= as.Date("2015-01-01"))
# One site in Whitsundays where they just took one more value post 2015. Remove
if (nrm.nm == "wtsndys")
{
  bad.nm <- names(which(table(nrm.analyte.pst.2015$SHORT_NAME) < 2))
  nrm.analyte.pst.2015 <- nrm.analyte.pst.2015[-which(nrm.analyte.pst.2015$SHORT_NAME == bad.nm), ]
}
# Post monitors for 1600 days so we take three snap shots of this length
# Some pre 2015 scenarios of about equal length in terms of time
pst.days <- (max(nrm.analyte.pst.2015$Date) - min(nrm.analyte.pst.2015$Date))
date.1   <- min(nrm.analyte.pre.2015$Date) + pst.days
date.2   <- max(nrm.analyte.pre.2015$Date) - pst.days
date.3   <- as.Date("2008-02-01") + pst.days
nrm.analyte.pre.2015.st1 <- filter(nrm.analyte.pre.2015, Date <= date.1)
nrm.analyte.pre.2015.st2 <- filter(nrm.analyte.pre.2015, Date >= date.2)
nrm.analyte.pre.2015.st3 <- filter(nrm.analyte.pre.2015, Date >= as.Date("2008-02-01") & 
                                     Date <= date.3)
max(nrm.analyte.pre.2015.st1$Date) - min(nrm.analyte.pre.2015.st1$Date)
max(nrm.analyte.pre.2015.st2$Date) - min(nrm.analyte.pre.2015.st2$Date)
max(nrm.analyte.pre.2015.st3$Date) - min(nrm.analyte.pre.2015.st3$Date)
# For each of the sets compute the minimum date and the seasonal terms
# relative to a Jan 1st date from the first year. This code is terrible 
# Pre all
T <- 365.25
nrm.analyte.pre.2015$Date_2  <- as.numeric(nrm.analyte.pre.2015$Date - min(nrm.analyte.pre.2015$Date))
nrm.analyte.pre.2015$Date_3  <- as.numeric(nrm.analyte.pre.2015$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pre.2015$Date)),"-01-01")))
nrm.analyte.pre.2015$SEAS_C1 <- cos(2 * pi * nrm.analyte.pre.2015$Date_3 / T)
nrm.analyte.pre.2015$SEAS_C2 <- sin(2 * pi * nrm.analyte.pre.2015$Date_3 / T)
# Pre S1
nrm.analyte.pre.2015.st1$Date_2  <- as.numeric(nrm.analyte.pre.2015.st1$Date - min(nrm.analyte.pre.2015.st1$Date))
nrm.analyte.pre.2015.st1$Date_3  <- as.numeric(nrm.analyte.pre.2015.st1$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pre.2015.st1$Date)),"-01-01")))
nrm.analyte.pre.2015.st1$SEAS_C1 <- cos(2 * pi * nrm.analyte.pre.2015.st1$Date_3 / T)
nrm.analyte.pre.2015.st1$SEAS_C2 <- sin(2 * pi * nrm.analyte.pre.2015.st1$Date_3 / T)
# Pre S2
nrm.analyte.pre.2015.st2$Date_2  <- as.numeric(nrm.analyte.pre.2015.st2$Date - min(nrm.analyte.pre.2015.st2$Date))
nrm.analyte.pre.2015.st2$Date_3  <- as.numeric(nrm.analyte.pre.2015.st2$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pre.2015.st2$Date)),"-01-01")))
nrm.analyte.pre.2015.st2$SEAS_C1 <- cos(2 * pi * nrm.analyte.pre.2015.st2$Date_3 / T)
nrm.analyte.pre.2015.st2$SEAS_C2 <- sin(2 * pi * nrm.analyte.pre.2015.st2$Date_3 / T)
# Pre S3
nrm.analyte.pre.2015.st3$Date_2  <- as.numeric(nrm.analyte.pre.2015.st3$Date - min(nrm.analyte.pre.2015.st3$Date))
nrm.analyte.pre.2015.st3$Date_3  <- as.numeric(nrm.analyte.pre.2015.st3$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pre.2015.st3$Date)),"-01-01")))
nrm.analyte.pre.2015.st3$SEAS_C1 <- cos(2 * pi * nrm.analyte.pre.2015.st3$Date_3 / T)
nrm.analyte.pre.2015.st3$SEAS_C2 <- sin(2 * pi * nrm.analyte.pre.2015.st3$Date_3 / T)
# Post
nrm.analyte.pst.2015$Date_2  <- as.numeric(nrm.analyte.pst.2015$Date - min(nrm.analyte.pst.2015$Date))
nrm.analyte.pst.2015$Date_3  <- as.numeric(nrm.analyte.pst.2015$Date) - as.numeric(as.Date(paste0(min(year(nrm.analyte.pst.2015$Date)),"-01-01")))
nrm.analyte.pst.2015$SEAS_C1 <- cos(2 * pi * nrm.analyte.pst.2015$Date_3 / T)
nrm.analyte.pst.2015$SEAS_C2 <- sin(2 * pi * nrm.analyte.pst.2015$Date_3 / T)
nrm.analyte.pst.2015$SEAS_C3 <- cos(4 * pi * nrm.analyte.pst.2015$Date_3 / T)
nrm.analyte.pst.2015$SEAS_C4 <- sin(4 * pi * nrm.analyte.pst.2015$Date_3 / T)
# --------------
# Now plot again
# --------------
if (plots == 1)
{
  p.pre.2015 <-  ggplot(nrm.analyte.pre.2015, 
                        aes_string(x = "Date", 
                                   y = "ANALYTE_MN_DPTH", 
                                   group = "PROJECT", 
                                   col = "PROJECT")) +
    geom_smooth(method="lm", se = F) + 
    facet_wrap(.~SHORT_NAME, ncol = 2) +
    geom_point() + ylab(analyte) + 
    xlab("Year")
  png(filename = paste0(fig.pth, "/", analyte, "_averaged_rep_depth_2005_2015.png"), 
      bg = "white", width =  14, height = 6, units = 'in', res = 300)
  print(p.pre.2015)
  dev.off()
  p.pst.2015 <-  ggplot(nrm.analyte.pst.2015, 
                        aes_string(x = "Date", 
                                   y = "ANALYTE_MN_DPTH", 
                                   group = "PROJECT", 
                                   col = "PROJECT")) +
    geom_smooth(method="lm", se = F) + 
    facet_wrap(.~SHORT_NAME, ncol = 2) +
    geom_point() + ylab(analyte) + 
    xlab("Year")
  png(filename = paste0(fig.pth, "/", analyte, "_averaged_rep_depth_2015_2019.png"), 
      bg = "white", width =  14, height = 6, units = 'in', res = 300)
  print(p.pst.2015)
  dev.off()
}
# ============================================================
# JUMP HERE TO nox_extra_investigation.R for power comparison
# by minimum adding 
# ============================================================
# ---------------------------------------
# Set the bootstrap simulation parameters
# ---------------------------------------
deltas   <- seq(-0.2, 0.2, 0.01) 
nSim     <- 1000
# ===========================================
# Site scale optional
# ===========================================
if (do.site == 1)
{
  # ============================================================
  # # Let's try for nrm pre and post
  # ============================================================
  nrm.analyte.pre.2015.nrm.st     <- filter(nrm.analyte.pre.2015,     SHORT_NAME == site.sel)
  nrm.analyte.pre.2015.nrm.st.st1 <- filter(nrm.analyte.pre.2015.st1, SHORT_NAME == site.sel)
  nrm.analyte.pre.2015.nrm.st.st2 <- filter(nrm.analyte.pre.2015.st2, SHORT_NAME == site.sel)
  nrm.analyte.pre.2015.nrm.st.st3 <- filter(nrm.analyte.pre.2015.st3, SHORT_NAME == site.sel)
  nrm.analyte.pst.2015.nrm.st     <- filter(nrm.analyte.pst.2015,     SHORT_NAME == site.sel)
  # --------------------------------
  # Run the linear models - for nrm
  # --------------------------------
  # ------------------
  # Run the models
  # ------------------
  fit.nrm.st.pre           <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2,  data = nrm.analyte.pre.2015.nrm.st)
  fit.nrm.st.pre.st1       <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2,  data = nrm.analyte.pre.2015.nrm.st.st1)
  fit.nrm.st.pre.st2       <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2,  data = nrm.analyte.pre.2015.nrm.st.st2)
  fit.nrm.st.pre.st3       <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2,  data = nrm.analyte.pre.2015.nrm.st.st3)
  fit.nrm.st.pst           <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2,  data = nrm.analyte.pst.2015.nrm.st)
  summary(fit.nrm.st.pre)
  summary(fit.nrm.st.pst)
  qqPlot(fit.nrm.st.pre$residuals)
  qqPlot(fit.nrm.st.pst$residuals)
  # -----
  # PRE
  # -----
  nrm.st.pre <- sapply(deltas, function(x) powerBoot2(fit.nrm.st.pre, nSim, x, "Date_2", "SEAS_C1", "SEAS_C2", 1, 1))
  colMeans(nrm.st.pre < 0.05)
  # Set 1
  nrm.st.pre.st1 <- sapply(deltas, function(x) powerBoot2(fit.nrm.st.pre.st1, nSim, x, "Date_2", "SEAS_C1", "SEAS_C2", 1, 1))
  colMeans(nrm.st.pre.st1 < 0.05)
  # Set 2
  nrm.st.pre.st2 <- sapply(deltas, function(x) powerBoot2(fit.nrm.st.pre.st2, nSim, x, "Date_2", "SEAS_C1", "SEAS_C2", 1, 1))
  colMeans(nrm.st.pre.st2 < 0.05)
  # Set 3
  nrm.st.pre.st3 <- sapply(deltas, function(x) powerBoot2(fit.nrm.st.pre.st3, nSim, x, "Date_2", "SEAS_C1", "SEAS_C2", 1, 1))
  colMeans(nrm.st.pre.st3 < 0.05)
  # -----
  # Post
  # -----
  nrm.st.pst <- sapply(deltas, function(x) powerBoot2(fit.nrm.st.pst, nSim, x, "Date_2", "SEAS_C1", "SEAS_C2", 1, 1))
  colMeans(nrm.st.pst < 0.05)
  # -------------------------------------
  # Write out the power sets
  # -------------------------------------
  colnames(nrm.st.pre)     <- deltas
  write.csv(nrm.st.pre, paste0(pwr.pth, "/nrm_", site.sel ,"_", analyte, "_power_pre_2015.csv"))
  colnames(nrm.st.pre.st1) <- deltas
  write.csv(nrm.st.pre.st1, paste0(pwr.pth, "/nrm_", site.sel ,"_", analyte, "_power_pre_2015_scenario_1.csv"))
  colnames(nrm.st.pre.st2) <- deltas
  write.csv(nrm.st.pre.st3, paste0(pwr.pth, "/nrm_", site.sel ,"_", analyte, "_power_pre_2015_scenario_2.csv"))
  colnames(nrm.st.pre.st3) <- deltas
  write.csv(nrm.st.pre.st3, paste0(pwr.pth, "/nrm_", site.sel ,"_", analyte, "_power_pre_2015_scenario_3.csv"))
  colnames(nrm.st.pst)     <- deltas
  write.csv(nrm.st.pst, paste0(pwr.pth, "/nrm_", site.sel ,"_", analyte, "_power_post_2015.csv"))
  # -------------------------------------
  # Make a data frame and plot the curves
  # -------------------------------------
  nrm.st.res <- rbind(data.frame(deltas, scenario = paste0(site.sel,"_Pre"), power = colMeans(nrm.st.pre < 0.05)),
                      data.frame(deltas, scenario = paste0(site.sel,"_S1"),  power = colMeans(nrm.st.pre.st1 < 0.05)),
                      data.frame(deltas, scenario = paste0(site.sel,"_S2"),  power = colMeans(nrm.st.pre.st2 < 0.05)),
                      data.frame(deltas, scenario = paste0(site.sel,"_S3"),  power = colMeans(nrm.st.pre.st3 < 0.05)),
                      data.frame(deltas, scenario = paste0(site.sel,"_Pst"), power = colMeans(nrm.st.pst < 0.05)))
  p.nrm.st.res <-  ggplot(nrm.st.res, 
                          aes_string(x = "deltas", 
                                     y = "power", 
                                     group = "scenario", 
                                     col = "scenario")) +
    geom_line(lwd = 2) + 
    ylab("Power") + 
    xlab("Fraction year-on-year change in linear trend")
  png(filename = paste0(fig.pth, "/", site.sel, "_", analyte, "_power_mod_res.png"), 
      bg = "white", width =  10, height = 5, units = 'in', res = 300)
  print(p.nrm.st.res)
  dev.off()
}
# ===========================================
# Region scale
# ===========================================
if (plots == 1)
{
  # Plot it first
  p.pre.2015 <-  ggplot(nrm.analyte.pre.2015, 
                        aes_string(x = "Date", 
                                   y = "ANALYTE_MN_DPTH", 
                                   group = "SHORT_NAME", 
                                   col = "SHORT_NAME")) +
    geom_smooth(method="lm", se = F) + 
    geom_point() + ylab(analyte) + 
    xlab("Year")
  # Not a whole lot of difference in intercepts with some difference in the slope
  # Ignoring differences in project now.
  png(filename = paste0(fig.pth, "/", analyte, "_averaged_rep_depth_2005_2015_region.png"), 
      bg = "white", width =  14, height = 6, units = 'in', res = 300)
  print(p.pre.2015)
  dev.off()
  p.pst.2015 <-  ggplot(nrm.analyte.pst.2015, 
                        aes_string(x = "Date", 
                                   y = "ANALYTE_MN_DPTH", 
                                   group = interaction(nrm.analyte.pst.2015$SHORT_NAME, nrm.analyte.pst.2015$PROJECT), 
                                   col = "SHORT_NAME")) +
    geom_smooth(method="lm", se = F) + 
    geom_point() + ylab(analyte) + 
    xlab("Year")
  png(filename = paste0(fig.pth, "/", analyte, "_averaged_rep_depth_2015_current_region.png"), 
      bg = "white", width =  14, height = 6, units = 'in', res = 300)
  print(p.pst.2015)
  dev.off()
}
# =====================
# Seasonal linear model
# =====================
# --------
# PRE 2015
# --------
# LOG LM
if (nrm.nm != "tully")
{
  fit.nrm.pre <- lm(log(ANALYTE_MN_DPTH) ~ SHORT_NAME + Date_2 + SEAS_C1 + SEAS_C2, 
                     data   = nrm.analyte.pre.2015)
} else {
  fit.nrm.pre <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2, 
                     data   = nrm.analyte.pre.2015)
}
summary(fit.nrm.pre)
png(filename = paste0(fig.pth, "/regional_", analyte, "_pre2015_qqplot_log_lm.png"), 
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
qqPlot(fit.nrm.pre$residuals)
dev.off()
# ----------
# Scenario 1
# ----------
if (nrm.nm != "tully")
{
  fit.nrm.pre.s1                <- lm(log(ANALYTE_MN_DPTH) ~ SHORT_NAME + Date_2 + SEAS_C1 + SEAS_C2, 
                                      data = nrm.analyte.pre.2015.st1)
} else {
  fit.nrm.pre.s1                <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2, 
                                      data = nrm.analyte.pre.2015.st1)
}
summary(fit.nrm.pre.s1)
png(filename = paste0(fig.pth, "/regional_", analyte, "_pre2015_s1_qqplot_log_lm.png"), 
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
qqPlot(fit.nrm.pre.s1$residuals)
dev.off()
# ----------
# Scenario 2
# ----------
if (nrm.nm != "tully")
{
  fit.nrm.pre.s2                <- lm(log(ANALYTE_MN_DPTH) ~ SHORT_NAME + Date_2 + SEAS_C1 + SEAS_C2, 
                                       data = nrm.analyte.pre.2015.st2)
} else {
  fit.nrm.pre.s2                <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2, 
                                       data = nrm.analyte.pre.2015.st2)
}
summary(fit.nrm.pre.s2)
png(filename = paste0(fig.pth, "/regional_", analyte, "_pre2015_s2_qqplot_log_lm.png"), 
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
qqPlot(fit.nrm.pre.s2$residuals)
dev.off()
# ----------
# Scenario 3
# ----------
if (nrm.nm != "tully")
{
  fit.nrm.pre.s3                <- lm(log(ANALYTE_MN_DPTH) ~ SHORT_NAME + Date_2 + SEAS_C1 + SEAS_C2, 
                                      data = nrm.analyte.pre.2015.st3)
} else {
  fit.nrm.pre.s3                <- lm(log(ANALYTE_MN_DPTH) ~ Date_2 + SEAS_C1 + SEAS_C2, 
                                      data = nrm.analyte.pre.2015.st3)
}
summary(fit.nrm.pre.s3)
png(filename = paste0(fig.pth, "/regional_", analyte, "_pre2015_s3_qqplot_log_lm.png"), 
    bg = "white", width =  10, height = 5, units = 'in', res = 300)
qqPlot(fit.nrm.pre.s3$residuals)
dev.off()
# ---------
# Post 2015
# ---------
if ( analyte == "NOX_QAQC" | nrm.nm == "wtsndys")
{
  # ONly one project for NOX - AIMS
  fit.nrm.pst          <- lm(log(ANALYTE_MN_DPTH) ~ SHORT_NAME + Date_2 + SEAS_C1 + SEAS_C2,
                              data = nrm.analyte.pst.2015)
  summary(fit.nrm.pst)
  png(filename = paste0(fig.pth, "/regional_", analyte, "_post2015_qqplot_log_lm.png"), bg = "white", width =  10, height = 5, units = 'in', res = 300)
  qqPlot(fit.nrm.pst$residuals)
  dev.off()
} else {
  fit.nrm.pst          <- lm(log(ANALYTE_MN_DPTH) ~ PROJECT + SHORT_NAME + Date_2 + SEAS_C1 + SEAS_C2,
                              data = nrm.analyte.pst.2015)
  summary(fit.nrm.pst)
  png(filename = paste0(fig.pth, "/regional_", analyte, "_post2015_qqplot_log_lm.png"), bg = "white", width =  10, height = 5, units = 'in', res = 300)
  qqPlot(fit.nrm.pst$residuals)
  dev.off()
} 
# Interaction models
if (analyte != "NOX_QAQC" & nrm.nm != "wtsndys")
{
  # Look at trend interaction by project or site 
  fit.nrm.pst.int.prj  <- lm(log(ANALYTE_MN_DPTH) ~  SHORT_NAME + Date_2 * PROJECT  + Date_2 + SEAS_C1 + SEAS_C2,
                              data = nrm.analyte.pst.2015)
  fit.nrm.pst.int.st   <- lm(log(ANALYTE_MN_DPTH) ~  Date_2 * SHORT_NAME + Date_2 * PROJECT  + Date_2 + SEAS_C1 + SEAS_C2,
                              data = nrm.analyte.pst.2015)
}
# ------------------------------------------
# Write out models and informative model output
# ------------------------------------------
if (scenario == "pre")
{
  save(fit.nrm.pre,    file = paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_pre_2015_model_fit_log_lm.Rdata"))
  sink(paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_pre_2015_model_fit_log_lm.txt"))
  print(summary(fit.nrm.pre))
  sink()
  
} else if (scenario == "s1") {
  
  save(fit.nrm.pre.s1, file = paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_pre_2015_model_fit_s1_log_lm.Rdata"))
  sink(paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_pre_2015_model_fit_s1_log_lm.txt"))
  print(summary(fit.nrm.pre.s1))
  sink()
  
} else if (scenario == "s2") {
  
  save(fit.nrm.pre.s2, file = paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_pre_2015_model_fit_s2_log_lm.Rdata"))
  sink(paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_pre_2015_model_fit_s2_log_lm.txt"))
  print(summary(fit.nrm.pre.s2))
  sink()
  
} else if (scenario == "s3") { 
  
  save(fit.nrm.pre.s3, file = paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_pre_2015_model_fit_s3_log_lm.Rdata"))
  sink(paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_pre_2015_model_fit_s3_log_lm.txt"))
  print(summary(fit.nrm.pre.s3))
  sink()
  
} else if (scenario == "post") {
  save(fit.nrm.pst,    file = paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_post_2015_model_fit_log_lm.Rdata"))
  sink(paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_post_2015_model_fit_log_lm.txt"))
  print(summary(fit.nrm.pst))
  sink()
  # Interaction models - non fitted for NOX
  # ------------------
  if ( analyte != "NOX_QAQC" & nrm.nm != "wtsndys")
  {
    
    save(fit.nrm.pst.int.prj, file = paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_post_2015_model_fit_int_prj_log_lm.Rdata"))
    save(fit.nrm.pst.int.st,  file = paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_post_2015_model_fit_int_site_log_lm.Rdata"))
    
    sink(paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_post_2015_model_fit_int_prj_log_lm.txt"))
    print(summary(fit.nrm.pst.int.prj))
    sink()
    
    sink(paste0(mdl.pth, "/", nrm.nm,"_", analyte, "_post_2015_model_fit_int_site_log_lm.txt"))
    print(summary(fit.nrm.pst.int.st))
    sink()
  }
}
# =========================
# Let's do the power 
# =========================
if (scenario == "pre")
{
  
  nrm.pre   <- sapply(deltas, function(x) powerBoot2Regional(fit.nrm.pre,    nSim, x, "Date_2", 1))
  colnames(nrm.pre)   <- deltas
  write.csv(nrm.pre, paste0(pwr.pth, "/", nrm.nm, "_",   analyte, "_power_pre_2015_log_lm.csv"))
  
} else if (scenario == "s1") {
  
  nrm.pre.1 <- sapply(deltas, function(x) powerBoot2Regional(fit.nrm.pre.s1, nSim, x, "Date_2", 1))
  colnames(nrm.pre.1) <- deltas
  write.csv(nrm.pre.1, paste0(pwr.pth, "/", nrm.nm, "_", analyte, "_power_pre_2015_scenario_1_log_lm.csv"))
  
} else if (scenario == "s2") {
  
  nrm.pre.2 <- sapply(deltas, function(x) powerBoot2Regional(fit.nrm.pre.s2, nSim, x, "Date_2", 1))
  colnames(nrm.pre.2) <- deltas
  write.csv(nrm.pre.2, paste0(pwr.pth, "/", nrm.nm, "_", analyte, "_power_pre_2015_scenario_2_log_lm.csv"))
  
} else if (scenario == "s3") {
  
  nrm.pre.3 <- sapply(deltas, function(x) powerBoot2Regional(fit.nrm.pre.s3, nSim, x, "Date_2", 1))
  colnames(nrm.pre.3) <- deltas
  write.csv(nrm.pre.3, paste0(pwr.pth, "/", nrm.nm, "_", analyte, "_power_pre_2015_scenario_3_log_lm.csv"))
  
} else if (scenario == "post") {
  
  nrm.pst   <- sapply(deltas, function(x) powerBoot2Regional2(fit.nrm.pst,    nSim, x, "Date_2", 1))
  colnames(nrm.pst)   <- deltas
  write.csv(nrm.pst, paste0(pwr.pth, "/", nrm.nm, "_",   analyte, "_power_post_2015_log_lm.csv"))
  
}





