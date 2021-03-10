# ======================================================
# Try and draw crayon plots for all the analytes and 
# regions analysed
# Author: Luke Lloyd-Jones
# Date started: 07/12/2020
# Date updated: 23/02/2020
# ======================================================
library(zoo)
library(reshape)
library(ggplot2)
library(dplyr)
library(ggwebthemes)
#devtools::install_gitlab("peterbar/ggwebthemes")
library(ggwebthemes)
system("mkdir github/plotting/results/results_hpc_loglm/global_figures/")
# Some path orientation
top.path <- "github/plotting/results/results_hpc_loglm/"
mdl      <- "_log_lm"
# ----------------------------------
# Read in the regional power results 
# ----------------------------------
regions  <- c("burdekin", "russ_mull", "wtsndys", "tully")
analytes <- c("CHL_QAQC", "SS_QAQC", "SECCHI_DEPTH", "PP_QAQC", "PN_SHIM_QAQC", "NOX_QAQC")
deltas   <- seq(-0.2, 0.2, 0.01)
for (region in regions)
{
  for (analyte in analytes)
  {
    print(c(region, analyte))
    data.in.pre <- paste0(top.path, 
                          region, "_phase_3_analysis/power_sets/", analyte,
                          "/", region,"_", analyte, "_", "power_pre_2015", mdl, ".csv")
    data.in.pre.s1 <- paste0(top.path, 
                             region, "_phase_3_analysis/power_sets/", analyte,
                             "/", region,"_", analyte, "_", "power_pre_2015_scenario_1", mdl, ".csv")
    data.in.pre.s2 <- paste0(top.path, 
                             region, "_phase_3_analysis/power_sets/", analyte,
                             "/", region,"_", analyte, "_", "power_pre_2015_scenario_2", mdl, ".csv")
    data.in.pre.s3 <- paste0(top.path, 
                             region, "_phase_3_analysis/power_sets/", analyte,
                             "/", region,"_", analyte, "_", "power_pre_2015_scenario_3", mdl, ".csv")
    data.in.pst <- paste0(top.path, 
                          region, "_phase_3_analysis/power_sets/", analyte,
                          "/", region,"_", analyte, "_", "power_post_2015", mdl, ".csv")
    # Read in the data
    pre.2015     <- read.csv(data.in.pre)
    pre.2015.s1  <- read.csv(data.in.pre.s1)
    pre.2015.s2  <- read.csv(data.in.pre.s2)
    pre.2015.s3  <- read.csv(data.in.pre.s3)
    pst.2015     <- read.csv(data.in.pst)
    # Pre versus post
    deltas       <- seq(-0.2, 0.2, 0.01) 
    # Compute the power and then the AUC
    pre.2015.pwr    <- colMeans(pre.2015[, -1]    < 0.05)
    pre.2015.pwr.s1 <- colMeans(pre.2015.s1[, -1] < 0.05)
    pre.2015.pwr.s2 <- colMeans(pre.2015.s2[, -1] < 0.05)
    pre.2015.pwr.s3 <- colMeans(pre.2015.s3[, -1] < 0.05)
    pst.2015.pwr    <- colMeans(pst.2015[, -1]    < 0.05)
    # Load up the results
    if (analyte == "CHL_QAQC")
    {
      results <- data.frame(analyte = "CHL_QAQC", delta = deltas, power = pre.2015.pwr, scenario = "Pre-2015")
      results <- rbind(results, data.frame(analyte = "CHL_QAQC", delta = deltas, power = pre.2015.pwr.s1, scenario = "Pre-2015-S1"))
      results <- rbind(results, data.frame(analyte = "CHL_QAQC", delta = deltas, power = pre.2015.pwr.s2, scenario = "Pre-2015-S2"))
      results <- rbind(results, data.frame(analyte = "CHL_QAQC", delta = deltas, power = pre.2015.pwr.s3, scenario = "Pre-2015-S3"))
      results <- rbind(results, data.frame(analyte = "CHL_QAQC", delta = deltas, power = pst.2015.pwr, scenario = "Post-2015"))
    } else {
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = pre.2015.pwr, scenario = "Pre-2015"))
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = pre.2015.pwr.s1, scenario = "Pre-2015-S1"))
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = pre.2015.pwr.s2, scenario = "Pre-2015-S2"))
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = pre.2015.pwr.s3, scenario = "Pre-2015-S3"))
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = pst.2015.pwr, scenario = "Post-2015"))
    }
  }
  if (region == "burdekin")
  {
    results.all        <- data.frame(results, NRM = "Burdekin")
  } else {
    results.all        <- rbind(results.all,  data.frame(results, NRM = region))
  }
  rm(results)
}
# ---------------------
# Plot the power curves
# ---------------------
rep(c("Chl-a", "TSS", "Secchi", "PP", "PN", "NOx"), times = 4)
results.all$analyte  <- gsub("CHL_QAQC", "Chl-a", results.all$analyte)
results.all$analyte  <- gsub("SS_QAQC",  "TSS",   results.all$analyte)
results.all$analyte  <- gsub("SECCHI_DEPTH", "Secchi", results.all$analyte)
results.all$analyte  <- gsub("PN_SHIM_QAQC", "PN", results.all$analyte)
results.all$analyte  <- gsub("NOX_QAQC", "NOx", results.all$analyte)
results.all$analyte  <- gsub("PP_QAQC", "PP", results.all$analyte)
# Replace NRMs
results.all$NRM  <- gsub("russ_mull", "Wet Tropics - Russell-Mulgrave", results.all$NRM)
results.all$NRM  <- gsub("wtsndys", "Mackay Whitsunday", results.all$NRM)
results.all$NRM  <- gsub("tully", "Wet Tropics - Tully", results.all$NRM)
colnames(results.all)[4] <- "Scenario"
results.all <- results.all[-which(results.all$Scenario == "Pre-2015"), ]
results.all$delta <- round(results.all$delta, 2)
# Write out the plot
png(filename = paste0("results_hpc_loglm/global_figures/crayon_all_pre_post.png"), 
    bg = "white", width =  13, height = 13, units = 'in', res = 300)
ggplot(results.all, aes(x = delta, y = power, group = Scenario, col = Scenario)) +
       geom_line(stat = "identity", lwd = 0.8) + theme_web_bw() + facet_grid(NRM~analyte) +  
       geom_hline(aes(yintercept = 0.8), col = "grey65") + geom_hline(aes(yintercept = 0.05), col = "grey65") +
       ylab("Power") +
       xlab("Fractional year-on-year linear change") +
       theme(axis.text.y = element_text(size = 15, face = "bold"),
             axis.title.y = element_text(size = 20, face = "bold"),
             text  = element_text(size = 25, face = "bold"),
       legend.key.size = unit("5", "mm")) + scale_colour_brewer(name = "Scenario", palette = "Set1") 
dev.off()


# =======================================================
# Sites vs samples
# =======================================================
top.path <- "github/plotting/results/sites_vs_samples/"
dir.create(file.path(top.path, "global_figures"))
# ----------------------------------
# Read in the regional power results 
# ----------------------------------
# LEave PN out as it's a bit weird with JCU
regions   <- c("burdekin", "russ_mull", "wtsndys", "tully")
analytes <- c("CHL_QAQC", "SS_QAQC", "SECCHI_DEPTH", "PP_QAQC", "PN_SHIM_QAQC")
deltas   <- seq(-0.2, 0.2, 0.01)
for (region in regions)
{
  for (analyte in analytes)
  {
    data.in.null    <- paste0(top.path, 
                              region, "_phase_3_analysis/power_sets/", analyte,
                              "/", analyte, "_", "site_vs_samples_null.csv")
    data.in.sites   <- paste0(top.path, 
                              region, "_phase_3_analysis/power_sets/", analyte,
                              "/", analyte, "_", "site_vs_samples_null_plus_sites.csv")
    data.in.samples <- paste0(top.path, 
                              region, "_phase_3_analysis/power_sets/", analyte,
                              "/",  analyte, "_", "site_vs_samples_nul_plus_samples.csv")
    data.in.all     <- paste0(top.path, 
                              region, "_phase_3_analysis/power_sets/", analyte,
                              "/", analyte, "_", "site_vs_samples_all.csv")
    # Read in the data
    null    <- read.csv(data.in.null)
    sites   <- read.csv(data.in.sites)
    samples <- read.csv(data.in.samples)
    all     <- read.csv(data.in.all)
    # Pre versus post
    deltas       <- seq(-0.2, 0.2, 0.01) 
    # Compute the power and then the AUC
    null.pwr    <- colMeans(null[, -1]    < 0.05)
    sites.pwr   <- colMeans(sites[, -1]   < 0.05)
    samples.pwr <- colMeans(samples[, -1] < 0.05)
    all.pwr     <- colMeans(all[, -1]     < 0.05)
    # Load up the results
    if (analyte == "CHL_QAQC")
    {
      results <- data.frame(analyte = "CHL_QAQC", delta = deltas, power = null.pwr, scenario = "Null")
      results <- rbind(results, data.frame(analyte = "CHL_QAQC", delta = deltas, power = sites.pwr, scenario = "Null + sites"))
      results <- rbind(results, data.frame(analyte = "CHL_QAQC", delta = deltas, power = samples.pwr, scenario = "Null + samples"))
      results <- rbind(results, data.frame(analyte = "CHL_QAQC", delta = deltas, power = all.pwr, scenario = "Current sampling"))
    } else {
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = null.pwr, scenario = "Null"))
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = sites.pwr, scenario = "Null + sites"))
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = samples.pwr, scenario = "Null + samples"))
      results <- rbind(results, data.frame(analyte = analyte, delta = deltas, power = all.pwr, scenario = "Current sampling"))
    }
  }
  if (region == "burdekin")
  {
    results.all        <- data.frame(results, NRM = "Burdekin")
  } else {
    results.all        <- rbind(results.all,  data.frame(results, NRM = region))
  }
  rm(results)
}

# ----------
# Do a plot
# ----------


results.all$analyte  <- gsub("CHL_QAQC", "Chl-a", results.all$analyte)
results.all$analyte  <- gsub("SS_QAQC",  "TSS",   results.all$analyte)
results.all$analyte  <- gsub("SECCHI_DEPTH", "Secchi", results.all$analyte)
results.all$analyte  <- gsub("PN_SHIM_QAQC", "PN", results.all$analyte)
results.all$analyte  <- gsub("NOX_QAQC", "NOx", results.all$analyte)
results.all$analyte  <- gsub("PP_QAQC", "PP", results.all$analyte)
# Replace NRMs
results.all$NRM  <- gsub("russ_mull", "Wet Tropics - Russell-Mulgrave", results.all$NRM)
results.all$NRM  <- gsub("wtsndys", "Mackay Whitsunday", results.all$NRM)
results.all$NRM  <- gsub("tully", "Wet Tropics - Tully", results.all$NRM)
colnames(results.all)[4] <- "Scenario"
# Write out the plot
png(filename = paste0(top.path, "global_figures/crayon_all_sites_samples.png"), 
    bg = "white", width =  13, height = 13, units = 'in', res = 300)
ggplot(results.all, aes(x = delta, y = power, group = Scenario, col = Scenario)) + geom_hline(aes(yintercept = 0.8), col = "grey65") + geom_hline(aes(yintercept = 0.05), col = "grey65") +
  geom_line(stat = "identity", lwd = 0.8) + theme_web_bw() + facet_grid(NRM~analyte)  +
  ylab("Power") +
  xlab("Fractional year-on-year linear change") +
  theme(axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        text  = element_text(size = 25, face = "bold"),
        legend.key.size = unit("5", "mm")) + scale_colour_brewer(name = "Scenario", palette = "Set1") 
dev.off()
