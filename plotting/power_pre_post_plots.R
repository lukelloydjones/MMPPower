# =========================================================
# Compute the power at 0.1 fractional change for pre v post
# A measure for plotting over many regions and constituents
# and scenarios
# Author: Luke Lloyd-Jones
# Date started: 09/11/2020
# Date updated: 01/12/2020
# =========================================================
library(zoo)
library(reshape)
library(ggplot2)
library(dplyr)
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
for (region in regions)
{
  results     <- data.frame(matrix(0, nrow = 6, ncol = 5))
  results.avg <- data.frame(matrix(0, nrow = 6, ncol = 5))
  results.avg.pw.pt1 <- data.frame(matrix(0, nrow = 6, ncol = 5))
  rownames(results) <- analytes
  i <- 0
  for (analyte in analytes)
  {
    print(c(region, analyte))
    i <- i + 1
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
    # AUC
    AUC.pre    <- sum(diff(deltas) * rollmean(pre.2015.pwr, 2))
    AUC.pre.s1 <- sum(diff(deltas) * rollmean(pre.2015.pwr.s1, 2))
    AUC.pre.s2 <- sum(diff(deltas) * rollmean(pre.2015.pwr.s2, 2))
    AUC.pre.s3 <- sum(diff(deltas) * rollmean(pre.2015.pwr.s3, 2))
    AUC.pst    <- sum(diff(deltas) * rollmean(pst.2015.pwr, 2))
    # AVG power
    avg.pre    <- mean(pre.2015.pwr)
    avg.pre.s1 <- mean(pre.2015.pwr.s1)
    avg.pre.s2 <- mean(pre.2015.pwr.s2)
    avg.pre.s3 <- mean(pre.2015.pwr.s3)
    avg.pst    <- mean(pst.2015.pwr)
    # Power avg at plus minus 0.1
    avg.pre.pt1     <- mean(c(pre.2015.pwr["X.0.1"], pre.2015.pwr["X0.1"]))
    avg.pre.s1.pt1  <- mean(c(pre.2015.pwr.s1["X.0.1"], pre.2015.pwr.s1["X0.1"]))
    avg.pre.s2.pt1  <- mean(c(pre.2015.pwr.s2["X.0.1"], pre.2015.pwr.s2["X0.1"]))
    avg.pre.s3.pt1  <- mean(c(pre.2015.pwr.s3["X.0.1"], pre.2015.pwr.s3["X0.1"]))
    avg.pst.pt1     <- mean(c(pst.2015.pwr["X.0.1"], pst.2015.pwr["X0.1"]))
    # Load up the results
    results[i, ] <- c(AUC.pre, AUC.pre.s1, AUC.pre.s2, AUC.pre.s3, AUC.pst)
    results.avg[i, ] <- c(avg.pre, avg.pre.s1, avg.pre.s2, avg.pre.s3, avg.pst)
    results.avg.pw.pt1[i, ] <- c(avg.pre.pt1, avg.pre.s1.pt1, avg.pre.s2.pt1, avg.pre.s3.pt1, avg.pst.pt1)
  }
  if (region == "burdekin")
  {
    results.all        <- data.frame(results, NRM = region)
    results.avg.all    <- data.frame(results.avg, NRM = region)
    results.pw.pt1.all <- data.frame(results.avg.pw.pt1, NRM = region)
  } else {
    results.all        <- rbind(results.all,        data.frame(results, NRM = region))
    results.avg.all    <- rbind(results.avg.all,    data.frame(results.avg, NRM = region))
    results.pw.pt1.all <- rbind(results.pw.pt1.all, data.frame(results.avg.pw.pt1, NRM = region))
  }
}
# --------------------------
# Tidy up the results names
# --------------------------
results.pw.pt1.all$Analytes  <- rep(c("Chl-a", "TSS", "Secchi", "PP", "PN", "NOx"), times = 4)
colnames(results.pw.pt1.all) <- c("Pre-2015", "Pre-2015-S1", "Pre-2015-S2", "Pre-2015-S3", "Post-2015", "NRM", "Analytes")
results.pw.pt1.all$NRM       <- rep(c("Burdekin", "Wet Tropics - Russell-Mulgrave", "Mackay Whitsunday", "Wet Tropics - Tully"), each = 6)
results.pw.pt1.all <- results.pw.pt1.all[, -which(colnames(results.pw.pt1.all)  == "Pre-2015")]
# ------------------------------------------
# Plot a power Barplot with all time periods
# Supp Figure 6 
# ------------------------------------------
res.m <- melt(results.pw.pt1.all, id.vars = c("Analytes", "NRM"))
# Write out the plot
png(filename = paste0("results_hpc_loglm/global_figures/all_regions_global_power_at_pt1_pre_vs_post.png"), 
    bg = "white", width =  13, height = 13, units = 'in', res = 300)
ggplot(res.m, aes(x = as.factor(variable), y = value, fill = as.factor(variable))) +  geom_hline(aes(yintercept = 0.8), col = "grey65") +
       geom_bar(stat = "identity") + theme_web_bw() + facet_grid(NRM~Analytes) +
       ylab("Power at 0.1 fractional change") +
       xlab("") +
       theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank(), 
             legend.position = c(0.24, 0.16),
             axis.text.y = element_text(size = 15, face = "bold"),
             axis.title.y = element_text(size = 20, face = "bold"),
             text  = element_text(size = 25, face = "bold"),
             legend.key.size = unit("5", "mm")) + scale_fill_brewer(name = "Scenario", palette = "Set1") 
dev.off()

# ------------------------------------------------------
# The post with the average over the three time periods
# Figure 3 in main text
# ------------------------------------------------------
res.m <- melt(results.pw.pt1.all, id.vars  = c("Analytes", "NRM"))
#res.m <- res.m[-which(res.m$variable == "Pre-2015"), ]
res.m$NV <- rep("Pre-2015", dim(res.m)[1])
res.m$NV[which(res.m$variable == "Post-2015")] <- "Post-2015"
res.m <- res.m[, -3]
res.m.2 <- res.m %>% group_by(Analytes, NRM, NV) %>% summarise(mn = mean(value))
res.m.2$NV <- factor(res.m.2$NV , levels = c("Pre-2015", "Post-2015"))
# Write out the plot
png(filename = paste0("results_hpc_loglm/global_figures/all_regions_global_power_at_pt1_pre_vs_post_two_way.png"), 
    bg = "white", width =  12, height = 14, units = 'in', res = 300)
ggplot(res.m.2, aes(x = as.factor(NV), y = mn, fill = as.factor(NV))) +  geom_hline(aes(yintercept = 0.8), col = "grey65") +
  geom_bar(stat = "identity") + theme_web_bw() + facet_grid(NRM~Analytes) +
  ylab("Power at 0.1 fractional change") +
  xlab("") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = c(0.24, 0.65),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        text  = element_text(size = 25, face = "bold"),
        legend.key.size = unit("5", "mm")) + scale_fill_brewer(name = "Scenario", palette = "Set1") 
dev.off()
