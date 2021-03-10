# ======================================================
# Compute the area under the power curve as a single 
# measure for plotting over many regions and constituents
# and scenarios. Site vs sampling
# Author: Luke Lloyd-Jones
# Date started: 16/11/2020
# Date updated: 23/02/2020
# ======================================================
library(zoo)
library(reshape)
library(ggplot2)
#devtools::install_gitlab("peterbar/ggwebthemes")
library(ggwebthemes)
region   <- "wtsndys"
top.pth <- "github/plotting/results/sites_vs_samples/"
dir.create(file.path(top.pth, "global_figures"))
# ----------------------------------
# Read in the regional power results 
# ----------------------------------
# LEave PN out as it's a bit weird with JCU
regions   <- c("burdekin", "russ_mull", "wtsndys", "tully")
analytes <- c("CHL_QAQC", "SS_QAQC", "SECCHI_DEPTH", "PP_QAQC", "PN_SHIM_QAQC")
for (region in regions)
{
  results  <- data.frame(matrix(0, nrow = 5, ncol = 8))
  results.avg        <- data.frame(matrix(0, nrow = 5, ncol = 8))
  results.avg.pw.pt1 <- data.frame(matrix(0, nrow = 5, ncol = 8))
  rownames(results) <- analytes
  i <- 0
  for (analyte in analytes)
  {
    print(c(region, analyte))
    i <- i + 1
    data.in.null    <- paste0(top.pth, 
                              region, "_phase_3_analysis/power_sets/", analyte,
                              "/", analyte, "_", "site_vs_samples_null.csv")
    data.in.sites   <- paste0(top.pth, 
                              region, "_phase_3_analysis/power_sets/", analyte,
                              "/", analyte, "_", "site_vs_samples_null_plus_sites.csv")
    data.in.samples <- paste0(top.pth, 
                              region, "_phase_3_analysis/power_sets/", analyte,
                              "/",  analyte, "_", "site_vs_samples_nul_plus_samples.csv")
    data.in.all     <- paste0(top.pth, 
                              region, "_phase_3_analysis/power_sets/", analyte,
                              "/", analyte, "_", "site_vs_samples_all.csv")
    # Model in
    mod.in.null    <- paste0(top.pth, 
                              region, "_phase_3_analysis/model_fits/", analyte,
                              "/", analyte, "_", "site_vs_samples_null.Rdata")
    mod.in.sites   <- paste0(top.pth, 
                              region, "_phase_3_analysis/model_fits/", analyte,
                              "/", analyte, "_", "site_vs_samples_null_plus_sites.Rdata")
    mod.in.samples <- paste0(top.pth, 
                              region, "_phase_3_analysis/model_fits/", analyte,
                              "/",  analyte, "_", "site_vs_samples_nul_plus_samples.Rdata")
    mod.in.all     <- paste0(top.pth, 
                              region, "_phase_3_analysis/model_fits/", analyte,
                              "/", analyte, "_", "site_vs_samples_all.Rdata")
    # Read in the data
    null    <- read.csv(data.in.null)
    sites   <- read.csv(data.in.sites)
    samples <- read.csv(data.in.samples)
    all     <- read.csv(data.in.all)
    # Read in the models
    null.mod    <- get(load(mod.in.null))
    sites.mod   <- get(load(mod.in.sites))
    samples.mod <- get(load(mod.in.samples))
    all.mod     <- get(load(mod.in.all))
    # Pre versus post
    deltas       <- seq(-0.2, 0.2, 0.01) 
    # Compute the power and then the AUC
    null.pwr    <- colMeans(null[, -1]    < 0.05)
    sites.pwr   <- colMeans(sites[, -1]   < 0.05)
    samples.pwr <- colMeans(samples[, -1] < 0.05)
    all.pwr     <- colMeans(all[, -1]     < 0.05)
    # AUC
    AUC.null    <- sum(diff(deltas) * rollmean(null.pwr, 2))
    AUC.sites   <- sum(diff(deltas) * rollmean(sites.pwr, 2))
    AUC.samples <- sum(diff(deltas) * rollmean(samples.pwr, 2))
    AUC.all     <- sum(diff(deltas) * rollmean(all.pwr, 2))
    # AVG power
    avg.null    <- mean(null.pwr)
    avg.sites   <- mean(sites.pwr)
    avg.samples <- mean(samples.pwr)
    avg.all     <- mean(all.pwr)
    # Power avg at plus minus 0.1
    avg.null.pt1    <- mean(c(null.pwr["X.0.1"],    null.pwr["X0.1"]))
    avg.sites.pt1   <- mean(c(sites.pwr["X.0.1"],   sites.pwr["X0.1"]))
    avg.samples.pt1 <- mean(c(samples.pwr["X.0.1"], samples.pwr["X0.1"]))
    avg.all.pt1     <- mean(c(all.pwr["X.0.1"],     all.pwr["X0.1"]))
    # Load up the results
    mod.ns <- c(length(null.mod$residuals), length(sites.mod$residuals), length(samples.mod$residuals), length(all.mod$residuals))
    results[i, ]            <- c(AUC.null, AUC.sites, AUC.samples, AUC.all, mod.ns)
    results.avg[i, ]        <- c(avg.null, avg.sites, avg.samples, avg.all, mod.ns)
    results.avg.pw.pt1[i, ] <- c(avg.null.pt1, avg.sites.pt1, avg.samples.pt1, avg.all.pt1, mod.ns)
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


# ----------------------
# Plot a power Barplot
# Figure 4 in main text
# ---------------------

results.pw.pt1.all$Analytes  <- rep(c("Chl-a", "TSS", "Secchi", "PP", "PN"), times = 4)
colnames(results.pw.pt1.all) <- c("Null", "Null + sites", "Null + samples",  "Current sampling", "N-Null", "N-Null + sites", "N-Null + samples",  "N-Current sampling", "NRM", "Analytes")
results.pw.pt1.all$NRM       <- rep(c("Burdekin", "Wet Tropics - Russell-Mulgrave", "Mackay Whitsunday", "Wet Tropics - Tully"), each = 5)


results.pw.pt1.all.sub <- results.pw.pt1.all[, c("Null", "Null + sites", "Null + samples",  "Current sampling", "Analytes", "NRM")]
res.m <- melt(results.pw.pt1.all.sub, id.vars = c("Analytes", "NRM"))

results.avg.all.sub.n <- results.pw.pt1.all[, c("Analytes", "NRM", "N-Null", "N-Null + sites", "N-Null + samples",  "N-Current sampling")]
colnames(results.avg.all.sub.n) <- c( "Analytes", "NRM", "Null", "Null + sites", "Null + samples",  "Current sampling")
res.m.n              <- melt(results.avg.all.sub.n, id.vars = c("Analytes", "NRM"))
res.m.n$pwr <- res.m$value

# Write out the plot
png(filename = paste0(top.path, "/global_figures/all_regions_global_power_at_pt1_sites_vs_samples.png"), 
    bg = "white", width =  13, height = 13, units = 'in', res = 300)
ggplot(res.m, aes(x = as.factor(variable), y = value, fill = as.factor(variable))) + geom_hline(aes(yintercept = 0.8), col = "grey65") +
  geom_bar(stat = "identity") + theme_web_bw() + facet_grid(NRM~Analytes) +
  ylab("Power at 0.1 fractional change") +
  xlab("") +
  geom_text(data = res.m.n, aes(label = value, y = pwr + 0.05), size = 5) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = c(0.91, 0.68),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        text  = element_text(size = 25, face = "bold"),
        legend.key.size = unit("5", "mm")) + scale_fill_brewer(name = "Scenario", palette = "Set1") 
dev.off()

