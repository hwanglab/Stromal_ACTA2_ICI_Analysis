library("survminer");
library("survival")

#-------------------------------------------------------------------------------
# Figure 2. ACTA2 expression was independently associated with overall survival 
# in gastric cancer patients. 

# The following codes reproduce the Figure 2-(B, C and D) which were generated using the pooled data 

# The pooled data was combined from cohorts previously published by The Cancer Genome Atlas Project (TCGA),
# Asian Cancer Research Group (ACRG), Sohn et al,[28] Kim et al, For the pooled data (n=974), 
# we used "ComBat" in "sva" package (R version 4.1.1) to remove possible batch effects in the expression values across the data sets. 

# TCGA: Both gene expression values and the clinical variables
#       were downloaded from CBioPortal (https://www.cbioportal.org/)
# ACRG: Cristescu, R., et al., Molecular analysis of gastric cancer identifies 
#       subtypes associated with distinct clinical outcomes. Nature medicine, 2015. 21(5): p. 449-456. 
# Sohn et al.: Sohn, B.H., et al., Clinical Significance of Four Molecular 
#              Subtypes of Gastric Cancer Identified by The Cancer Genome Atlas Project. 
#              Clinical Cancer Research, 2017. 23(15): p. 4441-4449 

# ACTA2 expression values are obtained after the pooled data was normalized by ComBat
#-------------------------------------------------------------------------------

# the current working directory 
# setwd("D:\\Research\\codes\\00_distributions\\Stromal_ACTA2_ICI_Analysis\\codes")

#- load data
data_f <- read.table("../data/pooled_data_ACTA2_MSI_EBV.csv", sep=",", header=TRUE);
data_f$Stage = factor(data_f$Stage)

#-------------------------------------------------------------------------------
#- Figure 2-B: Multivariable analysis of the pooled cohort. 
fit_mv <- coxph(formula = Surv(surv_time_m, vital_sign_1_d)
                ~ Age + Sex_ + Stage + ACTA2, data = data_f)

# result summary
summary(fit_mv)

beta <- exp(coef(fit_mv))
CI   <- round(exp(confint(fit_mv)), 5)

coeffs <- coef(summary(fit_mv))
p    <- as.matrix(coeffs[,5])

dir.create("../results/", showWarnings = FALSE)

fileConn<-file("../results/Figure2B_additional.txt", open = 'w')
for (mn_i in 1:length(beta)){
  str_line = sprintf("%10s:: %.2f (%.2f, %.2f), pval = %.5f", 
                     rownames(coeffs)[mn_i], beta[mn_i], CI[mn_i,1], CI[mn_i,2], p[mn_i])
  writeLines(str_line, fileConn)
}
close(fileConn)

# Visualize the results from the Cox model using forest plot
pdf_fn = file.path("../results/Figure2B.pdf")
message(pdf_fn)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file=pdf_fn,width=5,height=5)
ggforest(fit_mv, data = data_f)
dev.off()
message("Done.")



#-------------------------------------------------------------------------------
#- Figure 2 - (C & D): Kaplan-Meier plots
# data_f$MSI[301:567] = NA

data_f_sub = data_f[!(is.na(data_f$MSI) | is.na(data_f$surv_time_m)),]

#- Microsatellite instability-high (MSI-H) patients 
data_f_MSI = data_f_sub[data_f_sub$MSI==1, ]

# cutoff for MSI only samples
mr_thres = quantile(data_f_MSI$ACTA2, c(0.25))
print(paste0('Cutoff (MSI-H)  ', sprintf('%.4f', mr_thres)))

x = rep(NA, rep=nrow(data_f_MSI));
x[data_f_MSI$ACTA2 <= mr_thres] = 1
x[data_f_MSI$ACTA2  > mr_thres] = 2

data_f_MSI = cbind(data_f_MSI, ACTA2_bin = factor(x))

surv_MSI = data.frame(data_f_MSI$surv_time_m,  data_f_MSI$vital_sign_1_d, data_f_MSI$ACTA2_bin)
colnames(surv_MSI) <- c("time", "status", "gInfo")


#- Microsatellite stable (MSS) patients
data_f_MSS = data_f_sub[data_f_sub$MSI==0, ]

# cutoff for MSS only samples
mr_thres = quantile(data_f_MSS$ACTA2, c(0.25))
print(paste0('Cutoff (MSS)  ', sprintf('%.4f', mr_thres)))

x = rep(NA, rep=nrow(data_f_MSI));
x[data_f_MSS$ACTA2 <= mr_thres] = 1
x[data_f_MSS$ACTA2  > mr_thres] = 2

data_f_MSS = cbind(data_f_MSS, ACTA2_bin = factor(x))

surv_MSS = data.frame(data_f_MSS$surv_time_m,  data_f_MSS$vital_sign_1_d, data_f_MSS$ACTA2_bin)
colnames(surv_MSS) <- c("time", "status", "gInfo")
write.table(surv_MSS, file = "../results/data_Figure2D_ACTA2_low_high_in_MSS.txt", sep = "\t")

# Figure 2 - C: Kaplan-Meier survival analysis of microsatellite instability-high patients 
# stratified by ACTA2 expression 
fit_MSI <- survfit(Surv(time, status) ~ gInfo, surv_MSI)
p <- ggsurvplot(
  fit_MSI,                     # survfit object with calculated statistics.
  data = surv_MSI,
  surv.median.line = "hv", # Add medians survival
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("black", "red"),
  xlim = c(0, 120),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Months",   # customize X axis label.
  censor.shape = "|",
  censor.size = 2,
  legend.labs = c("low", "high")
) 

pdf_fn = file.path("../results/Figure2C.pdf")
message(pdf_fn)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file=pdf_fn,width=5,height=5)
p
dev.off()
message("Done.")

print(sprintf('[Figure 2-C]  ACTA2-Low: %3d vs ACTA2-High: %3d', sum(data_f_MSI$ACTA2_bin==1), sum(data_f_MSI$ACTA2_bin==2)))
print(fit_MSI) # print some information including the median survival times

# Figure 2 - D: Kaplan-Meier survival analysis of microsatellite stable patients 
# stratified by ACTA2 expression
fit_MSS <- survfit(Surv(time, status) ~ gInfo, surv_MSS)
p <- ggsurvplot(
  fit_MSS,                     # survfit object with calculated statistics.
  data = surv_MSS,             # data used to fit survival curves. 
  surv.median.line = "hv", # Add medians survival
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("black", "red"),
  xlim = c(0, 120),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Months",   # customize X axis label.
  censor.shape = "|",
  censor.size = 2,
  legend.labs = c("low", "high")
) 

pdf_fn = file.path("../results/Figure2D.pdf")
message(pdf_fn)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file=pdf_fn,width=5,height=5)
p
dev.off()
message("Done.")

print(sprintf('[Figure 2-D]  ACTA2-Low: %3d vs ACTA2-High: %3d', sum(data_f_MSS$ACTA2_bin==1), sum(data_f_MSS$ACTA2_bin==2)))
print(fit_MSS) # print some information including the median survival times
