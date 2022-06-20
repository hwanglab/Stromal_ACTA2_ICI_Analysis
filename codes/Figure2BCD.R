library("survminer");
library("survival")

#-------------------------------------------------------------------------------
# Figure 2. ACTA2 expression was independently associated with overall survival 
# in gastric cancer patients. 

# The following codes reproduce the Figure 2-(B, C and D) which were generated from the pooled data 

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

#- load data
data_f <- read.table("../data/pooled_data_MSI_ACTA2.csv", sep=",", header=TRUE);
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

fileConn<-file("../results/Figure2B_additional.txt", open = 'a')
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
#- Figure 2 - (C & D): ) Kaplan-Meier plots

mr_thres = quantile(data_f$ACTA2[!is.na(data_f$ACTA2)], c(0.25))

x = rep(NA, rep=nrow(data_f));
x[data_f$ACTA2 <= mr_thres] = 0
x[data_f$ACTA2  > mr_thres] = 1

data_f_tmp = cbind(data_f, MSI_ACTA2 = interaction(data_f$MSI, factor(x)))
data_f_sub = data_f_tmp[!(is.na(data_f_tmp$MSI_ACTA2) | is.na(data_f_tmp$surv_time_m)),]

# microsatellite instability-high patients 
data_f_MSI = data_f_sub[data_f_sub$MSI==1, ]
surv_MSI = data.frame(data_f_MSI$surv_time_m,  data_f_MSI$vital_sign_1_d, data_f_MSI$MSI_ACTA2)
colnames(surv_MSI) <- c("time", "status", "gInfo")

# microsatellite stable patients
data_f_MSS = data_f_sub[data_f_sub$MSI==0, ]
surv_MSS = data.frame(data_f_MSS$surv_time_m,  data_f_MSS$vital_sign_1_d, data_f_MSS$MSI_ACTA2)
colnames(surv_MSS) <- c("time", "status", "gInfo")

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

print(fit_MSS) # print some information including the median survival times