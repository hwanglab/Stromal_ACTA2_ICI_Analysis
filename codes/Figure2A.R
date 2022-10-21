library("survminer")
library("survival")

#-------------------------------------------------------------------------------
# Figure 2: ACTA2 expression was independently associated with overall survival 
# in gastric cancer patients. 

# The following codes reproduce the Figure 2-A which was generated using the Yonsei cohort 

# (A) Multivariable analysis of the Yonsei cohort (samples from 567 gastric adenocarcinoma patients)
#-------------------------------------------------------------------------------

# the current working directory 
# setwd("D:\\Research\\codes\\00_distributions\\Stromal_ACTA2_ICI_Analysis\\codes")

#- load data
my_data <- read.table("../data/Yonsei_567_OS_ACTA2_MSI_EBV_Subtypes.txt", sep="\t", header=TRUE);

#- Stage
Stage_CODE = factor(x = my_data$Stage)

#- Sex
Sex_CODE = factor(my_data$Sex, levels = c("Female", "Male"));

# Age
x1 = rep(NA, rep=nrow(my_data));
x1[my_data$Age == "?60"] = '60below';
x1[my_data$Age == ">60"] = '60older';

Age_CODE = factor(x1, levels = c('60below', '60older'));

#- Survival time (months)
Surv_Months = as.numeric(my_data$'Overall.survival..months.')

#- Vital sign
Surv_Sign = as.numeric(my_data$'Vital.status..0..Alive.1..Dead.')

#- ACTA2 expression values
data_f = data.frame(Stage_CODE, Sex_CODE, Age_CODE, 
                    Surv_Months, Surv_Sign, my_data$ACTA2)
colnames(data_f) <- c("Stage", "Sex_", "Age", "surv_time_m", "vital_sign_1_d", "ACTA2");


#-------------------------------------------------------------------------------
#- Figure 2-A: Multivariable Cox proportional hazards model 
fit_mv <- coxph(formula = Surv(surv_time_m, vital_sign_1_d) 
              ~ Age + Sex_ + Stage + ACTA2, data = data_f)

# result summary
summary(fit_mv)

beta <- exp(coef(fit_mv))
CI   <- round(exp(confint(fit_mv)), 5)

coeffs <- coef(summary(fit_mv))
p    <- as.matrix(coeffs[,5])

dir.create("../results/", showWarnings = FALSE)

fileConn<-file("../results/Figure2A_additional.txt", open = 'w')
for (mn_i in 1:length(beta)){
  str_line = sprintf("%10s:: %.2f (%.2f, %.2f), pval = %.5f", 
                     rownames(coeffs)[mn_i], beta[mn_i], CI[mn_i,1], CI[mn_i,2], p[mn_i])
  writeLines(str_line, fileConn)
}
close(fileConn)

# Visualize the results from the Cox model using forest plot
pdf_fn = file.path("../results/Figure2A.pdf")
message(pdf_fn)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file=pdf_fn,width=5,height=5)
ggforest(fit_mv, data = data_f)
dev.off()
message("Done.")
