# install.packages('pracma')
library(pracma)

#-------------------------------------------------------------------------------
# Figure 3. ACTA2 expression predicted response to immune checkpoint inhibitors. 
# Abbreviations: CR, complete response; PR, partial response; 
#                SD, stable disease; PD, progressive disease.
#
# The data is from 
#-------------------------------------------------------------------------------

# the current working directory 
setwd("D:\\Research\\codes\\ACTA2_Analysis\\paper_codes\\codes\\")

#--- data loading
my_data_immun <- read.table(file = '../data/data_figure_3.txt')


# the cutoff value for the ACTA2 expression values (high vs low): 25% quantile
ACTA2_thres = quantile(my_data_immun$ACTA2, c(0.25))

x1 = rep(NA, rep=length(my_data_immun$ACTA2))
x1[my_data_immun$ACTA2 <= ACTA2_thres] = 1;
x1[my_data_immun$ACTA2 >  ACTA2_thres] = 2;

my_data_immun$ACTA2_bin = factor(x1);

# Count patients that belong to each combination (ACTA2 high or low X patients' treatment responses)  
mat_counts =  matrix(0, nrow = 2, ncol = 2);
mat_counts[1, 1] = sum((my_data_immun$ACTA2_bin == 1) & (my_data_immun$res_b1 == 0))
mat_counts[2, 1] = sum((my_data_immun$ACTA2_bin == 1) & (my_data_immun$res_b1 == 1))
mat_counts[1, 2] = sum((my_data_immun$ACTA2_bin == 2) & (my_data_immun$res_b1 == 0))
mat_counts[2, 2] = sum((my_data_immun$ACTA2_bin == 2) & (my_data_immun$res_b1 == 1))

print(mat_counts)

# Fisher's exact test
df_p_val <- data.frame(
  group1 = "Low",
  group2 = "High",
  label = sprintf('p=%.5f', fisher.test(mat_counts)$p.value),
  x = 0.5,
  y.position = sum(mat_counts[,2]) + 2
)


# generate a bar plot
mat_res <- matrix(0, nrow = 4, ncol = 3);

mat_res[1:2, 1] = c("Low")
mat_res[3:4, 1] = c("High")

mat_res[c(1,3), 2] = c("SD/PD")
mat_res[c(2,4), 2] = c("CR/PR")

mat_res[, 3] = Reshape(mat_counts, 4, 1)

df_res = data.frame(mat_res);
colnames(df_res) = c("ACTA2_bin", "response", "patient_count")

df_res$patient_count <- as.numeric(df_res$patient_count)
df_res$ACTA2_bin <- factor(df_res$ACTA2_bin, levels=c("Low", "High"))

library(plyr)
df_res <- ddply(df_res, .(ACTA2_bin), 
              transform, pos = cumsum(patient_count) - (0.5 * patient_count))

ggplot(df_res, aes(x = ACTA2_bin, y = patient_count)) +
  geom_bar(aes(fill = response), stat="identity") +
  xlab("ACTA2 expression") + ylab("Number of patients") +
  geom_text(aes(label = patient_count, y = pos), size = 3) +
  add_pvalue(df_p_val) # bracket

