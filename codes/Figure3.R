library("ggplot2")
library(plyr)

#-------------------------------------------------------------------------------
# Figure 3. ACTA2 expression predicted response to immune checkpoint inhibitors. 
# Abbreviations: CR, complete response; PR, partial response; 
#                SD, stable disease; PD, progressive disease.
#-------------------------------------------------------------------------------

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
tst <- chisq.test(mat_counts)

# generate a bar plot
mat_res <- matrix(0, nrow = 4, ncol = 3);

mat_res[1:2, 1] = c("Low")
mat_res[3:4, 1] = c("High")

mat_res[c(1,3), 2] = c("SD/PD")
mat_res[c(2,4), 2] = c("CR/PR")

# mat_res[, 3] = Reshape(mat_counts, 4, 1)
mat_res[, 3] = c(mat_counts[1,1], mat_counts[2,1], mat_counts[1,2], mat_counts[2,2])

df_res = data.frame(mat_res);
colnames(df_res) = c("ACTA2_bin", "response", "patient_count")

df_res$patient_count <- as.numeric(df_res$patient_count)
df_res$ACTA2_bin <- factor(df_res$ACTA2_bin, levels=c("Low", "High"))

df_res <- ddply(df_res, .(ACTA2_bin), 
              transform, pos = cumsum(patient_count) - (0.5 * patient_count))

p <- ggplot(df_res, aes(x = ACTA2_bin, y = patient_count)) +
  geom_bar(aes(fill = response), stat="identity") +
  xlab("ACTA2 expression") + ylab("Number of patients") +
  geom_text(aes(label = patient_count, y = pos), size = 3) + 
  geom_text(aes(x=1.5, y=83, label=sprintf('p-value=%.3f', tst$p.value)), size = 3) +
  geom_text(aes(x=1.5, y=80, label=sprintf('|-------------------------|', tst$p.value)), size = 3) 

dir.create("../results/", showWarnings = FALSE)

pdf_fn = file.path("../results/Figure3.pdf")
message(pdf_fn)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file=pdf_fn,width=5,height=5)
p
dev.off()
message("Done.")
