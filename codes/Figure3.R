library("ggplot2")
library("plyr")

#-------------------------------------------------------------------------------
# Figure 3. ACTA2 expression predicted response to immune checkpoint inhibitors. 
# Abbreviations: CR, complete response; PR, partial response; 
#                SD, stable disease; PD, progressive disease.
#-------------------------------------------------------------------------------

#--- data loading
my_data_immun <- read.table(file = '../data/data_figure_3_ACTA2_MSI_EBV.txt', header = TRUE, sep= "\t")

mstr_save_path = "../results/"
dir.create(mstr_save_path, showWarnings = FALSE)

Response_boxplot <- function (data_in, str_filenname) {
  # the cutoff value for the ACTA2 expression values (low vs high): 25% quantile
  str_thres = "25_quantile"
  ACTA2_thres = quantile(data_in$ACTA2, c(0.25))
  
  print(sprintf('dataset=%s, threshold=%s (%.3f)', str_filenname, str_thres, ACTA2_thres))
  
  x5 = rep(NA, rep=length(data_in$ACTA2))
  x5[data_in$ACTA2 <= ACTA2_thres] = 1;
  x5[data_in$ACTA2 >  ACTA2_thres] = 2;
  
  data_display = data.frame(data_in, ACTA2_bin=factor(x5));
  
  # Count patients that belong to each combination (ACTA2 high or low X patients' treatment responses)  
  mat_counts =  matrix(0, nrow = 2, ncol = 2);
  mat_counts[1, 1] = sum((data_display$ACTA2_bin == 1) & (data_display$res_b1 == 0))
  mat_counts[2, 1] = sum((data_display$ACTA2_bin == 1) & (data_display$res_b1 == 1))
  mat_counts[1, 2] = sum((data_display$ACTA2_bin == 2) & (data_display$res_b1 == 0))
  mat_counts[2, 2] = sum((data_display$ACTA2_bin == 2) & (data_display$res_b1 == 1))
  
  tst <- fisher.test(mat_counts)
  
  df_p_val <- data.frame(
    group1 = "Low",
    group2 = "High",
    label = sprintf('p=%.5f', tst$p.value),
    x = 0.5,
    y.position = sum(mat_counts[,2]) + 2
  )
  
  # generate a bar plot
  mat_res <- matrix(0, nrow = 4, ncol = 3);
  
  mat_res[1:2, 1] = c("Low")
  mat_res[3:4, 1] = c("High")
  
  mat_res[c(1,3), 2] = c("SD/PD")
  mat_res[c(2,4), 2] = c("CR/PR")
  
  # mat_res[, 3] = Reshape(mat_counts, 4, 1);
  mat_res[, 3] = c(mat_counts[1,1], mat_counts[2,1], mat_counts[1,2], mat_counts[2,2])
  
  df_res = data.frame(mat_res);
  colnames(df_res) = c("ACTA2_bin", "response", "pat_count")
  
  df_res$pat_count <- as.numeric(df_res$pat_count)
  df_res$ACTA2_bin <- factor(df_res$ACTA2_bin, levels=c("Low", "High"))
  
  df_res <- ddply(df_res, .(ACTA2_bin), 
                  transform, pos = cumsum(pat_count) - (0.5 * pat_count))
  
  p <- ggplot(df_res, aes(x = ACTA2_bin, y = pat_count)) +
    geom_bar(aes(fill = response), stat="identity") +
    xlab("ACTA2 expression") + ylab("Number of patients") +
    geom_text(aes(label = pat_count, y = pos), size = 3) +
    geom_text(aes(x=1.5, y=sum(mat_counts[, 2])+3, label=sprintf('p-value=%.3f', tst$p.value)), size = 3) +
    geom_text(aes(x=1.5, y=sum(mat_counts[, 2])+2, label='|------------------------------------------------|'), size = 3) 
    # + add_pvalue(df_p_val) # bracket
  
  str_full_path = paste0(mstr_save_path, 'Figure3', str_filenname, '.pdf')
  message(str_full_path)
  p
  ggsave(str_full_path, p, width = 10, height = 8, dpi = 150, units = "in")
}

Response_boxplot(my_data_immun, 'A-All_patients')
Response_boxplot(my_data_immun[(my_data_immun$MSI==1) & !is.na(my_data_immun$MSI),], 'B-MSI_H')
Response_boxplot(my_data_immun[(my_data_immun$MSI==0) & !is.na(my_data_immun$MSI),], 'C-MSS')
