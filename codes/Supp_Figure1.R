library(ggplot2)
library(ggpubr)

#-------------------------------------------------------------------------------
# Supplementary 1. ACTA2 expression was elevated in the molecular subtype with 
# the worse prognosis (Group 4). The bar in the middle of the box represents the median. 
# Student's t-test was used to compare groups. Abbreviation: ****, p <0.0001.
#-------------------------------------------------------------------------------

# ACTA2
df_input <- read.table("../data/Yonsei_567_Clinical_ACTA2_Group.txt", header = TRUE, sep = "\t");
df_input$group_id = factor(df_input$group_id, levels=c(1,2,3,4))
df_input_sub = df_input[,c(1, 14, 15)]

my_comparisons <- list(c("4", "1"), c("4", "2"), c("4", "3"))

p <- ggplot(df_input_sub, aes(x=group_id, y=ACTA2, color=group_id)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.85), size=1.8, alpha=0.7, 
             aes(shape=group_id, group=group_id)) +
  theme(legend.position="right",  panel.background = element_rect(fill='transparent')) +
  labs(title="ACTA2 values in Group 4 vs others", 
       x="group", y = "ACTA2 Expression") +
  stat_compare_means(method='t.test', vjust=0.8, label = "p.signif", comparisons=my_comparisons) 

dir.create("../results/", showWarnings = FALSE)

pdf_fn = file.path("../results/Supp_Figure1.pdf")
message(pdf_fn)
pdf.options(reset = TRUE, onefile = FALSE)
pdf(file=pdf_fn,width=5,height=5)
p
dev.off()
message("Done.")
