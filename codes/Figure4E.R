source('./lib_project_lite.r')

#input parameters
args <- data.table(ici_spgc_geomx_rds="../data/ici_spgc_geomx.rds",
									 goi='ACTA2',
									 outd='../results/mf04e')

#create output directory
create_dir_if_not_exist(args$outd)

#####################
stopifnot(file.exists(args$ici_spgc_geomx_rds))
geomx.lst = readRDS(args$ici_spgc_geomx_rds)
geomx.lst$smeta$sbarcode = rownames(geomx.lst$smeta)

##############
resp_levels <- c("R","NR")
goi_exprj.dt <- geomx.lst$cnt %>%
  reshape2::melt(varnames=c("gene","sbarcode"),value.name = "q3") %>%
  dplyr::filter(gene==args$goi) %>%
  left_join(y=geomx.lst$smeta,by="sbarcode") %>%
  dplyr::filter(tumor_status=="Tumor") %>%
  dplyr::rename(msi_status=MSI) %>%
  mutate(response:=factor(response,levels=resp_levels)) %>%
  mutate(q3_log:=log(q3)) %>%
  as.data.table()
group.colors = c("#00BFC4","#F8766D")
names(group.colors) = resp_levels
norm_method = "Q3"
p = goi_exprj.dt %>%
  group_by(SegmentLabel,scanlabel)%>%
  summarize(avg_q3=unique(median(q3_log)), response=unique(response)) %>%
  mutate(response:=factor(response,levels=resp_levels)) %>%
  ggplot(aes(x=response,fill=response,y=avg_q3)) +
  scale_fill_manual(values=group.colors) +
  geom_boxplot(outlier.shape=NA,alpha=0.5) +
  geom_point(aes(shape=response),size=1.5, position = position_jitterdodge(jitter.width = 0.3)) +
  scale_shape_manual(values=c(1,1))+
  # ggrepel::geom_text_repel(aes(label=scanlabel),size=3,min.segment.length = 0.5, box.padding=0.5)+
  stat_compare_means(method="wilcox.test",
                     hide.ns=T,
                     vjust=0.8,
                     label = "p.format",
                     comparisons = list(resp_levels)) +
  facet_grid(cols=vars(SegmentLabel)) +
  labs(y="ACTA2 Expression Level [log(Q3)]")+
  theme_classic2()

pdf_fn = file.path(args$outd,"mf04e_goi_by_resp_SPark2022.pdf")
message(pdf_fn)
pdf(file=pdf_fn,width=6,height=4)
plot(p)
dev.off()
message("Done.")