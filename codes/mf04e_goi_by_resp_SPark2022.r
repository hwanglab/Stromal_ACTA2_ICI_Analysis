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
	melt(varnames=c("gene","sbarcode"),value.name = "q3") %>%
	dplyr::filter(gene==args$goi) %>%
	left_join(y=geomx.lst$smeta,by="sbarcode") %>%
	dplyr::filter(tumor_status=="Tumor") %>%
	rename(msi_status=MSI) %>%
	mutate(response:=factor(response,levels=resp_levels)) %>%
	as.data.table()

group.colors = c("#00bfc4","#f8766d")
names(group.colors) = resp_levels
norm_method = "Q3"

p = ggplot(goi_exprj.dt,aes(x=response,y=q3,fill=response)) +
	geom_violin(alpha=0.2)+
	geom_quasirandom(alpha = 0.5, dodge.width = 0.2, size=.5) +
	scale_fill_manual(values=group.colors) +
	stat_summary(fun = mean, na.rm = TRUE, 
							 geom = "point", color = "red", 
							 size = 3, shape = "diamond") +
	stat_compare_means(method="wilcox.test",
										 hide.ns=T,
										 vjust=0.8,
										 label = "p.format",
										 # label = "p.signif",
										 comparisons = list(resp_levels)) +
	facet_grid(cols=vars(SegmentLabel)) +
	# theme(axis.text.x = element_blank(),legend.position="right") +
	labs(title=sprintf("%s in %s (GeoMx DSP) Wilcoxon rank-sum test\nSt.Mary & Yonsei patient tumor ROIs",args$goi,norm_method), y="gene expression [Q3]") + 
	theme_bw() + 
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf_fn = file.path(args$outd,"mf04e_goi_by_resp_SPark2022.pdf")
message(pdf_fn)
pdf(file=pdf_fn,width=6,height=4)
plot(p)
dev.off()
message("Done.")