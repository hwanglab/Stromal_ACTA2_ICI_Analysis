graphics.off()

source('./lib_project_lite.r')

args=data.table(goi = "ACTA2",
	exptag = "kwon2021_ACTA2",
	# seu_rds = "../data/kwon2021_GC_MSIH-all_samples-filtered.rds",
	acta2_gexpr_rds="../data/kwon2021_GC_MSIH-ACTA2.rds",
	outd="../results/mf04a")

create_dir_if_not_exist(args$outd)

if (F) {
	seu = readRDS(args$seu_rds)
	
	gexpr <- GetAssayData(seu,assay = "RNA", slot="data")[args$goi,]
	
	seu <- AddMetaData(seu,
										 metadata = gexpr,
										 col.name = 'expr')
	
	goi_expr.dt <- seu@meta.data %>%
		as.data.table(keep.rownames = TRUE)
	
	saveRDS(goi_expr.dt, file=args$acta2_gexpr_rds)
} else {
	stopifnot(file.exists(args$acta2_gexpr_rds))
	goi_expr.dt = readRDS(args$acta2_gexpr_rds)
}

goi_expr.dt[,TimePoint:=factor(TimePoint,levels=c("Pre","Post"))]

goi_expr.dt[,patid:=gsub("-[BF]$","",orig.ident)]

# goi_exprj.dt = goi_expr.dt[annot_2nd=="Stromal",]
goi_exprj.dt = goi_expr.dt[TimePoint=="Pre",]

goi_exprj.dt[annot_2nd=="Stromal",annot_group:="Stromal"]
goi_exprj.dt[annot_2nd=="Epithelial",annot_group:="Epithelial"]
goi_exprj.dt[annot_2nd %in% c("T","NK","B","Myeloid"),annot_group:="Immune"]

a=goi_exprj.dt[,.(mean=mean(expr),sd=sd(expr)),by=c("annot_group","patid")]

tsv_fn = file.path(args$outd,sprintf("%s_lognorm_mean.tsv",args$exptag))
message(tsv_fn)
fwrite(a,file=tsv_fn,sep="\t")

##################
setnames(goi_exprj.dt,'Response','response')
goi_exprj.dt[response=="Responder",response:="R"]
goi_exprj.dt[response=="NonResponder",response:="NR"]

goi_exprj.dt[,response:=factor(response,levels=c("R","NR"))]
upaid = goi_exprj.dt[,sort(unique(patid))]
goi_exprj.dt[,patid:=factor(patid,levels=upaid)]

my_comps = list(c("EP-75","EP-72"),c("EP-76","EP-72"),c("EP-77","EP-72"),c("EP-78","EP-72"))

p = ggplot(goi_exprj.dt[annot_group=="Stromal",],aes(x=patid,y=expr,fill=response)) +
	geom_violin(alpha=0.2)+
	geom_quasirandom(alpha = 0.5, dodge.width = 0.2, size=.5) +
	stat_summary(fun = mean, na.rm = TRUE,
							 geom = "point", color = "red",
							 size = 3, shape = "diamond") +
	stat_compare_means(method="wilcox.test",
										 hide.ns=T,
										 vjust=0.8,
										 label = "p.signif",
										 comparisons = my_comps) +
	labs(title="ACTA2 in pre-treatment samples\nKwon et.al., Wilcoxon Rank-sum test", y="Log Normalized ACTA2 Expression Level") +
	theme_bw() +
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf_fn = file.path(args$outd,sprintf("%s_violin_by_patid.pdf",args$goi))
message(pdf_fn)
pdf(file=pdf_fn,width=10,height=10)
plot(p)
dev.off()
message("Done.")