# You can use install.packages("package name") to install the packages.

library(data.table)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(plyr)
library(dplyr)
library(reshape2)

create_dir_if_not_exist <- function(newd,active=TRUE) {
	if (!file.exists(newd)){
		if (active==T) {
			message(sprintf("creating %s ...",newd))
			dir.create(newd,showWarnings = F,recursive = T)
		}
	}
	message(newd)
	newd
}