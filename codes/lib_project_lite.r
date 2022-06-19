library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggbeeswarm)

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